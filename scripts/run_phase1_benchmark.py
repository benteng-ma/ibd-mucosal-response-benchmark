#!/usr/bin/env python3
"""Run the frozen published-program benchmark without model fitting."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import random
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXPR = ROOT / "data" / "processed" / "phase1_v1"
MANIFEST = ROOT / "metadata" / "phase1"
OUT = ROOT / "results" / "performance" / "phase1_v1"


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(rows[0]) if rows else []
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def percentile(values: list[float], probability: float) -> float:
    values = sorted(values)
    if not values:
        return math.nan
    index = (len(values) - 1) * probability
    lo = math.floor(index)
    hi = math.ceil(index)
    if lo == hi:
        return values[lo]
    return values[lo] * (hi - index) + values[hi] * (index - lo)


def average_ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i + 1
        while j < len(order) and values[order[j]] == values[order[i]]:
            j += 1
        rank = ((i + 1) + j) / 2
        for k in range(i, j):
            ranks[order[k]] = rank
        i = j
    return ranks


def auroc(labels: list[int], scores: list[float]) -> float:
    n_pos = sum(labels)
    n_neg = len(labels) - n_pos
    if n_pos == 0 or n_neg == 0:
        return math.nan
    ranks = average_ranks(scores)
    rank_sum = sum(rank for rank, label in zip(ranks, labels) if label == 1)
    return (rank_sum - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)


def average_precision(labels: list[int], scores: list[float]) -> float:
    n_pos = sum(labels)
    if n_pos == 0:
        return math.nan
    groups: dict[float, list[int]] = defaultdict(list)
    for label, score in zip(labels, scores):
        groups[score].append(label)
    tp = 0
    fp = 0
    ap = 0.0
    previous_recall = 0.0
    for score in sorted(groups, reverse=True):
        group = groups[score]
        tp += sum(group)
        fp += len(group) - sum(group)
        recall = tp / n_pos
        precision = tp / (tp + fp)
        ap += (recall - previous_recall) * precision
        previous_recall = recall
    return ap


def bootstrap_auc(labels: list[int], scores: list[float], replicates: int, seed: int) -> tuple[float, float]:
    rng = random.Random(seed)
    positives = [(y, s) for y, s in zip(labels, scores) if y == 1]
    negatives = [(y, s) for y, s in zip(labels, scores) if y == 0]
    estimates: list[float] = []
    for _ in range(replicates):
        sample = [rng.choice(positives) for _ in positives] + [rng.choice(negatives) for _ in negatives]
        estimates.append(auroc([x[0] for x in sample], [x[1] for x in sample]))
    return percentile(estimates, 0.025), percentile(estimates, 0.975)


def permutation_p(labels: list[int], scores: list[float], replicates: int, seed: int) -> float:
    observed = abs(auroc(labels, scores) - 0.5)
    rng = random.Random(seed)
    permuted = labels[:]
    exceed = 0
    for _ in range(replicates):
        rng.shuffle(permuted)
        if abs(auroc(permuted, scores) - 0.5) >= observed - 1e-15:
            exceed += 1
    return (exceed + 1) / (replicates + 1)


def stable_seed(signature_id: str, task_id: str) -> int:
    token = hashlib.sha256(f"20260831::{signature_id}::{task_id}".encode()).hexdigest()[:8]
    return int(token, 16)


def holm_adjust(rows: list[dict]) -> None:
    eligible = [(i, float(r["permutation_p_value"])) for i, r in enumerate(rows) if r["permutation_p_value"]]
    eligible.sort(key=lambda x: x[1])
    running = 0.0
    m = len(eligible)
    for rank, (index, p_value) in enumerate(eligible):
        adjusted = min(1.0, p_value * (m - rank))
        running = max(running, adjusted)
        rows[index]["holm_adjusted_p_value"] = f"{running:.6g}"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    signature_inventory = {
        r["signature_id"]: r
        for r in read_csv(ROOT / "literature" / "signature_inventory.csv")
        if r["main_grid_eligible"].upper() == "TRUE"
    }
    signature_genes: dict[str, list[dict]] = defaultdict(list)
    for row in read_csv(ROOT / "literature" / "signature_gene_table.csv"):
        if row["signature_id"] in signature_inventory:
            signature_genes[row["signature_id"]].append(row)

    summaries = read_csv(MANIFEST / "phase1_task_summary.csv")
    task_ids = {r["task_id"] for r in summaries if r["status"] == "EXECUTABLE"}
    evaluation_rows = [
        r for r in read_csv(ROOT / "metadata" / "evaluation_grid.csv")
        if r["task_id"] in task_ids and r["transportability_eligible"].upper() == "TRUE"
    ]

    task_manifests = {
        task_id: read_csv(MANIFEST / f"{task_id}_analysis_manifest.csv")
        for task_id in task_ids
    }
    dataset_samples: dict[str, set[str]] = defaultdict(set)
    for rows in task_manifests.values():
        for row in rows:
            dataset_samples[row["dataset_id"]].add(row["sample_id"])

    standardized: dict[str, dict[str, dict[str, float | None]]] = {}
    for dataset_id, sample_ids in dataset_samples.items():
        expression = {r["sample_id"]: r for r in read_csv(EXPR / f"{dataset_id}_signature_genes.csv") if r["sample_id"] in sample_ids}
        genes = [g for g in next(iter(expression.values())) if g != "sample_id"]
        standardized[dataset_id] = {sample_id: {} for sample_id in sample_ids}
        for gene in genes:
            values = [float(expression[s][gene]) for s in sample_ids if expression[s].get(gene, "") != ""]
            if len(values) != len(sample_ids):
                for sample_id in sample_ids:
                    standardized[dataset_id][sample_id][gene] = None
                continue
            mean = sum(values) / len(values)
            sd = math.sqrt(sum((x - mean) ** 2 for x in values) / len(values))
            for sample_id in sample_ids:
                standardized[dataset_id][sample_id][gene] = (float(expression[sample_id][gene]) - mean) / sd if sd > 0 else None

    result_rows: list[dict] = []
    score_rows: list[dict] = []
    for grid in evaluation_rows:
        signature_id = grid["signature_id"]
        task_id = grid["task_id"]
        dataset_id = grid["dataset_id"]
        genes = signature_genes[signature_id]
        available = [g for g in genes if all(standardized[dataset_id][s].get(g["gene_symbol"]) is not None for s in dataset_samples[dataset_id])]
        coverage = len(available) / len(genes)
        status = "EXECUTED" if coverage >= 0.8 else "SKIPPED_ACTUAL_COVERAGE_LT_0.80"
        base = {
            "signature_id": signature_id,
            "task_id": task_id,
            "dataset_id": dataset_id,
            "trial_id": grid["trial_id"],
            "evaluation_role": grid["evaluation_role"],
            "lineage_class": grid["lineage_class"],
            "same_therapy": grid["same_therapy"],
            "same_endpoint": grid["same_endpoint"],
            "same_disease": grid["same_disease"],
            "n_subjects": len(task_manifests[task_id]),
            "n_responder": sum(int(r["outcome_label"]) for r in task_manifests[task_id]),
            "n_nonresponder": sum(1 - int(r["outcome_label"]) for r in task_manifests[task_id]),
            "n_signature_genes": len(genes),
            "n_genes_used": len(available),
            "actual_gene_coverage": f"{coverage:.3f}",
            "execution_status": status,
            "score_definition": "published-coefficient L1-weighted z sum" if signature_inventory[signature_id]["coefficients_available"].upper() == "TRUE" else "signed directional z mean",
            "auroc": "",
            "auroc_ci_lower": "",
            "auroc_ci_upper": "",
            "average_precision": "",
            "outcome_prevalence": "",
            "prevalence_normalized_ap": "",
            "rank_biserial_separation": "",
            "permutation_p_value": "",
            "holm_adjusted_p_value": "",
            "transport_classification": "NOT_EVALUATED",
            "calibration_status": "NOT_APPLICABLE_UNCALIBRATED_PROGRAM_SCORE",
        }
        if status != "EXECUTED":
            result_rows.append(base)
            continue

        weights: dict[str, float] = {}
        if signature_inventory[signature_id]["coefficients_available"].upper() == "TRUE":
            raw_weights = {g["gene_symbol"]: float(g["coefficient"]) for g in available}
            denominator = sum(abs(x) for x in raw_weights.values())
            weights = {gene: weight / denominator for gene, weight in raw_weights.items()}
        else:
            raw_weights = {g["gene_symbol"]: 1.0 if g["direction"] == "R_HIGH" else -1.0 for g in available}
            weights = {gene: weight / len(raw_weights) for gene, weight in raw_weights.items()}

        labels: list[int] = []
        scores: list[float] = []
        for sample in task_manifests[task_id]:
            score = sum(weights[g] * float(standardized[dataset_id][sample["sample_id"]][g]) for g in weights)
            label = int(sample["outcome_label"])
            labels.append(label)
            scores.append(score)
            score_rows.append({
                "signature_id": signature_id,
                "task_id": task_id,
                "dataset_id": dataset_id,
                "trial_id": grid["trial_id"],
                "sample_id": sample["sample_id"],
                "subject_id": sample["subject_id"],
                "outcome_label": label,
                "score": f"{score:.10g}",
                "score_orientation": "higher predicts responder",
                "n_genes_used": len(available),
            })

        auc = auroc(labels, scores)
        ap = average_precision(labels, scores)
        prevalence = sum(labels) / len(labels)
        replicates = 2000 if grid["evaluation_role"] == "SAME_THERAPY_SAME_ENDPOINT_EXTERNAL" else 1000
        seed = stable_seed(signature_id, task_id)
        ci_low, ci_high = bootstrap_auc(labels, scores, replicates, seed)
        classification = "SUPPORTED_DIRECTION" if ci_low > 0.5 else "REVERSED_DIRECTION" if ci_high < 0.5 else "UNCERTAIN_OR_NONTRANSPORTABLE"
        base.update({
            "auroc": f"{auc:.6f}",
            "auroc_ci_lower": f"{ci_low:.6f}",
            "auroc_ci_upper": f"{ci_high:.6f}",
            "average_precision": f"{ap:.6f}",
            "outcome_prevalence": f"{prevalence:.6f}",
            "prevalence_normalized_ap": f"{ap / prevalence:.6f}",
            "rank_biserial_separation": f"{2 * auc - 1:.6f}",
            "transport_classification": classification,
        })
        if grid["evaluation_role"] == "SAME_THERAPY_SAME_ENDPOINT_EXTERNAL":
            base["permutation_p_value"] = f"{permutation_p(labels, scores, 10000, seed + 1):.6g}"
        result_rows.append(base)

    holm_adjust(result_rows)
    fields = list(result_rows[0])
    write_csv(OUT / "study_specific_evaluations.csv", result_rows, fields)
    write_csv(OUT / "subject_level_scores.csv", score_rows, list(score_rows[0]))
    write_csv(OUT / "strict_external_evaluations.csv", [r for r in result_rows if r["evaluation_role"] == "SAME_THERAPY_SAME_ENDPOINT_EXTERNAL"], fields)
    write_csv(OUT / "transportability_matrix.csv", [r for r in result_rows if r["execution_status"] == "EXECUTED"], fields)
    write_csv(
        OUT / "negative_transportability_matrix.csv",
        [r for r in result_rows if r["transport_classification"] in {"UNCERTAIN_OR_NONTRANSPORTABLE", "REVERSED_DIRECTION"}],
        fields,
    )

    summary = {
        "version": "phase1_v1",
        "n_frozen_executable_grid_cells": len(evaluation_rows),
        "n_executed_cells": sum(r["execution_status"] == "EXECUTED" for r in result_rows),
        "n_skipped_actual_coverage": sum(r["execution_status"] != "EXECUTED" for r in result_rows),
        "n_strict_external_cells": sum(r["evaluation_role"] == "SAME_THERAPY_SAME_ENDPOINT_EXTERNAL" and r["execution_status"] == "EXECUTED" for r in result_rows),
        "n_cross_therapy_stress_cells": sum(r["evaluation_role"] == "CROSS_THERAPY_TRANSPORT_STRESS_TEST" and r["execution_status"] == "EXECUTED" for r in result_rows),
        "transport_classifications": dict(sorted({k: sum(r["transport_classification"] == k for r in result_rows) for k in {r["transport_classification"] for r in result_rows}}.items())),
        "model_fitting": False,
        "feature_selection": False,
        "cross_cohort_merging": False,
        "calibrated_probability_claim": False,
    }
    (OUT / "benchmark_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
