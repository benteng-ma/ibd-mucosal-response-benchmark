#!/usr/bin/env python3
"""Run the frozen Phase 3 exploratory robustness analyses.

The script never modifies Phase 1 outputs. It recomputes the frozen score only as an
integrity check, then writes exploratory diagnostics to results/phase3_robustness/v1.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import math
import random
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import cohen_kappa_score, roc_auc_score, roc_curve


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "phase3_robustness" / "v1"
SEED = 20260831
N_BOOT = 2000
N_CONTEXT_BOOT = 5000
N_STRICT_PERM = 10000


def stable_seed(*tokens: str, offset: int = 0) -> int:
    token = "::".join([str(SEED), *tokens])
    value = int(hashlib.sha256(token.encode()).hexdigest()[:8], 16)
    return (value + offset) % (2**32 - 1)


def parse_bool(value: object) -> bool:
    return str(value).strip().upper() == "TRUE"


def percentile(values: list[float], probability: float) -> float:
    return float(np.quantile(np.asarray(values, dtype=float), probability))


def classify(ci_low: float, ci_high: float) -> str:
    if ci_low > 0.5:
        return "SUPPORTED_DIRECTION"
    if ci_high < 0.5:
        return "REVERSED_DIRECTION"
    return "UNCERTAIN_OR_NONTRANSPORTABLE"


def raw_program_weights(signature_id: str, genes: list[str], inventory: pd.DataFrame, gene_table: pd.DataFrame) -> dict[str, float]:
    row = inventory.loc[signature_id]
    subset = gene_table[(gene_table["signature_id"] == signature_id) & gene_table["gene_symbol"].isin(genes)]
    if parse_bool(row["coefficients_available"]):
        weights = {r.gene_symbol: float(r.coefficient) for r in subset.itertuples()}
    else:
        weights = {r.gene_symbol: (1.0 if r.direction == "R_HIGH" else -1.0) for r in subset.itertuples()}
    denominator = sum(abs(value) for value in weights.values())
    return {gene: value / denominator for gene, value in weights.items()}


def score_matrix(matrix: pd.DataFrame, weights: dict[str, float]) -> np.ndarray:
    genes = list(weights)
    vector = np.asarray([weights[gene] for gene in genes], dtype=float)
    return matrix[genes].to_numpy(dtype=float) @ vector


def build_scaled_matrices(expression: pd.DataFrame) -> dict[str, pd.DataFrame]:
    matrix = expression.set_index("sample_id").astype(float)
    means = matrix.mean(axis=0)
    population_sd = matrix.std(axis=0, ddof=0)
    if (population_sd <= 0).any():
        raise ValueError(f"Zero population SD for {list(population_sd[population_sd <= 0].index)}")
    z = (matrix - means) / population_sd

    medians = matrix.median(axis=0)
    mad = (matrix - medians).abs().median(axis=0) * 1.4826
    iqr_scale = (matrix.quantile(0.75) - matrix.quantile(0.25)) / 1.349
    robust_scale = mad.where(mad > 0, iqr_scale)
    robust_scale = robust_scale.where(robust_scale > 0, population_sd)
    if (robust_scale <= 0).any():
        raise ValueError(f"No positive robust dispersion for {list(robust_scale[robust_scale <= 0].index)}")
    robust = (matrix - medians) / robust_scale

    rank = matrix.rank(axis=0, method="average", pct=True) - 0.5
    return {"z": z, "mad": robust, "rank": rank}


def stratified_bootstrap_direction_probability(labels: np.ndarray, scores: np.ndarray, seed: int) -> tuple[float, float, float]:
    rng = np.random.default_rng(seed)
    positive = np.flatnonzero(labels == 1)
    negative = np.flatnonzero(labels == 0)
    estimates = np.empty(N_BOOT, dtype=float)
    for index in range(N_BOOT):
        sampled = np.concatenate(
            [rng.choice(positive, size=len(positive), replace=True), rng.choice(negative, size=len(negative), replace=True)]
        )
        estimates[index] = roc_auc_score(labels[sampled], scores[sampled])
    return float(np.mean(estimates > 0.5)), float(np.quantile(estimates, 0.025)), float(np.quantile(estimates, 0.975))


def jackknife_auc(labels: np.ndarray, scores: np.ndarray) -> np.ndarray:
    estimates = []
    for index in range(len(labels)):
        keep = np.arange(len(labels)) != index
        estimates.append(roc_auc_score(labels[keep], scores[keep]))
    return np.asarray(estimates, dtype=float)


def two_way_context_bootstrap(cells: pd.DataFrame) -> dict[str, float | int]:
    observed = float(spearmanr(cells["alignment_score"], cells["auroc"]).statistic)
    programs = sorted(cells["signature_id"].unique())
    tasks = sorted(cells["task_id"].unique())
    rng = np.random.default_rng(stable_seed("context", "two_way_bootstrap"))
    estimates: list[float] = []
    for _ in range(N_CONTEXT_BOOT):
        sampled_programs = rng.choice(programs, size=len(programs), replace=True)
        sampled_tasks = rng.choice(tasks, size=len(tasks), replace=True)
        program_counts = pd.Series(sampled_programs).value_counts()
        task_counts = pd.Series(sampled_tasks).value_counts()
        expanded_x: list[float] = []
        expanded_y: list[float] = []
        for row in cells.itertuples():
            multiplicity = int(program_counts.get(row.signature_id, 0) * task_counts.get(row.task_id, 0))
            if multiplicity:
                expanded_x.extend([row.alignment_score] * multiplicity)
                expanded_y.extend([row.auroc] * multiplicity)
        if len(set(expanded_x)) > 1 and len(set(expanded_y)) > 1:
            estimate = spearmanr(expanded_x, expanded_y).statistic
            if np.isfinite(estimate):
                estimates.append(float(estimate))
    return {
        "observed_spearman_rho": observed,
        "bootstrap_replicates_requested": N_CONTEXT_BOOT,
        "bootstrap_replicates_valid": len(estimates),
        "bootstrap_ci_lower": percentile(estimates, 0.025),
        "bootstrap_ci_upper": percentile(estimates, 0.975),
        "interpretation": "descriptive two-way program/task bootstrap; cells are dependent",
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    frozen = pd.read_csv(ROOT / "results" / "performance" / "phase1_v1" / "study_specific_evaluations.csv")
    frozen = frozen[frozen["execution_status"] == "EXECUTED"].copy()
    if len(frozen) != 46:
        raise AssertionError(f"Expected 46 executed cells, found {len(frozen)}")
    frozen["auroc"] = frozen["auroc"].astype(float)
    frozen["auroc_ci_lower"] = frozen["auroc_ci_lower"].astype(float)
    frozen["auroc_ci_upper"] = frozen["auroc_ci_upper"].astype(float)
    for column in ["same_therapy", "same_endpoint", "same_disease"]:
        frozen[column] = frozen[column].map(parse_bool)

    inventory = pd.read_csv(ROOT / "literature" / "signature_inventory.csv").set_index("signature_id")
    inventory = inventory[inventory["main_grid_eligible"].map(parse_bool)]
    gene_table = pd.read_csv(ROOT / "literature" / "signature_gene_table.csv")
    gene_table = gene_table[gene_table["signature_id"].isin(inventory.index)].copy()

    task_ids = sorted(frozen["task_id"].unique())
    task_manifests = {
        task_id: pd.read_csv(ROOT / "metadata" / "phase1" / f"{task_id}_analysis_manifest.csv") for task_id in task_ids
    }
    dataset_sample_ids: dict[str, set[str]] = {}
    for manifest in task_manifests.values():
        dataset_id = str(manifest["dataset_id"].iloc[0])
        dataset_sample_ids.setdefault(dataset_id, set()).update(manifest["sample_id"].astype(str))

    scaled_by_dataset: dict[str, dict[str, pd.DataFrame]] = {}
    for dataset_id, sample_ids in dataset_sample_ids.items():
        expression = pd.read_csv(ROOT / "data" / "processed" / "phase1_v1" / f"{dataset_id}_signature_genes.csv")
        expression["sample_id"] = expression["sample_id"].astype(str)
        expression = expression[expression["sample_id"].isin(sample_ids)].copy()
        if len(expression) != len(sample_ids):
            raise AssertionError(f"{dataset_id}: expected {len(sample_ids)} samples, found {len(expression)}")
        scaled_by_dataset[dataset_id] = build_scaled_matrices(expression)

    cell_rows: list[dict] = []
    loo_rows: list[dict] = []
    sign_rows: list[dict] = []
    bootstrap_rows: list[dict] = []

    for cell in frozen.sort_values(["signature_id", "task_id"]).itertuples(index=False):
        signature_id = cell.signature_id
        task_id = cell.task_id
        dataset_id = cell.dataset_id
        manifest = task_manifests[task_id].copy()
        manifest["sample_id"] = manifest["sample_id"].astype(str)
        manifest = manifest.set_index("sample_id")
        sample_order = list(manifest.index)
        labels = manifest.loc[sample_order, "outcome_label"].astype(int).to_numpy()
        available_genes = [
            gene
            for gene in gene_table.loc[gene_table["signature_id"] == signature_id, "gene_symbol"].tolist()
            if gene in scaled_by_dataset[dataset_id]["z"].columns
            and scaled_by_dataset[dataset_id]["z"].loc[sample_order, gene].notna().all()
        ]
        weights = raw_program_weights(signature_id, available_genes, inventory, gene_table)
        method_scores = {
            method: score_matrix(matrices.loc[sample_order], weights)
            for method, matrices in scaled_by_dataset[dataset_id].items()
        }
        method_auc = {method: float(roc_auc_score(labels, scores)) for method, scores in method_scores.items()}
        if abs(method_auc["z"] - float(cell.auroc)) > 5e-6:
            raise AssertionError(
                f"Frozen AUROC mismatch for {signature_id}/{task_id}: {method_auc['z']} versus {cell.auroc}"
            )

        probability, boot_low, boot_high = stratified_bootstrap_direction_probability(
            labels, method_scores["z"], stable_seed(signature_id, task_id, "direction_bootstrap")
        )
        bootstrap_rows.append(
            {
                "signature_id": signature_id,
                "task_id": task_id,
                "replicates": N_BOOT,
                "probability_auc_gt_0_5": probability,
                "phase3_bootstrap_ci_lower": boot_low,
                "phase3_bootstrap_ci_upper": boot_high,
            }
        )

        jackknife = jackknife_auc(labels, method_scores["z"])
        loo_auc_values: list[float] = []
        most_influential_gene = ""
        most_influential_delta = math.nan
        if len(available_genes) >= 2:
            for omitted_gene in available_genes:
                retained = [gene for gene in available_genes if gene != omitted_gene]
                loo_weights = raw_program_weights(signature_id, retained, inventory, gene_table)
                loo_score = score_matrix(scaled_by_dataset[dataset_id]["z"].loc[sample_order], loo_weights)
                loo_auc = float(roc_auc_score(labels, loo_score))
                loo_auc_values.append(loo_auc)
                delta = loo_auc - method_auc["z"]
                loo_rows.append(
                    {
                        "signature_id": signature_id,
                        "task_id": task_id,
                        "dataset_id": dataset_id,
                        "omitted_gene": omitted_gene,
                        "n_genes_full": len(available_genes),
                        "n_genes_retained": len(retained),
                        "full_auroc": method_auc["z"],
                        "leave_one_gene_out_auroc": loo_auc,
                        "delta_auroc": delta,
                        "direction_retained": loo_auc > 0.5,
                    }
                )
            influence = max(loo_rows[-len(available_genes):], key=lambda row: abs(float(row["delta_auroc"])))
            most_influential_gene = str(influence["omitted_gene"])
            most_influential_delta = float(influence["delta_auroc"])

        published_signs = tuple(1 if weights[gene] > 0 else -1 for gene in available_genes)
        absolute_weights = np.asarray([abs(weights[gene]) for gene in available_genes], dtype=float)
        absolute_weights = absolute_weights / absolute_weights.sum()
        z_matrix = scaled_by_dataset[dataset_id]["z"].loc[sample_order, available_genes].to_numpy(dtype=float)
        orientation_aucs: list[float] = []
        for pattern_index, signs in enumerate(itertools.product([-1, 1], repeat=len(available_genes))):
            orientation_score = z_matrix @ (absolute_weights * np.asarray(signs, dtype=float))
            orientation_auc = float(roc_auc_score(labels, orientation_score))
            orientation_aucs.append(orientation_auc)
            sign_rows.append(
                {
                    "signature_id": signature_id,
                    "task_id": task_id,
                    "dataset_id": dataset_id,
                    "pattern_index": pattern_index,
                    "sign_pattern": ";".join(f"{gene}:{sign:+d}" for gene, sign in zip(available_genes, signs)),
                    "is_published_orientation": signs == published_signs,
                    "auroc": orientation_auc,
                }
            )
        orientation_array = np.asarray(orientation_aucs)
        midpoint_percentile = 100 * (
            np.sum(orientation_array < method_auc["z"]) + 0.5 * np.sum(np.isclose(orientation_array, method_auc["z"]))
        ) / len(orientation_array)

        cell_rows.append(
            {
                "signature_id": signature_id,
                "task_id": task_id,
                "dataset_id": dataset_id,
                "trial_id": cell.trial_id,
                "evaluation_role": cell.evaluation_role,
                "same_therapy": bool(cell.same_therapy),
                "same_endpoint": bool(cell.same_endpoint),
                "same_disease": bool(cell.same_disease),
                "alignment_score": int(cell.same_therapy) + int(cell.same_endpoint) + int(cell.same_disease),
                "n_subjects": int(cell.n_subjects),
                "n_genes_used": len(available_genes),
                "frozen_auroc": float(cell.auroc),
                "z_recomputed_auroc": method_auc["z"],
                "mad_auroc": method_auc["mad"],
                "rank_auroc": method_auc["rank"],
                "scaling_auroc_min": min(method_auc.values()),
                "scaling_auroc_max": max(method_auc.values()),
                "scaling_max_abs_delta": max(abs(value - method_auc["z"]) for value in method_auc.values()),
                "all_scalings_same_direction": all(value > 0.5 for value in method_auc.values())
                or all(value < 0.5 for value in method_auc.values()),
                "frozen_classification": cell.transport_classification,
                "bootstrap_probability_auc_gt_0_5": probability,
                "jackknife_min_auroc": float(jackknife.min()),
                "jackknife_max_auroc": float(jackknife.max()),
                "jackknife_max_abs_delta": float(np.max(np.abs(jackknife - method_auc["z"]))),
                "loo_gene_min_auroc": min(loo_auc_values) if loo_auc_values else math.nan,
                "loo_gene_max_auroc": max(loo_auc_values) if loo_auc_values else math.nan,
                "loo_gene_direction_retention_fraction": float(np.mean(np.asarray(loo_auc_values) > 0.5))
                if loo_auc_values
                else math.nan,
                "most_influential_gene": most_influential_gene,
                "most_influential_gene_delta_auroc": most_influential_delta,
                "orientation_patterns": len(orientation_array),
                "published_orientation_percentile": midpoint_percentile,
                "orientation_null_min_auroc": float(orientation_array.min()),
                "orientation_null_median_auroc": float(np.median(orientation_array)),
                "orientation_null_max_auroc": float(orientation_array.max()),
            }
        )

    cells = pd.DataFrame(cell_rows)
    cells.to_csv(OUT / "phase3_cell_robustness.csv", index=False)
    pd.DataFrame(loo_rows).to_csv(OUT / "phase3_leave_one_gene_out.csv", index=False)
    pd.DataFrame(sign_rows).to_csv(OUT / "phase3_orientation_enumeration.csv", index=False)
    pd.DataFrame(bootstrap_rows).to_csv(OUT / "phase3_bootstrap_direction_stability.csv", index=False)

    context_group_rows: list[dict] = []
    for score, group in cells.groupby("alignment_score", sort=True):
        context_group_rows.append(
            {
                "alignment_score": int(score),
                "n_cells": len(group),
                "median_auroc": float(group["frozen_auroc"].median()),
                "q1_auroc": float(group["frozen_auroc"].quantile(0.25)),
                "q3_auroc": float(group["frozen_auroc"].quantile(0.75)),
                "n_supported": int((group["frozen_classification"] == "SUPPORTED_DIRECTION").sum()),
                "support_fraction": float((group["frozen_classification"] == "SUPPORTED_DIRECTION").mean()),
            }
        )
    pd.DataFrame(context_group_rows).to_csv(OUT / "phase3_context_alignment_groups.csv", index=False)
    context_trend = two_way_context_bootstrap(cells.rename(columns={"frozen_auroc": "auroc"}))
    (OUT / "phase3_context_alignment_trend.json").write_text(json.dumps(context_trend, indent=2), encoding="utf-8")

    contrast_rows: list[dict] = []
    for signature_id, group in cells.groupby("signature_id"):
        for variable in ["same_therapy", "same_endpoint", "same_disease"]:
            if group[variable].nunique() == 2:
                aligned = float(group.loc[group[variable], "frozen_auroc"].median())
                mismatched = float(group.loc[~group[variable], "frozen_auroc"].median())
                contrast_rows.append(
                    {
                        "signature_id": signature_id,
                        "context_component": variable,
                        "aligned_n": int(group[variable].sum()),
                        "mismatched_n": int((~group[variable]).sum()),
                        "aligned_median_auroc": aligned,
                        "mismatched_median_auroc": mismatched,
                        "aligned_minus_mismatched": aligned - mismatched,
                    }
                )
    pd.DataFrame(contrast_rows).to_csv(OUT / "phase3_within_program_context_contrasts.csv", index=False)

    # Program architecture.
    program_ids = sorted(inventory.index)
    all_genes = sorted(gene_table.loc[gene_table["signature_id"].isin(program_ids), "gene_symbol"].unique())
    architecture = pd.DataFrame(0.0, index=program_ids, columns=all_genes)
    for signature_id in program_ids:
        genes = gene_table.loc[gene_table["signature_id"] == signature_id, "gene_symbol"].tolist()
        weights = raw_program_weights(signature_id, genes, inventory, gene_table)
        for gene, weight in weights.items():
            architecture.loc[signature_id, gene] = weight
    architecture.index.name = "signature_id"
    architecture.reset_index().to_csv(OUT / "phase3_program_gene_weight_matrix.csv", index=False)

    pair_rows: list[dict] = []
    for first, second in itertools.combinations(program_ids, 2):
        first_genes = set(architecture.columns[architecture.loc[first] != 0])
        second_genes = set(architecture.columns[architecture.loc[second] != 0])
        shared = sorted(first_genes & second_genes)
        union = first_genes | second_genes
        directional_agreement = (
            float(np.mean([np.sign(architecture.loc[first, gene]) == np.sign(architecture.loc[second, gene]) for gene in shared]))
            if shared
            else math.nan
        )
        pair_rows.append(
            {
                "signature_1": first,
                "signature_2": second,
                "n_genes_1": len(first_genes),
                "n_genes_2": len(second_genes),
                "n_shared_genes": len(shared),
                "shared_genes": ";".join(shared),
                "jaccard_similarity": len(shared) / len(union),
                "shared_gene_direction_agreement": directional_agreement,
            }
        )
    pd.DataFrame(pair_rows).to_csv(OUT / "phase3_program_pair_gene_similarity.csv", index=False)

    subject_scores = pd.read_csv(ROOT / "results" / "performance" / "phase1_v1" / "subject_level_scores.csv")
    score_correlation_rows: list[dict] = []
    for task_id, group in subject_scores.groupby("task_id"):
        pivot = group.pivot(index="subject_id", columns="signature_id", values="score")
        correlation = pivot.corr(method="spearman")
        for first, second in itertools.combinations(sorted(correlation.columns), 2):
            score_correlation_rows.append(
                {
                    "task_id": task_id,
                    "signature_1": first,
                    "signature_2": second,
                    "n_subjects": int(pivot[[first, second]].dropna().shape[0]),
                    "spearman_rho": float(correlation.loc[first, second]),
                }
            )
    score_correlations = pd.DataFrame(score_correlation_rows)
    score_correlations.to_csv(OUT / "phase3_task_score_correlations.csv", index=False)
    median_correlations = (
        score_correlations.groupby(["signature_1", "signature_2"], as_index=False)
        .agg(n_shared_tasks=("task_id", "nunique"), median_spearman_rho=("spearman_rho", "median"))
    )
    median_correlations.to_csv(OUT / "phase3_program_pair_median_score_correlations.csv", index=False)

    performance_rows: list[dict] = []
    performance_pivot = cells.pivot(index="task_id", columns="signature_id", values="frozen_auroc")
    for first, second in itertools.combinations(sorted(performance_pivot.columns), 2):
        shared = performance_pivot[[first, second]].dropna()
        if len(shared) >= 3:
            rho = float(spearmanr(shared[first], shared[second]).statistic)
        else:
            rho = math.nan
        performance_rows.append(
            {
                "signature_1": first,
                "signature_2": second,
                "n_shared_tasks": len(shared),
                "performance_profile_spearman_rho": rho,
            }
        )
    pd.DataFrame(performance_rows).to_csv(OUT / "phase3_program_performance_profile_correlations.csv", index=False)

    # UNIFI shared-subject endpoint drift.
    mh = task_manifests["T_UNIFI_MH_W8"][["subject_id", "outcome_label"]].rename(columns={"outcome_label": "mucosal_healing"})
    remission = task_manifests["T_UNIFI_REM_W8"][["subject_id", "outcome_label"]].rename(columns={"outcome_label": "clinical_remission"})
    unifi = mh.merge(remission, on="subject_id", how="inner")
    contingency = pd.crosstab(unifi["mucosal_healing"], unifi["clinical_remission"], dropna=False)
    contingency_rows = []
    for mh_label in [0, 1]:
        for remission_label in [0, 1]:
            contingency_rows.append(
                {
                    "mucosal_healing": mh_label,
                    "clinical_remission": remission_label,
                    "n_subjects": int(contingency.reindex(index=[0, 1], columns=[0, 1], fill_value=0).loc[mh_label, remission_label]),
                }
            )
    pd.DataFrame(contingency_rows).to_csv(OUT / "phase3_unifi_endpoint_contingency.csv", index=False)
    unifi_summary = {
        "overlapping_subjects": len(unifi),
        "raw_agreement": float(np.mean(unifi["mucosal_healing"] == unifi["clinical_remission"])),
        "cohen_kappa": float(cohen_kappa_score(unifi["mucosal_healing"], unifi["clinical_remission"])),
        "mucosal_healing_prevalence": float(unifi["mucosal_healing"].mean()),
        "clinical_remission_prevalence": float(unifi["clinical_remission"].mean()),
    }
    (OUT / "phase3_unifi_endpoint_summary.json").write_text(json.dumps(unifi_summary, indent=2), encoding="utf-8")
    unifi_auc = cells[cells["task_id"].isin(["T_UNIFI_MH_W8", "T_UNIFI_REM_W8"])].pivot(
        index="signature_id", columns="task_id", values="frozen_auroc"
    ).dropna()
    unifi_auc = unifi_auc.rename(
        columns={"T_UNIFI_MH_W8": "mucosal_healing_auroc", "T_UNIFI_REM_W8": "clinical_remission_auroc"}
    )
    unifi_auc["clinical_remission_minus_mucosal_healing"] = (
        unifi_auc["clinical_remission_auroc"] - unifi_auc["mucosal_healing_auroc"]
    )
    unifi_auc.reset_index().to_csv(OUT / "phase3_unifi_paired_program_aurocs.csv", index=False)

    # Strict-cell visual source data and exact reproduction of frozen permutation sequences.
    strict = frozen[frozen["evaluation_role"] == "SAME_THERAPY_SAME_ENDPOINT_EXTERNAL"]
    strict_roc_rows: list[dict] = []
    strict_permutation_rows: list[dict] = []
    strict_score_rows: list[dict] = []
    for cell in strict.itertuples(index=False):
        group = subject_scores[(subject_scores["signature_id"] == cell.signature_id) & (subject_scores["task_id"] == cell.task_id)].copy()
        labels = group["outcome_label"].astype(int).tolist()
        scores = group["score"].astype(float).tolist()
        fpr, tpr, thresholds = roc_curve(labels, scores)
        for index, (x, y, threshold) in enumerate(zip(fpr, tpr, thresholds)):
            strict_roc_rows.append(
                {
                    "signature_id": cell.signature_id,
                    "task_id": cell.task_id,
                    "point_index": index,
                    "false_positive_rate": x,
                    "true_positive_rate": y,
                    "threshold": threshold,
                }
            )
        for row in group.itertuples(index=False):
            strict_score_rows.append(
                {
                    "signature_id": cell.signature_id,
                    "task_id": cell.task_id,
                    "subject_id": row.subject_id,
                    "outcome_label": int(row.outcome_label),
                    "score": float(row.score),
                }
            )
        rng = random.Random(stable_seed(cell.signature_id, cell.task_id, offset=1))
        # Match the Phase 1 seed definition exactly, including its token format.
        phase1_token = hashlib.sha256(f"20260831::{cell.signature_id}::{cell.task_id}".encode()).hexdigest()[:8]
        rng = random.Random(int(phase1_token, 16) + 1)
        permuted = labels[:]
        for replicate in range(N_STRICT_PERM):
            rng.shuffle(permuted)
            strict_permutation_rows.append(
                {
                    "signature_id": cell.signature_id,
                    "task_id": cell.task_id,
                    "replicate": replicate + 1,
                    "permuted_auroc": float(roc_auc_score(permuted, scores)),
                }
            )
    pd.DataFrame(strict_roc_rows).to_csv(OUT / "phase3_strict_roc_curves.csv", index=False)
    pd.DataFrame(strict_score_rows).to_csv(OUT / "phase3_strict_subject_scores.csv", index=False)
    pd.DataFrame(strict_permutation_rows).to_csv(OUT / "phase3_strict_permutation_nulls.csv", index=False)

    summary = {
        "version": "phase3_robustness_v1",
        "source_submission_commit": "66d31f95f8cbb7eafaaf0e32e7aac8132f049990",
        "plan_commit": "8bb1e6c2c79079abd1348202dd1e83929fb5c4db",
        "executed_cells": len(cells),
        "multi_gene_leave_one_out_rows": len(loo_rows),
        "orientation_enumeration_rows": len(sign_rows),
        "all_z_aurocs_match_frozen": bool(np.allclose(cells["frozen_auroc"], cells["z_recomputed_auroc"], atol=5e-6)),
        "all_three_scalings_same_direction_cells": int(cells["all_scalings_same_direction"].sum()),
        "median_scaling_max_abs_delta": float(cells["scaling_max_abs_delta"].median()),
        "max_scaling_max_abs_delta": float(cells["scaling_max_abs_delta"].max()),
        "median_jackknife_max_abs_delta": float(cells["jackknife_max_abs_delta"].median()),
        "max_jackknife_max_abs_delta": float(cells["jackknife_max_abs_delta"].max()),
        "context_alignment": context_trend,
        "unifi_endpoint": unifi_summary,
        "new_confirmatory_tests": 0,
        "phase1_outputs_modified": False,
    }
    (OUT / "phase3_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
