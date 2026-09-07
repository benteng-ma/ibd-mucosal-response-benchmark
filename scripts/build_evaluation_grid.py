#!/usr/bin/env python3
"""Freeze the Phase 0 signature-by-task evaluation grid without outcomes analysis."""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def therapy_match(signature_therapy: str, task_therapy: str) -> bool:
    sig = signature_therapy.upper()
    task = task_therapy.upper()
    if "ANTI_TNF" in sig and task.startswith("ANTI_TNF"):
        return True
    return any(token.strip() == task for token in sig.split(";") if token.strip())


def disease_match(signature_disease: str, task_disease: str) -> bool:
    if signature_disease.upper() == "IBD":
        return True
    sig = set(signature_disease.upper().split(";"))
    task = set(task_disease.upper().split(";"))
    return bool(sig & task)


def lineage_class(signature: dict, dataset: dict) -> str:
    accession = dataset["accession"].upper()
    source_text = ";".join([
        signature["discovery_datasets"],
        signature["datasets_used_in_original_paper"],
    ]).upper()
    if accession != "NONE LOCATED" and accession in source_text:
        return "REUSED_ORIGINAL_DATA"
    related_groups = [
        {"GSE12251", "GSE23597"},
        {"GSE14580", "GSE16879"},
        {"GSE112366", "GSE207022", "GSE207465"},
    ]
    for group in related_groups:
        if accession in group and any(member in source_text for member in group):
            return "SAME_TRIAL_OR_REUPLOAD_RELATED"
    return "INDEPENDENT_TRIAL_OR_COHORT"


def main() -> None:
    signatures = read_csv(ROOT / "literature" / "signature_inventory.csv")
    tasks = read_csv(ROOT / "metadata" / "task_manifest.csv")
    datasets = {r["dataset_id"]: r for r in read_csv(ROOT / "metadata" / "dataset_manifest.csv")}
    coverage = {
        (r["signature_id"], r["dataset_id"]): r
        for r in read_csv(ROOT / "metadata" / "signature_coverage_matrix.csv")
    }

    rows: list[dict] = []
    for sig in signatures:
        for task in tasks:
            dataset = datasets[task["dataset_id"]]
            cov = coverage[(sig["signature_id"], task["dataset_id"])]
            lineage = lineage_class(sig, dataset)
            same_therapy = therapy_match(sig["therapy"], task["therapy_class"])
            same_endpoint = sig["endpoint_family"] == task["endpoint_family"]
            same_disease = disease_match(sig["disease"], task["disease"])
            same_tissue = sig["sample_compartment"] == "mucosa" and task["sample_compartment"] == "mucosa"
            adequate_arms = int(task["n_responder"]) >= 15 and int(task["n_nonresponder"]) >= 15
            executable_task = (
                task["baseline_or_on_treatment"] == "baseline"
                and task["primary_or_secondary"] == "primary"
                and dataset["outcome_available"].upper() == "TRUE"
                and adequate_arms
                and same_tissue
            )
            signature_ready = (
                sig["main_grid_eligible"].upper() == "TRUE"
                and sig["reconstruction_status"] == "COMPLETE"
                and cov["primary_coverage_pass"] == "TRUE"
            )
            independent = lineage == "INDEPENDENT_TRIAL_OR_COHORT"
            primary_validation = all([executable_task, signature_ready, independent, same_therapy, same_endpoint, same_disease])
            transportability = all([executable_task, signature_ready, independent])
            if not independent:
                role = "REUSED_OR_RELATED_DATA_AUDIT_ONLY"
            elif primary_validation:
                role = "SAME_THERAPY_SAME_ENDPOINT_EXTERNAL"
            elif transportability and same_therapy:
                role = "WITHIN_THERAPY_ENDPOINT_OR_DISEASE_TRANSPORT"
            elif transportability:
                role = "CROSS_THERAPY_TRANSPORT_STRESS_TEST"
            else:
                role = "NOT_EXECUTABLE_OR_SECONDARY"
            exclusions: list[str] = []
            if not adequate_arms: exclusions.append("ARM_LT_15")
            if task["primary_or_secondary"] != "primary": exclusions.append("SECONDARY_TASK")
            if task["baseline_or_on_treatment"] != "baseline": exclusions.append("NOT_BASELINE")
            if dataset["outcome_available"].upper() != "TRUE": exclusions.append("OUTCOME_MAPPING_UNAVAILABLE")
            if not same_tissue: exclusions.append("COMPARTMENT_MISMATCH")
            if sig["main_grid_eligible"].upper() != "TRUE": exclusions.append("SIGNATURE_NOT_MAIN_GRID")
            if sig["reconstruction_status"] != "COMPLETE": exclusions.append("SIGNATURE_INCOMPLETE")
            if cov["primary_coverage_pass"] != "TRUE": exclusions.append("COVERAGE_NOT_CONFIRMED")
            if not independent: exclusions.append(lineage)
            rows.append({
                "signature_id": sig["signature_id"],
                "task_id": task["task_id"],
                "dataset_id": task["dataset_id"],
                "trial_id": task["trial_id"],
                "lineage_class": lineage,
                "same_therapy": str(same_therapy).upper(),
                "same_endpoint": str(same_endpoint).upper(),
                "same_disease": str(same_disease).upper(),
                "same_tissue_compartment": str(same_tissue).upper(),
                "coverage_fraction": cov["coverage_fraction"],
                "coverage_status": cov["mapping_status"],
                "adequate_task_arms": str(adequate_arms).upper(),
                "evaluation_role": role,
                "primary_validation_eligible": str(primary_validation).upper(),
                "transportability_eligible": str(transportability).upper(),
                "exclusion_reasons": ";".join(exclusions),
                "phase0_note": "Frozen eligibility only; no scores, predictions, metrics, intervals or hypothesis tests.",
            })

    path = ROOT / "metadata" / "evaluation_grid.csv"
    write_csv(path, rows)
    summary = {
        "n_grid_rows": len(rows),
        "n_primary_same_therapy_same_endpoint_external_opportunities": sum(r["primary_validation_eligible"] == "TRUE" for r in rows),
        "n_transportability_opportunities": sum(r["transportability_eligible"] == "TRUE" for r in rows),
        "n_cross_therapy_transport_opportunities": sum(r["evaluation_role"] == "CROSS_THERAPY_TRANSPORT_STRESS_TEST" for r in rows),
        "n_independent_trials_in_primary_external_grid": len({r["trial_id"] for r in rows if r["primary_validation_eligible"] == "TRUE"}),
        "n_independent_trials_in_transportability_grid": len({r["trial_id"] for r in rows if r["transportability_eligible"] == "TRUE"}),
        "n_signatures_in_transportability_grid": len({r["signature_id"] for r in rows if r["transportability_eligible"] == "TRUE"}),
        "performance_computed": False,
        "phase0_only": True,
    }
    output = ROOT / "results" / "evaluation_grid"
    output.mkdir(parents=True, exist_ok=True)
    (output / "evaluation_grid_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    with (output / "eligible_primary_opportunities.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        eligible = [r for r in rows if r["primary_validation_eligible"] == "TRUE"]
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(eligible)


if __name__ == "__main__":
    main()
