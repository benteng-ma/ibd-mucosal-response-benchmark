#!/usr/bin/env python3
"""Create frozen subject-level Phase 1 task manifests from public metadata."""

from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
GEO = ROOT / "data" / "raw" / "geo_metadata"
AE = ROOT / "data" / "raw" / "arrayexpress_metadata"
EXPR = ROOT / "data" / "processed" / "phase1_v1"
OUT = ROOT / "metadata" / "phase1"


def read_csv(path: Path, delimiter: str = ",") -> list[dict]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def characteristics(value: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in value.split(" || "):
        if ": " in item:
            key, val = item.split(": ", 1)
            result[key.strip().casefold()] = val.strip()
    return result


def yn(value: str) -> int | None:
    value = value.strip().casefold()
    if value in {"y", "yes", "responder"}:
        return 1
    if value in {"n", "no", "non-responder", "nonresponder"}:
        return 0
    return None


def expression_ids(dataset_id: str) -> set[str]:
    path = EXPR / f"{dataset_id}_signature_genes.csv"
    return {r["sample_id"] for r in read_csv(path)}


def make_row(task: dict, sample_id: str, subject_id: str, label: int,
             endpoint_original: str, disease: str, therapy: str,
             tissue: str, source_accession: str) -> dict:
    return {
        "task_id": task["task_id"],
        "dataset_id": task["dataset_id"],
        "trial_id": task["trial_id"],
        "source_accession": source_accession,
        "sample_id": sample_id,
        "subject_id": subject_id,
        "outcome_label": label,
        "outcome_text": "RESPONDER" if label == 1 else "NONRESPONDER",
        "endpoint_original": endpoint_original,
        "endpoint_family": task["endpoint_family"],
        "endpoint_time": task["endpoint_time"],
        "therapy": therapy,
        "disease": disease,
        "tissue": tissue,
        "baseline_active": "TRUE",
        "expression_file": f"data/processed/phase1_v1/{task['dataset_id']}_signature_genes.csv",
        "inclusion_reason": "frozen baseline active-treatment task with public individual outcome and expression",
    }


def geo_task(task: dict) -> tuple[list[dict], list[dict]]:
    accession = {
        "GSE92415_PURSUIT": "GSE92415",
        "GSE212849_PROGECT": "GSE212849",
        "GSE112366_UNITI2_ILEUM": "GSE112366",
        "GSE206285_UNIFI": "GSE206285",
    }[task["dataset_id"]]
    ids = expression_ids(task["dataset_id"])
    rows: list[dict] = []
    exclusions: list[dict] = []
    for raw in read_csv(GEO / f"{accession}_samples.csv"):
        c = characteristics(raw["Sample_characteristics_ch1"])
        sample_id = raw["sample_id"]
        subject = ""
        label: int | None = None
        endpoint_original = ""
        eligible = False
        reason = "not the frozen baseline active-treatment task"
        disease = c.get("diagnosis", c.get("disease", ""))
        therapy = c.get("treatment", c.get("treatment_induction", ""))
        tissue = c.get("tissue", "")

        if task["dataset_id"] == "GSE92415_PURSUIT":
            subject = c.get("subject", "")
            endpoint_original = c.get("wk6response", "")
            label = yn(endpoint_original)
            eligible = c.get("visit", "").casefold() == "week 0" and c.get("treatment", "").casefold() == "golimumab" and label is not None
        elif task["dataset_id"] == "GSE212849_PROGECT":
            subject = c.get("donor id", "")
            key = "mucosal healing at week 6" if task["endpoint_family"] == "MUCOSAL_HEALING" else "clinical remission at week 6"
            endpoint_original = c.get(key, "")
            label = yn(endpoint_original)
            eligible = c.get("visit", "").casefold() == "week 0" and "golimumab" in c.get("treatment", "").casefold() and label is not None
        elif task["dataset_id"] == "GSE112366_UNITI2_ILEUM":
            subject = c.get("subject", "")
            endpoint_original = c.get("i-wk8 response", "")
            label = yn(endpoint_original)
            eligible = c.get("visit", "").upper() == "I-WK0" and c.get("treatment_induction", "").casefold() == "ust" and label is not None
        elif task["dataset_id"] == "GSE206285_UNIFI":
            subject = c.get("donor id", "")
            key = "mucosal healing at week 8" if task["endpoint_family"] == "MUCOSAL_HEALING" else "clinical remission at week 8"
            endpoint_original = c.get(key, "")
            label = yn(endpoint_original)
            eligible = c.get("visit", "").upper() == "WEEK I-0" and "ustekinumab" in c.get("treatment", "").casefold() and label is not None

        if sample_id not in ids:
            eligible = False
            reason = "expression sample absent from extracted repository matrix"
        if eligible:
            rows.append(make_row(task, sample_id, subject, int(label), endpoint_original, disease, therapy, tissue, accession))
        else:
            exclusions.append({"task_id": task["task_id"], "sample_id": sample_id, "subject_id": subject, "reason": reason})
    return rows, exclusions


def emtab_task(task: dict) -> tuple[list[dict], list[dict]]:
    ids = expression_ids(task["dataset_id"])
    rows: list[dict] = []
    exclusions: list[dict] = []
    sdrf = read_csv(AE / "E-MTAB-7604_E-MTAB-7604.sdrf.txt", delimiter="\t")
    for raw in sdrf:
        sample_id = raw["Derived Array Data File"].strip()
        subject_id = raw["Comment[BioSD_SAMPLE]"].strip()
        original = raw["Characteristics[clinical history]"].strip()
        label = yn(original)
        if sample_id in ids and label is not None:
            rows.append(make_row(
                task, sample_id, subject_id, label, original,
                raw["Characteristics[disease]"], raw["Characteristics[treatment]"],
                raw["Characteristics[organism part]"], "E-MTAB-7604"
            ))
        else:
            reason = "processed expression file absent from public archive" if sample_id not in ids else "outcome unavailable"
            exclusions.append({"task_id": task["task_id"], "sample_id": sample_id or raw["Source Name"], "subject_id": subject_id, "reason": reason})
    return rows, exclusions


def main() -> None:
    tasks = {r["task_id"]: r for r in read_csv(ROOT / "metadata" / "task_manifest.csv")}
    execute = [
        "T_PURSUIT_CR_W6",
        "T_PROGECT_MH_W6",
        "T_EMTAB7604_ER",
        "T_UNITI2_CR_W8",
        "T_UNIFI_MH_W8",
        "T_UNIFI_REM_W8",
    ]
    all_rows: list[dict] = []
    all_exclusions: list[dict] = []
    summaries: list[dict] = []
    fields = [
        "task_id", "dataset_id", "trial_id", "source_accession", "sample_id", "subject_id",
        "outcome_label", "outcome_text", "endpoint_original", "endpoint_family", "endpoint_time",
        "therapy", "disease", "tissue", "baseline_active", "expression_file", "inclusion_reason",
    ]
    for task_id in execute:
        task = tasks[task_id]
        rows, exclusions = emtab_task(task) if task["dataset_id"] == "E-MTAB-7604_ATNF" else geo_task(task)
        if len({r["subject_id"] for r in rows}) != len(rows):
            raise RuntimeError(f"Duplicate subject in task {task_id}")
        write_csv(OUT / f"{task_id}_analysis_manifest.csv", rows, fields)
        counts = Counter(r["outcome_text"] for r in rows)
        summaries.append({
            "task_id": task_id,
            "dataset_id": task["dataset_id"],
            "trial_id": task["trial_id"],
            "n_subjects": len(rows),
            "n_responder": counts["RESPONDER"],
            "n_nonresponder": counts["NONRESPONDER"],
            "n_excluded_repository_samples": len(exclusions),
            "status": "EXECUTABLE" if counts["RESPONDER"] >= 15 and counts["NONRESPONDER"] >= 15 else "SECONDARY_ARM_LT_15",
        })
        all_rows.extend(rows)
        all_exclusions.extend(exclusions)
    write_csv(OUT / "phase1_task_summary.csv", summaries, list(summaries[0]))
    write_csv(OUT / "phase1_sample_exclusions.csv", all_exclusions, ["task_id", "sample_id", "subject_id", "reason"])
    with (OUT / "phase1_manifest_summary.json").open("w", encoding="utf-8") as handle:
        json.dump({"version": "phase1_v1", "tasks": summaries}, handle, indent=2)
    print(json.dumps(summaries, indent=2), flush=True)


if __name__ == "__main__":
    main()
