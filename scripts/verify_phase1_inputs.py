#!/usr/bin/env python3
"""Verify Phase 1 input hashes, frozen files and subject/sample constraints."""

from __future__ import annotations

import csv
import hashlib
import json
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE_COMMIT = "e30bd8f60446c0a61744a3be529de3ed6d135cbc"


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    provenance = json.loads((ROOT / "results" / "preprocessing" / "phase1_expression_provenance.json").read_text(encoding="utf-8"))
    hash_checks = []
    for item in provenance["datasets"]:
        path = ROOT / item["source_file"]
        observed = sha256(path)
        hash_checks.append({
            "dataset_id": item["dataset_id"],
            "source_file": item["source_file"],
            "expected_sha256": item["source_sha256"],
            "observed_sha256": observed,
            "pass": observed == item["source_sha256"],
        })

    frozen_files = [
        "metadata/evaluation_grid.csv",
        "metadata/task_manifest.csv",
        "metadata/data_lineage_matrix.csv",
        "metadata/endpoint_dictionary.csv",
        "metadata/endpoint_mapping.csv",
        "literature/signature_inventory.csv",
        "literature/signature_gene_table.csv",
    ]
    expected = json.loads((ROOT / "provenance" / "frozen_input_sha256.json").read_text())
    diff = [name for name in frozen_files if sha256(ROOT / name) != expected[name]]

    task_summary = read_csv(ROOT / "metadata" / "phase1" / "phase1_task_summary.csv")
    manifest_checks = []
    for task in task_summary:
        rows = read_csv(ROOT / "metadata" / "phase1" / f"{task['task_id']}_analysis_manifest.csv")
        expression = {r["sample_id"] for r in read_csv(ROOT / rows[0]["expression_file"])}
        subjects = [r["subject_id"] for r in rows]
        samples = [r["sample_id"] for r in rows]
        manifest_checks.append({
            "task_id": task["task_id"],
            "n_rows": len(rows),
            "n_unique_subjects": len(set(subjects)),
            "all_samples_in_expression": all(s in expression for s in samples),
            "no_blank_subjects": all(bool(s) for s in subjects),
            "both_outcome_arms": {r["outcome_label"] for r in rows} == {"0", "1"},
            "pass": len(rows) == len(set(subjects)) and all(s in expression for s in samples) and all(bool(s) for s in subjects),
        })

    result = {
        "source_commit": SOURCE_COMMIT,
        "hash_checks": hash_checks,
        "frozen_files_checked": frozen_files,
        "frozen_files_modified_since_phase0": diff,
        "manifest_checks": manifest_checks,
        "all_hashes_pass": all(x["pass"] for x in hash_checks),
        "frozen_inputs_unchanged": not diff,
        "all_manifests_pass": all(x["pass"] for x in manifest_checks),
    }
    result["overall_pass"] = result["all_hashes_pass"] and result["frozen_inputs_unchanged"] and result["all_manifests_pass"]
    path = ROOT / "results" / "qc" / "phase1_input_verification.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2), flush=True)
    if not result["overall_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
