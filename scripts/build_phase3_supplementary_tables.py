#!/usr/bin/env python3
"""Create manuscript-facing Phase 3 supplementary tables from locked results."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "phase3_robustness" / "v1"
OUT = ROOT / "manuscript" / "supplement"


COPIES = {
    "table_s6_phase3_cell_robustness.csv": "phase3_cell_robustness.csv",
    "table_s7_leave_one_gene_out.csv": "phase3_leave_one_gene_out.csv",
    "table_s8_orientation_enumeration.csv": "phase3_orientation_enumeration.csv",
    "table_s9_context_alignment_groups.csv": "phase3_context_alignment_groups.csv",
    "table_s10_within_program_context_contrasts.csv": "phase3_within_program_context_contrasts.csv",
    "table_s11_program_gene_similarity.csv": "phase3_program_pair_gene_similarity.csv",
}


def copy_text_csv(source: Path, destination: Path) -> int:
    text = source.read_text(encoding="utf-8-sig")
    destination.write_text(text, encoding="utf-8", newline="")
    with destination.open(encoding="utf-8-sig", newline="") as handle:
        return sum(1 for _ in csv.DictReader(handle))


def write_combined_correlations(destination: Path) -> int:
    rows: list[dict[str, str]] = []
    sources = [
        ("subject_score", RESULTS / "phase3_program_pair_median_score_correlations.csv", "median_spearman_rho"),
        ("performance_profile", RESULTS / "phase3_program_performance_profile_correlations.csv", "performance_profile_spearman_rho"),
    ]
    for correlation_type, source, value_column in sources:
        with source.open(encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle):
                rows.append(
                    {
                        "correlation_type": correlation_type,
                        "signature_1": row["signature_1"],
                        "signature_2": row["signature_2"],
                        "n_shared_tasks": row["n_shared_tasks"],
                        "spearman_rho": row[value_column],
                    }
                )
    with destination.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)


def write_unifi_summary(destination: Path) -> int:
    summary = json.loads((RESULTS / "phase3_unifi_endpoint_summary.json").read_text(encoding="utf-8"))
    rows = [
        {"record_type": "endpoint_summary", "item": key, "value": value, "value_2": ""}
        for key, value in summary.items()
    ]
    with (RESULTS / "phase3_unifi_endpoint_contingency.csv").open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            rows.append(
                {
                    "record_type": "endpoint_contingency",
                    "item": f"MH={row['mucosal_healing']};CR={row['clinical_remission']}",
                    "value": row["n_subjects"],
                    "value_2": "",
                }
            )
    with (RESULTS / "phase3_unifi_paired_program_aurocs.csv").open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            rows.append(
                {
                    "record_type": "paired_program_auroc",
                    "item": row["signature_id"],
                    "value": row["mucosal_healing_auroc"],
                    "value_2": row["clinical_remission_auroc"],
                }
            )
    with destination.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["record_type", "item", "value", "value_2"])
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)


def write_figure_source_manifest(destination: Path) -> int:
    source_root = ROOT / "manuscript" / "figure_source_data" / "phase3"
    rows = []
    for path in sorted(source_root.glob("*")):
        if not path.is_file():
            continue
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        rows.append(
            {
                "file": path.name,
                "bytes": path.stat().st_size,
                "sha256": digest,
                "purpose": "source data for Phase 3 main or supplementary figures",
            }
        )
    with destination.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    counts = {}
    for destination, source in COPIES.items():
        counts[destination] = copy_text_csv(RESULTS / source, OUT / destination)
    counts["table_s12_program_pair_correlations.csv"] = write_combined_correlations(
        OUT / "table_s12_program_pair_correlations.csv"
    )
    counts["table_s13_unifi_endpoint_sensitivity.csv"] = write_unifi_summary(
        OUT / "table_s13_unifi_endpoint_sensitivity.csv"
    )
    counts["table_s14_figure_source_manifest.csv"] = write_figure_source_manifest(
        OUT / "table_s14_figure_source_manifest.csv"
    )
    print(json.dumps(counts, indent=2))


if __name__ == "__main__":
    main()
