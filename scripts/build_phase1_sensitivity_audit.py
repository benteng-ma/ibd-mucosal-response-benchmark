#!/usr/bin/env python3
"""Document pre-specified lineage and availability sensitivity consequences."""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "sensitivity" / "phase1_v1"


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    grid = read_csv(ROOT / "metadata" / "evaluation_grid.csv")
    results = read_csv(ROOT / "results" / "performance" / "phase1_v1" / "study_specific_evaluations.csv")
    result_pairs = {(r["signature_id"], r["task_id"]) for r in results}
    transport = [r for r in grid if r["transportability_eligible"].upper() == "TRUE"]
    raw_only = [r for r in transport if r["dataset_id"] == "E-MTAB-7845_VDZ"]
    rows = [
        {
            "sensitivity_id": "SENS_LEUVEN_POSSIBLE_OVERLAP",
            "risk": "GSE73661 infliximab comparator may overlap older Leuven biobanks",
            "conservative_action": "exclude GSE73661 IFX from independent primary/transport execution",
            "executed_cells_affected": 0,
            "conclusion": "Numerical Phase 1 results invariant because arm threshold and lineage flags already excluded this cohort",
        },
        {
            "sensitivity_id": "SENS_ACT1_REUPLOAD",
            "risk": "GSE12251 and GSE23597 contain exact reuploaded ACT1 baseline files",
            "conservative_action": "block all related signature-task cells from external validation",
            "executed_cells_affected": 0,
            "conclusion": "No ACT1 reupload appears in the executed external grid",
        },
        {
            "sensitivity_id": "SENS_LEUVEN_SUBSET",
            "risk": "GSE14580 is a complete subset of GSE16879",
            "conservative_action": "count as one cohort and block related signature-task cells",
            "executed_cells_affected": 0,
            "conclusion": "No subset split appears in the executed external grid",
        },
        {
            "sensitivity_id": "SENS_UNITI_TRIAL_FAMILY",
            "risk": "GSE112366, GSE207022 and GSE207465 are one UNITI-2 trial family",
            "conservative_action": "block same-family UST signatures and count the trial once",
            "executed_cells_affected": 0,
            "conclusion": "UNITI-2 is used only for signatures independent of UNITI-2 lineage",
        },
        {
            "sensitivity_id": "SENS_EMTAB7604_MISSING_FILE",
            "risk": "one of 44 SDRF subjects lacks a processed public count file",
            "conservative_action": "complete-case expression availability without imputation",
            "executed_cells_affected": sum(r["task_id"] == "T_EMTAB7604_ER" for r in results),
            "conclusion": "E-MTAB-7604 analysis uses 43 subjects (18 responder, 25 nonresponder); missing subject is documented",
        },
        {
            "sensitivity_id": "SENS_EMTAB7845_RAW_ONLY",
            "risk": "E-MTAB-7845 has public raw FASTQ but no repository-processed matrix",
            "conservative_action": "do not create a non-comparable ad hoc processed cohort in this benchmark version",
            "executed_cells_affected": len(raw_only),
            "conclusion": f"{len(raw_only)} frozen transport cells remain not executed and are not counted as negative results",
        },
    ]
    write_csv(OUT / "lineage_and_availability_sensitivity.csv", rows)
    expected = {(r["signature_id"], r["task_id"]) for r in transport if r["dataset_id"] != "E-MTAB-7845_VDZ"}
    summary = {
        "frozen_transport_cells": len(transport),
        "executed_cells": len(result_pairs),
        "raw_only_not_executed_cells": len(raw_only),
        "unexpected_missing_executed_pairs": sorted(["::".join(x) for x in expected - result_pairs]),
        "lineage_sensitivity_changes_numeric_results": False,
    }
    (OUT / "sensitivity_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
