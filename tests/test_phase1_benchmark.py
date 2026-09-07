#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import math
import sys
import unittest
from collections import defaultdict
from pathlib import Path

from sklearn.metrics import average_precision_score, roc_auc_score


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from run_phase1_benchmark import average_precision, auroc  # noqa: E402


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


class Phase1BenchmarkTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not (ROOT / "results/performance/phase1_v1/subject_level_scores.csv").exists():
            raise unittest.SkipTest("Subject-level inputs intentionally excluded; run the local reproduction pipeline first.")
        cls.result_dir = ROOT / "results" / "performance" / "phase1_v1"
        cls.results = read_csv(cls.result_dir / "study_specific_evaluations.csv")
        cls.scores = read_csv(cls.result_dir / "subject_level_scores.csv")

    def test_frozen_inputs_and_hashes(self) -> None:
        qc = json.loads((ROOT / "results" / "qc" / "phase1_input_verification.json").read_text(encoding="utf-8"))
        self.assertTrue(qc["overall_pass"])
        self.assertEqual(qc["frozen_files_modified_since_phase0"], [])

    def test_expected_grid_partition(self) -> None:
        summary = json.loads((self.result_dir / "benchmark_summary.json").read_text(encoding="utf-8"))
        self.assertEqual(summary["n_frozen_executable_grid_cells"], 46)
        self.assertEqual(summary["n_executed_cells"], 46)
        self.assertEqual(summary["n_strict_external_cells"], 3)
        self.assertEqual(summary["n_cross_therapy_stress_cells"], 29)
        self.assertFalse(summary["model_fitting"])
        self.assertFalse(summary["feature_selection"])
        self.assertFalse(summary["cross_cohort_merging"])

    def test_metrics_match_sklearn(self) -> None:
        grouped: dict[tuple[str, str], list[dict]] = defaultdict(list)
        for row in self.scores:
            grouped[(row["signature_id"], row["task_id"])].append(row)
        for result in self.results:
            rows = grouped[(result["signature_id"], result["task_id"])]
            labels = [int(r["outcome_label"]) for r in rows]
            scores = [float(r["score"]) for r in rows]
            self.assertAlmostEqual(auroc(labels, scores), roc_auc_score(labels, scores), places=10)
            self.assertAlmostEqual(average_precision(labels, scores), average_precision_score(labels, scores), places=10)
            self.assertAlmostEqual(float(result["auroc"]), roc_auc_score(labels, scores), places=6)
            self.assertAlmostEqual(float(result["average_precision"]), average_precision_score(labels, scores), places=6)

    def test_subject_level_and_coverage_constraints(self) -> None:
        for task in {r["task_id"] for r in self.results}:
            manifest = read_csv(ROOT / "metadata" / "phase1" / f"{task}_analysis_manifest.csv")
            self.assertEqual(len(manifest), len({r["subject_id"] for r in manifest}))
        self.assertTrue(all(float(r["actual_gene_coverage"]) >= 0.8 for r in self.results))
        self.assertTrue(all(r["score_orientation"] == "higher predicts responder" for r in self.scores))

    def test_strict_family_and_holm_values(self) -> None:
        strict = read_csv(self.result_dir / "strict_external_evaluations.csv")
        self.assertEqual(len(strict), 3)
        self.assertTrue(all(r["permutation_p_value"] for r in strict))
        self.assertTrue(all(r["holm_adjusted_p_value"] for r in strict))
        self.assertTrue(all(float(r["holm_adjusted_p_value"]) >= float(r["permutation_p_value"]) for r in strict))
        self.assertTrue(all(float(r["holm_adjusted_p_value"]) >= 0.05 for r in strict))

    def test_negative_matrix_is_explicit(self) -> None:
        negative = read_csv(self.result_dir / "negative_transportability_matrix.csv")
        self.assertEqual(len(negative), 27)
        self.assertTrue(all(r["transport_classification"] in {"UNCERTAIN_OR_NONTRANSPORTABLE", "REVERSED_DIRECTION"} for r in negative))
        self.assertTrue(all(r["calibration_status"] == "NOT_APPLICABLE_UNCALIBRATED_PROGRAM_SCORE" for r in self.results))

    def test_hard_exclusions(self) -> None:
        datasets = {r["dataset_id"] for r in self.results}
        self.assertNotIn("E-MTAB-7845_VDZ", datasets)
        self.assertFalse(any(r["lineage_class"] != "INDEPENDENT_TRIAL_OR_COHORT" for r in self.results))
        self.assertFalse(any(math.isnan(float(r["auroc"])) for r in self.results))


if __name__ == "__main__":
    unittest.main(verbosity=2)
