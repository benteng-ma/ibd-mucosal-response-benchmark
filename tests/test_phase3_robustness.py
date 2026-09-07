from __future__ import annotations

import csv
import json
import unittest
from pathlib import Path

from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "phase3_robustness" / "v1"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


class Phase3RobustnessTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.summary = json.loads((RESULTS / "phase3_summary.json").read_text(encoding="utf-8"))
        cls.cells = read_csv(RESULTS / "phase3_cell_robustness.csv")

    def test_phase1_numerical_lock_is_reproduced(self) -> None:
        frozen = read_csv(ROOT / "results" / "performance" / "phase1_v1" / "study_specific_evaluations.csv")
        expected = {(row["signature_id"], row["task_id"]): float(row["auroc"]) for row in frozen}
        observed = {(row["signature_id"], row["task_id"]): float(row["z_recomputed_auroc"]) for row in self.cells}
        self.assertEqual(set(observed), set(expected))
        for key in expected:
            self.assertAlmostEqual(observed[key], expected[key], places=6)
        self.assertTrue(self.summary["all_z_aurocs_match_frozen"])
        self.assertFalse(self.summary["phase1_outputs_modified"])
        self.assertEqual(self.summary["new_confirmatory_tests"], 0)

    def test_complete_exploratory_result_shapes(self) -> None:
        self.assertEqual(len(self.cells), 46)
        self.assertEqual(len(read_csv(RESULTS / "phase3_leave_one_gene_out.csv")), 189)
        self.assertEqual(len(read_csv(RESULTS / "phase3_orientation_enumeration.csv")), 1920)
        self.assertEqual(self.summary["all_three_scalings_same_direction_cells"], 43)

    def test_context_and_endpoint_results_are_stable(self) -> None:
        context = self.summary["context_alignment"]
        self.assertEqual(context["bootstrap_replicates_valid"], 5000)
        self.assertLess(context["bootstrap_ci_lower"], 0)
        self.assertGreater(context["bootstrap_ci_upper"], 0)
        unifi = self.summary["unifi_endpoint"]
        self.assertEqual(unifi["overlapping_subjects"], 358)
        self.assertAlmostEqual(unifi["raw_agreement"], 0.9441340782122905)
        self.assertAlmostEqual(unifi["cohen_kappa"], 0.771712791735748)

    def test_main_and_supplementary_figure_inventory(self) -> None:
        main_names = [
            "figure_1_multicohort_design_and_evidence_geometry",
            "figure_2_strict_external_evidence",
            "figure_3_complete_transportability_atlas",
            "figure_4_context_alignment_and_endpoint_drift",
            "figure_5_program_architecture_and_redundancy",
            "figure_6_multilevel_robustness_diagnostics",
        ]
        supp_names = [
            "figure_s1_all_cell_forest",
            "figure_s2_leave_one_gene_out_heatmap",
            "figure_s3_task_specific_score_correlations",
            "figure_s4_complete_robustness_matrices",
        ]
        for stem in main_names:
            for suffix in (".png",):
                self.assertTrue((ROOT / "manuscript" / "figures" / f"{stem}{suffix}").is_file())
            with Image.open(ROOT / "manuscript" / "figures" / f"{stem}.png") as image:
                self.assertGreaterEqual(image.width, 2500)
                self.assertGreaterEqual(image.height, 1800)
                self.assertGreaterEqual(min(image.info.get("dpi", (0, 0))), 295)
        for stem in supp_names:
            for suffix in (".png",):
                self.assertTrue((ROOT / "manuscript" / "supplement" / "figures" / f"{stem}{suffix}").is_file())

    def test_source_data_are_manifested(self) -> None:
        source_dir = ROOT / "manuscript" / "figure_source_data" / "phase3"
        self.assertGreaterEqual(len([path for path in source_dir.iterdir() if path.is_file()]), 20)
        self.assertTrue((source_dir / "phase3_cell_robustness.csv").is_file())
        self.assertTrue((source_dir / "phase3_summary.json").is_file())


if __name__ == "__main__":
    unittest.main(verbosity=2)
