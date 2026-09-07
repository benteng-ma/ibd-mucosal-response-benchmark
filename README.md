# IBD mucosal response benchmark

Reproducibility code for **Trial aware evaluation of mucosal response programs reveals limited external support in inflammatory bowel disease**.

This repository evaluates finite published mucosal programs without feature selection, model refitting, outcome-driven thresholds, cross-cohort expression merging, or cross-study batch correction. Clinical trial/cohort lineage, endpoint definitions and evaluation roles are explicit.

## Frozen study

- Scientific cutoff and analysis seed anchor: 2026-08-31 / 20260831.
- Five independent clinical units, 636 unique subjects, six tasks, nine programs, 46 executed cells.
- Three class- and endpoint-family-matched confirmatory cells: two bootstrap intervals above AUROC 0.50; zero Holm-significant results.
- Nineteen directionally supported and 27 uncertain descriptive cells. Eight unavailable raw-only evaluations are not negative results.
- Phase 3 sensitivity analyses are exploratory and do not create further confirmatory tests or establish clinical utility or causality.

## Install

Python 3.13 was used. From the repository root, create an isolated environment, then:

```text
python -m pip install -r requirements.txt
python -m unittest discover -s tests -p "test_*.py"
```

The dependency-free Phase 1 scoring implementation is in `scripts/run_phase1_benchmark.py`. Scientific validation and Phase 3 require the packages above. Tests that require locally regenerated subject-level files explicitly skip in an aggregate-only checkout; the independent aggregate, synthetic metric and integrity tests still run. Skipped tests are not evidence of a successful full reanalysis.

Executed release-candidate checks and untested conditions are recorded in `docs/VALIDATION.md`.

## Reproduction modes

1. **Audit the shipped snapshot, without downloads:** run the tests and `python scripts/verify_release.py`. Aggregate outputs, frozen definitions and publication PNGs are included.
2. **Reproduce Figures 1, 3–6 and S1–S4 from aggregate outputs:** run `python scripts/plot_public_snapshot.py`. Figure 2 requires local subject scores for its distribution panel and is not regenerated in this mode; its publication PNG and aggregate ROC/null data are included. Rendering may differ with fonts/package versions.
3. **Recompute the benchmark and all sensitivity analyses:** follow `docs/REPRODUCING.md`. Obtain the six hash-identified public source files, regenerate local subject metadata, then run `python scripts/reproduce.py --stage full` in a separate checkout. This writes new derived files. It does not access the authors' original workspace or private Git history.

Use `python scripts/prepare_public_inputs.py --check` to check required files, or `--download-geo` to retrieve the four GEO files. The two BioStudies files should be downloaded from the official accession page as explained in the reproduction guide. Source changes cause a checksum failure, not automatic substitution.

## Contents

- `scripts/`, `R/`, `config/`: analysis, audit, visualization and frozen rules.
- `literature/`, `metadata/`: program definitions and study/task-level provenance, not individual clinical records.
- `results/`: frozen aggregate performance and sensitivity results.
- `manuscript/figures/`, `manuscript/supplement/`: publication PNGs and machine-readable Tables S1–S14, not manuscript or submission files.
- `manuscript/figure_source_data/phase3/`: aggregate figure sources. Subject-level sources named in the historical S14 manifest are deliberately omitted and are regenerated locally.
- `docs/`: input acquisition, reporting corrections and reproduction scope.

## Release history and reporting corrections

This is a clean export, not a copy of the private working repository's Git history. Historical analysis commit identifiers are documented in `provenance/source_commits.json`; they are audit identifiers and may not resolve in this public repository. Frozen input hashes are independently checked without requiring those commits.

The published statistics use a **two-sided** permutation test based on `abs(AUROC - 0.50)` with 10,000 permutations and the add-one correction. Earlier manuscript wording incorrectly called it one-sided; the test and its P values were not changed. The plotting code now explicitly maps Figure 4A legend labels to colors and uses seed 20260905 for that figure's horizontal jitter. See `docs/REPORTING_CORRECTIONS.md`.

## Data protection, licensing and citation

Public availability of source data is not permission to redistribute everything. Third-party expression files, subject-level clinical/outcome/expression tables, downloaded articles, personal administrative material, manuscript drafts and account credentials are excluded. Inputs and regenerated individual-level products are ignored by Git. Do not force-add them after running the pipeline.

The existing MIT license applies to original software; it does not grant new rights over third-party data or publications. Software credit in `CITATION.cff` is separate from article authorship. Cite the exact software release/version DOI once available; do not invent a DOI or cite a moving branch as an immutable version.

Code archiving is not journal acceptance, ethics approval, prospective validation or a change in the frozen scientific conclusion.
