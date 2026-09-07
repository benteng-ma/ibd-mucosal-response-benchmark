# Reproducing the frozen study

## 1. Preserve the reference snapshot

Use a separate checkout for reanalysis: the original scripts intentionally write derived outputs into their documented `results/`, `metadata/` and `manuscript/` paths. The initial release checksums describe the shipped files, not files overwritten during reanalysis. No original public source data are distributed here.

The original study used Python 3.13.11, NumPy 2.5.1, pandas 2.3.3, SciPy 1.16.3, scikit-learn 1.8.0, matplotlib 3.10.8, seaborn 0.13.2 and Pillow 12.0.0. The core Phase 1 implementation uses only the standard library. The supplied requirements capture this tested environment; identical figure pixels also require the original Arial font and rendering environment.

## 2. Obtain the Phase 1 inputs

`provenance/required_public_inputs.json` records the exact local destinations, sizes, SHA256 values, accession landing pages and available download URLs. Downloading is optional and explicit; running a snapshot test never downloads data.

```text
python scripts/prepare_public_inputs.py --download-geo
```

This retrieves the four GEO family SOFT files for GSE92415, GSE212849, GSE112366 and GSE206285. These repository-normalized files include the platform annotation, public sample metadata and expression values required by the frozen extraction rules.

From https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-7604 download:

- `E-MTAB-7604.processed.1.zip` to `data/raw/representative_expression/E-MTAB-7604.processed.1.zip`.
- `E-MTAB-7604.sdrf.txt` to `data/raw/arrayexpress_metadata/E-MTAB-7604_E-MTAB-7604.sdrf.txt`.

The different local SDRF name is intentional and is expected by the original parser. If the accession's file links move, locate these exact filenames through the official record and verify their hashes. Do not substitute reprocessed counts, map missing subjects, or ignore changed hashes. The frozen archive contains 43 count files, whereas the SDRF has 44 subjects; the unavailable sample remains excluded by the original rule.

```text
python scripts/prepare_public_inputs.py --check
python scripts/reproduce.py --stage phase1
```

The Phase 1 wrapper checks raw-file and frozen-input hashes, parses the four GEO SOFT files locally, extracts frozen program genes, builds the six subject-level task manifests, runs the original benchmark and verifies the recreated aggregate output against the shipped reference. It does not require the private Git repository.

## 3. Recompute exploratory analyses and figures

```text
python scripts/reproduce.py --stage phase3
python scripts/reproduce.py --stage figures
python -m unittest discover -s tests -p "test_*.py"
```

Or use `python scripts/reproduce.py --stage full` after obtaining the inputs. The Phase 3 program first verifies the 46 frozen z-score AUROCs and then performs the original scaling, subject/gene omission, orientation enumeration, bootstrap, context and endpoint analyses. Confirmatory P values remain restricted to the three original cells. A change of manuscript title or journal does not alter this family.

Figure 4A has the documented color-legend correction. No biological or statistical result was changed. The source-data export is local and can recreate ignored subject-score tables; do not force-add those files to Git. The historical supplemental S14 table is retained unchanged as the analysis record; `FILE_MANIFEST_SHA256.csv` lists the actual public release payload.

## 4. Phase 0 provenance

Frozen study-level Phase 0 tables and original audit code are included so lineage and eligibility decisions can be inspected independently of outcome performance. The original R dependency lock is `renv.lock`. Install R/renv, restore with `Rscript scripts/bootstrap_phase0.R`, then run `Rscript tests/testthat.R` for the available audit checks.

A full historical Phase 0 reconstruction needs more public metadata and source files than the five executed clinical cohorts. These are listed in `metadata/file_manifest.csv` and `metadata/dataset_manifest.csv`, with the original search strategy and citation inventory in `literature/`. `scripts/query_geo_metadata.py` retrieves and parses GEO accession metadata. E-MTAB-7604/E-MTAB-7845 IDF/SDRF and BioStudies records must be obtained from their official records. Coverage additionally requires the source platforms and GSE234736 supplementary table described in `scripts/audit_signature_coverage.py`; the ACT1 raw-CEL comparison requires its source files described in `scripts/audit_act1_overlap.py`.

Do not run `build_phase0_metadata.py` with an incomplete source collection and treat the output as the original audit. The public checkout is designed to reproduce Phase 1/3 using the frozen Phase 0 definitions, while preserving the broader historical audit trail. It is not claimed to include a one-command re-download of every historical Phase 0 source. Database search hit counts and mutable repository metadata may change after the 2026-08-31 cutoff.

## Validation scope

Snapshot tests check aggregate invariants and synthetic metric behavior, not a full raw-data reanalysis. Subject-level validation tests require regenerated local inputs and explicitly skip when those files are absent. The public preparation audit records separately which code was executed, which outputs matched, and whether network acquisition was tested. Historical source commit strings alone do not constitute a publicly timestamped preregistration.
