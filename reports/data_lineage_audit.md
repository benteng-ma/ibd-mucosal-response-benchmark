# Phase 0 data-lineage audit

## Independence rule

The trial/cohort is the external-evidence unit and the subject is the sampling unit. Accessions, tissue segments, biopsies, time points, sorted fractions, arrays and sequencing libraries are never counted as independent patients or trials.

## Confirmed high-impact relationships

1. **GSE12251 ↔ GSE23597 (ACT1): exact reupload/subset.** All 23 GSE12251 CEL files have identical uncompressed SHA-256 hashes in GSE23597. GSE12251 represents 22 unique subjects because P13 has a technical duplicate. Those 22 subjects occur in GSE23597. The accessions must never be split across training/evaluation or counted twice.
2. **ACT1 label mismatch on the same baseline files.** GSE12251 uses the stricter `WK8RSPHM` mucosal/histologic-healing label; GSE23597 uses week-8 clinical response. The same baseline profile can therefore carry different labels. Labels are endpoint-specific, not interchangeable truth values.
3. **GSE14580 ⊂ GSE16879.** All 30 GSE14580 sample accessions occur in GSE16879. GSE14580 is a subset/republication, not an independent cohort.
4. **GSE112366/GSE207022/GSE207465 are UNITI-2 family records.** GSE112366 ileum and GSE207022 rectum share 115 public patient IDs; GSE207465 is blood with remapped identifiers but belongs to the same trial. They can supply different tissue/endpoint views but only one independent trial unit.
5. **GSE73661 contains two clinical cohorts.** The vedolizumab GEMINI arm and the infliximab comparator are separate patient cohorts within one accession and are tracked separately. Conversely, repeated weeks within the GEMINI arm are one longitudinal cohort.
6. **GSE100833 includes tissue and blood from a shared clinical program.** Its multiple samples do not create independent patient sets, and missing individual outcomes prevent evaluation.
7. **GSE171770 and GSE234736 contain repeated tissues/fractions/time points.** Counts collapse to subjects; GSE234736 has 96 libraries from 30 donors, not 96 patients.

## Unresolved overlap

The GSE73661 infliximab comparator may reuse subjects from older Leuven biobanks represented in GSE14580/GSE16879. Public identifiers and platforms cannot resolve this. The primary benchmark must either exclude that pair from simultaneous independence claims or report a sensitivity analysis in which all 23 GSE73661 infliximab subjects are removed. This creates the 836-subject conservative lower bound.

## Enforcement

- Known related accessions receive `REUSED_ORIGINAL_DATA` or `SAME_TRIAL_OR_REUPLOAD_RELATED` in the frozen grid.
- Related records are barred from external-validation eligibility for any signature derived from their lineage family.
- No accession-level random split, no biopsy-level split and no cross-cohort batch-corrected mega-cohort is permitted.
- Subject/file comparison tables and exact ACT1 hashes are retained in `results/lineage/`.

See `metadata/data_lineage_matrix.csv`, `metadata/potential_overlap_pairs.csv`, `metadata/trial_manifest.csv`, `metadata/subject_manifest.csv` and `metadata/sample_manifest.csv`.

