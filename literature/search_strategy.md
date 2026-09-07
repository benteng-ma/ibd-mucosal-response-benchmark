# Locked literature search strategy

Audit cutoff: **2026-08-31**. This search supports Phase 0 novelty, data-lineage and signature-reconstructability auditing only.

## Sources

- PubMed E-utilities
- Europe PMC REST API and full-text XML when available
- Crossref Works API
- GEO and BioStudies/ArrayExpress repository records
- Trial/publication pages reached from verified identifiers

## Query families

The executable log in `search_log.csv` contains 19 locked query strings run against all three bibliographic sources (57 source-query records). Query concepts combined IBD/UC/CD, treatment response/nonresponse, mucosal or biopsy transcriptomics, anti-TNF/infliximab/golimumab, vedolizumab, ustekinumab, biomarker/signature/model, external validation, transportability, leave-one-study-out, endpoint harmonization, data lineage, duplicate cohorts and BIOSTOP/withdrawal.

## Inclusion

1. Human IBD therapy-response studies with pretreatment or clearly labelled longitudinal transcriptomics.
2. Published finite gene programs, coefficients, directions or explicitly described gene sets.
3. Multi-cohort/multi-therapy analyses relevant to novelty collision.
4. Repository records that expose sample metadata, platforms, files or trial links.
5. BIOSTOP because treatment withdrawal/continuation creates a distinct endpoint-family test.

## Exclusion or secondary-only

- Non-human experiments without a linked human response cohort.
- Unverifiable gene lists inferred only from heatmaps or figures.
- Post-treatment-only outcome-labelled profiles for pretreatment prediction.
- Blood/PBMC/sorted-cell cohorts as primary bulk-mucosa evaluations.
- Records lacking sample-level outcome mapping as executable tasks.
- Same-trial tissues, subset accessions or exact reuploads as independent validation.

## Verification policy

Every main-grid signature requires a stable publication identifier, a finite gene list, explicit direction or coefficient, and a source location. Exact wording was not transcribed from inaccessible supplements. Unresolved panels remain in `signature_extraction_queue.csv` and are excluded from the main grid. Search counts alone are not evidence of eligibility; eligibility comes from full record/article and repository verification.

## Novelty decision rule

An exact collision requires the conjunction of: at least three therapy classes; a systematic published-signature inventory; formal trial/subject/sample lineage; endpoint harmonization; original external validation; same- and cross-therapy evaluation; leave-study-out logic; and explicit negative transportability. Partial overlap on individual components is recorded as near-duplicate risk, not treated as an exact duplicate.

