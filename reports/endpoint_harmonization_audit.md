# Phase 0 endpoint-harmonization audit

## Harmonization rule

Endpoints are harmonized to named families only when the original clinical definition and assessment time remain traceable. “Response”, “remission”, “mucosal healing”, “endoscopic response”, “withdrawal relapse” and “continued-treatment flare” are not collapsed into a single binary outcome.

## Frozen endpoint families

1. `INDUCTION_CLINICAL_RESPONSE`
2. `INDUCTION_CLINICAL_REMISSION`
3. `MUCOSAL_HEALING`
4. `ENDOSCOPIC_RESPONSE`
5. `ENDOSCOPIC_REMISSION`
6. `PRIMARY_NON_RESPONSE`
7. `RELAPSE_AFTER_WITHDRAWAL`
8. `FLARE_ON_CONTINUED_THERAPY`
9. `ON_TREATMENT_MOLECULAR_RESPONSE`
10. `HISTOLOGIC_RESPONSE`
11. `HISTOLOGIC_REMISSION`
12. `STEROID_FREE_CLINICAL_REMISSION`
13. `COMBINED_CLINICAL_ENDOSCOPIC_RESPONSE`
14. `COMBINED_CLINICAL_ENDOSCOPIC_REMISSION`
15. `OUTCOME_UNAVAILABLE`

## Critical mappings

- ACT1 GSE12251 mucosal/histologic healing and GSE23597 week-8 clinical response remain separate despite identical baseline files.
- GSE73661 vedolizumab mucosal healing at weeks 6, 12 and 52 remains three time-specific tasks sharing one subject pool.
- E-MTAB-7604 uses disease-specific endoscopic-remission timing: UC week 10 and CD week 24. Because SDRF lacks per-sample UC/CD subtype, sample-level timing cannot be assigned without an additional verified key.
- Ustekinumab clinical response, clinical remission and mucosal healing are separate task rows; UNITI-2 tissue accessions remain one trial.
- BIOSTOP withdrawal relapse and continued-treatment flare are separate estimands. Neither is an induction-response proxy.

## Permitted comparisons

- Primary external validation: same therapy, same endpoint family, compatible disease/compartment, independent lineage.
- Transport stress tests: endpoint, disease or therapy mismatch must be labelled explicitly in the grid.
- Broadly comparable endpoints may support sensitivity analyses, never silent pooling or label substitution.

The complete definitions and row-level mappings are in `metadata/endpoint_dictionary.csv` and `metadata/endpoint_mapping.csv`.

