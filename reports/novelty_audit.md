# Phase 0 novelty audit

Audit cutoff: **2026-08-31**  
Status: **no exact duplicate located; near-duplicate risk HIGH**

## Locked novelty claim

The defensible novelty is not “another IBD response signature.” It is a lineage-aware benchmark of **already published, finite mucosal treatment-response programs** across independently identified trials/cohorts, with endpoint-family labels kept explicit and with both successful and negative transport treated as reportable. The design unit is the clinical trial/cohort and subject, not the repository accession, biopsy, time point or expression file.

No located study satisfied the full conjunction of: at least three therapy classes; systematic published-signature inventory; formal trial/subject/sample/accession lineage; endpoint harmonization; original external evaluation; same- and cross-therapy grid; leave-study-out logic; and explicit negative transportability.

## Ten closest studies

1. **2026 Biology Direct multi-therapy fibroblast analysis** — closest scope collision because it spans four therapies and seven GEO records. Its unit accounting is not an audit-grade clinical lineage: samples are described as patients, ACT1 reuploads are not collapsed, and longitudinal/sorted libraries and same-trial UST accessions are treated as if independent. It does not benchmark a frozen inventory of published signatures across harmonized endpoints.
2. **2026 Frontiers anti-TNF inflammatory-remodelling state** — close on pretreatment mucosal state and reproducibility claims, but only anti-TNF. It labels GSE12251, GSE16879 and GSE23597 as independent cohorts even though GSE12251 is an exact ACT1 reupload within GSE23597 and GSE14580 is a subset of GSE16879.
3. **2026 JCI Insight IL11+ fibroblast study** — mechanistically close and reuses GSE16879/GSE12251/GSE23597/GSE212849, but integrates duplicated ACT1 data and uses a random split rather than an independent trial benchmark.
4. **Ando et al. 2026** — a seven-gene vedolizumab mucosal-healing model with external use of GSE73661. It is single-therapy, has ten discovery subjects, and the public GSE282580 record lacks sample-to-week-54 outcome mapping.
5. **2026 medRxiv vedolizumab T-cell interferon study** — five-cohort and leave-one-cohort-out work, but its primary compartment is blood/T cells, two cohorts are author-held, and it does not provide a bulk-mucosa cross-endpoint/cross-therapy benchmark.
6. **Gaujoux et al. 2019** — the nearest systematic anti-TNF signature comparison and cell-centred meta-analysis, but it is anti-TNF only and lacks a formal cross-therapy endpoint/lineage grid.
7. **Shanthamallu et al. 2021 seven-gene anti-TNF signature** — uses leave-one-dataset-out logic, but the included accessions are not independent clinical units and no clean independent holdout is created after lineage collapse.
8. **Verstockt et al. 2020 vedolizumab four-gene program** — provides discovery and validation and consumes the two principal public vedolizumab mucosal cohorts E-MTAB-7845 and GSE73661. It is a signature source, not a multi-therapy benchmark.
9. **2021 Scientific Reports UC four-gene cross-biologic program** — spans anti-TNF and vedolizumab and is therefore a direct near neighbour, but it lacks three therapy classes, formal endpoint harmonization, and a prospective trial/subject/sample lineage audit.
10. **2021 MIN score** — crosses anti-TNF and vedolizumab and reports validation, but the accessible record did not yield a fully reconstructable directional/calibrated implementation and it does not perform the locked lineage-aware benchmark.

BIOSTOP is a clinically important adjacent study rather than an exact collision: relapse after withdrawal and flare on continued treatment are distinct endpoint families, not induction response or mucosal healing.

## Collision interpretation

- Exact duplication: **not found**.
- Component-level novelty: **limited**; signature discovery, multi-cohort analysis, individual cross-biologic markers and leave-one-dataset-out studies already exist.
- Design-level novelty: **retained**, conditional on preserving the lineage and endpoint constraints frozen here.
- Novelty risk: **HIGH**, because several 2026 publications are close and because a future manuscript could easily overstate the contribution as a new universal predictor.

The full machine-readable comparison is `literature/novelty_matrix.csv`; the executed source-query log is `literature/search_log.csv`.

