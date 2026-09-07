# Reporting corrections incorporated into this release

These corrections reconcile the manuscript/figure presentation with the existing data and code; none changes program definitions, outcome labels, AUROCs, confidence intervals or the confirmatory multiplicity family.

1. **Permutation sidedness:** the original Phase 1 code tests absolute deviation of AUROC from 0.50, hence two-sided. An earlier manuscript called this one-sided. Label shuffling is within task across subjects, preserving class counts; stratification by response is used for bootstrap sampling. Stored nulls reconstruct raw P values 0.0425957404, 0.1791820818 and 0.0430956904, and Holm values 0.1277872213, 0.1791820818 and 0.1277872213, respectively for Arijs5/PROgECT, Anti-TNF7/PURSUIT and Cross-UC4/PROgECT.
2. **Figure 4A:** previous handles were taken in encountered group order but labels were assigned in a different order. Explicit handles now map strict external to navy, within-class transport to blue and cross-therapy to orange. Seed 20260905 fixes horizontal jitter for this figure; vertical AUROC values and grouping are unchanged.
3. **Unavailable evaluations:** eight E-MTAB-7845 cells comprise two within-class and six cross-therapy planned transfers, not eight class-matched validations. They remain unavailable, never negative performance results.
4. **Treatment matching:** matching uses treatment class, not necessarily the same molecule. Anti-TNF programs derived with infliximab can be class-matched to golimumab tasks; this does not establish clinical interchangeability.
5. **Legends and citations:** the iScience manuscript corrects malformed S11/S12 citations and describes the actual displays in Figures 3–6 and S3/S4. Figure 6D shows only directionally supported multigene cells, while complete omission results remain in S2/S7. S4 displays orientation percentile, bootstrap fraction, maximum scaling change and gene-omission direction fraction; it does not display four alternative AUROC matrices.

The central conclusion remains limited confirmatory support despite numerical stability under the examined perturbations. No causal explanation for incomplete transfer, clinical utility or numerical inflation caused by duplicate accessions was established.

Public-export housekeeping: the workstation-specific Python executable path in Table S5 is replaced with a descriptive placeholder. Software versions and all analytical values are unchanged. The original manuscript package is retained separately, not edited by this export.
