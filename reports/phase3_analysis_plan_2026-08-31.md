# Phase 3 exploratory robustness and figure-density plan

Frozen: 2026-08-31  
Source commit: `66d31f95f8cbb7eafaaf0e32e7aac8132f049990`  
Status: frozen before inspecting any Phase 3 result

## Rationale

The Phase 2 manuscript has a defensible audit design but only four single-panel main figures. Phase 3 increases analytical and visual density without changing the Phase 1 confirmatory family, selecting new genes, optimizing thresholds, merging expression datasets, or hiding negative results.

## Fixed analyses

1. Recalculate every executable program-task cell under the frozen z-score implementation and two non-outcome-informed alternative scalings: median/MAD and centered percentile ranks.
2. Quantify subject influence by leave-one-subject-out AUROC for all 46 cells.
3. Quantify component dependence by leaving out each measured gene from every multi-gene program-cell and renormalizing remaining weights.
4. Enumerate the full same-gene sign-orientation space for every cell and locate the published orientation within that descriptive null space.
5. Estimate stratified-bootstrap direction stability for all cells without creating a new multiplicity family.
6. Summarize context alignment from therapy, endpoint, and disease matches, including a two-way program/task bootstrap and within-program same- versus cross-therapy contrasts.
7. Map program architecture with signed gene membership, Jaccard similarity, shared-gene direction agreement, subject-score correlations, and cross-task performance-profile correlations.
8. Use the overlapping UNIFI subjects to quantify mucosal-healing versus clinical-remission endpoint concordance and paired program performance.
9. Expand the three strict external cells into score-distribution, ROC, forest, and empirical-null diagnostics while preserving the original permutation P values and Holm results.

## Frozen visual architecture

- Figure 1: study design and evidence geometry, four panels.
- Figure 2: strict external evidence, four panels.
- Figure 3: complete transportability atlas, four panels.
- Figure 4: context alignment and UNIFI endpoint drift, four panels.
- Figure 5: program architecture and redundancy, four panels.
- Figure 6: implementation, gene, subject, bootstrap, and orientation robustness, six panels.
- Supplementary Figures S1-S4: full-cell forest, detailed gene influence, task-specific score-correlation matrices, and complete orientation/robustness diagnostics.

All panels will use letter labels and external figure legends rather than narrative titles inside the plotting area.

## Interpretation lock

Phase 3 is exploratory. A stable AUROC under alternative scaling or leave-one-out analysis can demonstrate implementation robustness but cannot supply missing multiplicity-controlled confirmation. Orientation enumeration is not a genome-wide random-gene null. Context summaries are dependent because programs and tasks repeat. The manuscript must still reject universal-biomarker, calibrated-model, treatment-selection, formal-ranking, clinical-utility, and causal-mechanism claims.

## Pre-result implementation clarification

The first execution stopped before producing Phase 3 results because PIWIL1 had a zero median absolute deviation in one dataset. To retain the gene and every frozen cell, the robust-scaling algorithm now uses a deterministic dispersion fallback: MAD × 1.4826, then IQR/1.349 if MAD is zero, then population SD only if both robust measures are zero. This rule is applied gene-wise without outcome labels and was recorded before a successful result run.
