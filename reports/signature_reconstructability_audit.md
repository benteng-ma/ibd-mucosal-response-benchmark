# Phase 0 published-signature reconstructability audit

## Main-grid inventory

Nine finite published programs pass the Phase 0 reconstructability rule across anti-TNF, vedolizumab and ustekinumab sources:

1. `SIG_ATNF_ARJIS5`: TNFRSF11B, STC1, PTGS2, IL13RA2, IL11; nonresponder-high.
2. `SIG_ATNF_OSM1`: OSM; nonresponder-high.
3. `SIG_ATNF_7GENE`: WNK2/OCRL/ASB7 responder-high; PCBP3/AMPD2/FAM155A/IL13RA2 nonresponder-high.
4. `SIG_VDZ4`: PIWIL1, MAATS1, RGS13, DCHS2; remitter-high.
5. `SIG_VDZ_ANDO7`: ARMC2, BCL11B, CDH23, CRB2, MEGF6, NAPSB responder-high; BCAP31 nonresponder-high.
6. `SIG_UST4_COEF`: HSD3B1, MUC4, CFI, CCL11 with published coefficients. The article's model string prints `CF1`; its own univariate predictor table identifies `CFI`, so this is recorded as a transparent I/1 typographic correction. The intercept is unavailable, preventing reconstruction of an original calibrated probability.
7. `SIG_UST_HUB4`: MUC1, DUOX2, LCN2, PDZK1IP1; lower in responders.
8. `SIG_UST_MH2`: LCN2 and KDM5D; nonresponder-high.
9. `SIG_CROSS_UC4`: IGFBP5, SELE, STC1, VNN2; nonresponder-high.

## Secondary or excluded programs

- `SIG_ATNF_TREM1_BLOOD` is reconstructable but blood-based and excluded from the primary mucosal grid.
- `SIG_MIN6` membership is known (G0S2, S100A9, SELE, CHI3L1, MMP1, CXCL13), but exact direction/transform/cut point was not recovered from an accessible verified source; it remains outside the main grid.
- Six legacy panels, a 21-gene vedolizumab panel, GIMATS and BIOSTOP baseline candidates remain queued. No genes were inferred from heatmaps.

## Coverage audit

- Direct platform-symbol catalogs were parsed for GPL570, GPL6244 and GPL13158.
- E-MTAB-7604 processed counts contain direct gene symbols.
- GSE234736 and GSE171770 use Ensembl-keyed RNA-seq counts and require a versioned crosswalk.
- E-MTAB-7845 is raw-FASTQ-only in the repository; platform-level coverage can be anticipated but dataset-level processed availability remains a high-cost Phase 1 prerequisite.
- The frozen matrix contains 242 signature-dataset pairs; an 80% directional-gene coverage rule is required for primary eligibility.

Reconstructability means a future deterministic score can be specified without new feature selection. It does not mean the original calibration, threshold or predictive performance has been reproduced.

See `literature/signature_inventory.csv`, `literature/signature_gene_table.csv`, `literature/signature_extraction_queue.csv` and `metadata/signature_coverage_matrix.csv`.

