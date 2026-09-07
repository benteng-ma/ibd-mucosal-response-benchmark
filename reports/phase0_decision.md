# Phase 0 executive decision

## Decision

**GO_CROSS_THERAPY_ENDPOINT_BENCHMARK**

This is a go to prepare a strictly lineage-aware Phase 1 benchmark, not a claim that any program predicts treatment response. Phase 1 has not been executed.

## Gate evidence

- Three major therapy families have public baseline mucosal expression with sample-level outcomes: anti-TNF, vedolizumab and ustekinumab.
- Six anti-TNF, three vedolizumab and two ustekinumab independent trial/cohort units are represented at the manifest level; secondary/non-executable units remain labelled.
- Conservative public active-baseline subject lower bound is 836, above the locked 300 threshold.
- Each major therapy has at least one task with both outcome arms at or above 15 subjects.
- Nine finite directional/weighted main-grid signatures come from nine source papers and span three therapy families.
- The frozen grid contains 54 independent transportability opportunities across six independent trials/cohorts, including 35 cross-therapy stress-test cells.
- Three strict same-therapy/same-endpoint external opportunities exist across two independent trials. This passes the “at least one” strict opportunity gate but is an important evidence-width limitation.
- Multiple endpoint families can be harmonized without collapsing their clinical meaning.
- Known reuploads, subsets, same-trial tissues and repeated samples are explicitly blocked from independent validation.
- No Phase 0 decision requires cross-cohort merging or batch-corrected mega-analysis.

## Recommended scope

Evaluate only frozen published programs. Begin with independent same-therapy/same-endpoint tasks, then run pre-labelled within-therapy endpoint/disease transport and cross-therapy stress tests. Preserve trial/cohort boundaries, summarize within study first, and treat non-transport as an expected scientific result.

## Prohibited claims

- Universal or pan-biologic biomarker.
- Mechanistic or causal treatment-resistance program.
- Independent validation from a reupload, subset or another tissue from the same trial.
- Equivalence of induction response, mucosal healing, remission, withdrawal relapse or continued-treatment flare.
- Clinical utility or deployability from retrospective public cohorts.

## Unresolved high risks

Possible Leuven subject overlap, raw-only E-MTAB-7845 processing, RNA-seq identifier mapping, narrow strict external trial width, and absent public BIOSTOP mapping must remain visible in Phase 1.

The machine-readable decision is `reports/phase0_decision.json`.

