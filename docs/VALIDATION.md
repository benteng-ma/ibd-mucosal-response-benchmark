# Release-candidate validation — 2026-09-08

This record describes executed checks, not an assertion that all environments or future downloads will reproduce identically.

- **Full local reanalysis:** copied the clean public candidate into a separate validation directory; supplied the six required public-source files from the existing local archive; verified all six SHA256 values; ran the candidate's `scripts/reproduce.py --stage full`. Exit status 0, about 230 seconds in the tested Windows/Python environment.
- **Frozen outputs:** every aggregate file registered in `provenance/frozen_result_sha256.json` matched the original Phase 1 and Phase 3 reference byte-for-byte. Seven frozen scientific input files remained unchanged. This includes the original three-cell confirmatory family and unchanged multiplicity results.
- **Full Python checks:** 19 tests passed; no skips after local subject-level reconstruction.
- **Aggregate-only Python checks:** 12 tests passed; one subject-dependent test class explicitly skipped because its inputs are deliberately not distributed. This skip must not be counted as seven successful subject-level tests.
- **Aggregate-only figures:** `plot_public_snapshot.py` successfully regenerated Figures 1, 3–6 and S1–S4 without subject-level files. Figure 2 remains an included publication PNG in this mode. The full reanalysis regenerated all six main and four supplemental figures.
- **Corrected Figure 4:** the full reanalysis output has identical dimensions and pixels to the corrected iScience publication PNG.
- **R audit checks:** the nine available Phase 0 test files passed; the individual-level subject-mapping test explicitly skipped. Tests used the existing original renv library through R 4.6.1. Locale warnings were emitted at startup; no test failure occurred.
- **Payload screening:** selected credential patterns, absolute workstation paths and explicit individual-record identifier columns were screened. Raw sources, individual-level derived files, manuscript/submission documents, administrative declarations and private Git history were excluded. This is a scoped screening, not a guarantee that every possible disclosure pattern has been ruled out.

Not tested: a fresh network acquisition of all six source files; a clean dependency install/renv restore on a second computer; Linux/macOS execution; complete historical Phase 0 source reconstruction. The same tested software versions and frozen input hashes are necessary reference points, not proof of universal platform equivalence.

The public payload does not include the private validation checkout or its individual-level reconstructed data. `FILE_MANIFEST_SHA256.csv` describes the exact distributed snapshot. Reanalysis intentionally writes new derived outputs into a separate checkout.
