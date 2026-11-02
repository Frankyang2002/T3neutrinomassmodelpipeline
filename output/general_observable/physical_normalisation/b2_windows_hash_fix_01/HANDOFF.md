# B2 Windows Git-blob provenance fix (isolated correction)

- Reference commit: `0b210f18f6a2cfb10fff03059e2ec4a73d2195a4`.
- Corrected B2-only source validation: raw file Git blob comparison still preferred; if raw file bytes differ, perform independent `git hash-object --path=...` canonical comparison plus exact HEAD tree blob, index blob, HEAD revision, and zero Git diff. A CRLF checkout can therefore match canonical LF source without relaxing source-drift rejection.
- No production, CG audit, RGE, benchmark, EFT, A or B1 files changed.
- Synthetic local regression: `15 passed` in Python 3.13.5 using Git 2.47.3 and pytest 9.0.2. Includes CRLF-positive, actual-drift-negative, HEAD-drift-negative and missing-source-negative cases.
- Actual user-local source and Matchete/RGBeta execution were **not** available here; user rerun required.
- The physical source-to-source tensor/action dictionary remains `physical_source_certified=false` for all 16.
- This fix only promotes *source provenance*, not native SU(2) numeric RGBeta tensor equivalence or physical couplings.
- Run the pre-existing B2 script to a **new directory** `output/general_observable/physical_normalisation/b2_local_02`.
- Inspect source `verification` entries for the two previously blocked `.wl` files; expected `path_specific_git_clean_filter` *only if* the cause is normalized line endings. Otherwise remain blocked.
