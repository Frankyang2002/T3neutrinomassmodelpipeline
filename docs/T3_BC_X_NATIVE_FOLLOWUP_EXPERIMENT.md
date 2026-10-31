# B/C(-1): native X-only followup of the successful Y-only Matchete probe

**Status: experimental, not physical matching.** This is a read-only standalone fork; existing numerical programs, output studies, and the original Matchete source are not modified.

## Current verified inputs

User's Windows log (2026-10-31): **both B(-1) and C(-1)** returned
`IndependentCanonicalFullEFTCheckPassedEOMEquivalenceUnproved` followed by
`YOnlyIndependentCanonicalRawEFT2NotCertified`. Thus the **independently
reconstructed** full EFT1 passed ordinary native Matchete checks and produced
raw Y-only EFT2 results for both models. This is not a proof that the
independently reconstructed scalar quartics are exactly equivalent to the
F-first serialized EOMEFT. The original saved EOMEFT was demonstrated to fail
Hermiticity in four self-quartic coupling sectors.

## Implemented step

New complete Wolfram file:
`Lagrangian/RunT3BCXYNativeFollowup.experimental.wl`

1. Repeats the successful independently reconstructed full Y-only
   `RunT3Matching` at tree and one loop, preserving its raw operator exports.
2. Exports the exact Y-only `treeResult["EOMEFT"]` and
   `loopResult["EOMEFT"]` as separate InputForm artifacts.
3. Requires Y-only tree EOMEFT to have light X but no Y/F; checks the full
   `CheckLagrangian` and all detailed flags **in the X-light registry**.
4. Captures `FreeLag[NewScalarX]` while X is light, removes X, defines X as
   heavy (`MSX`), rebuilds `FreeLag[NewScalarX]`, replaces only the free
   Lagrangian, and checks the full resulting X-heavy EFT2. No operator deletion,
   Hermitian averaging, or fallback to a different EFT representation.
5. With `--check-x`, stops after the verified X-heavy EFT2 check.
6. With `--match-x`, requests tree and one-loop native Matchete X-only matching
   on **the Y-tree EOMEFT**. Exports X tree, X tree+loop, their difference, and
   the raw LLHH part; strictly requires X/Y/F to be absent afterward.
7. Both paths emit an `_x_followup.json` report with explicit success/failure
   status, source checks and all physical certification flags **false**.

The numerical sequence under study is

`UV -> F tree matching -> independently reconstructed EFT1 -> native Y tree
 -> [candidate EFT2] -> native X tree+loop -> SM + LLHH`.

This is **not** the complete sequential one-loop chain for four separate
reasons:

- Original F-first EOMEFT and independently reconstructed EFT1 are not proven
  equivalent; the original is not Hermitian as serialized.
- The X input is **Y-tree EOMEFT only**, not Y-tree plus inherited Y one-loop
  matching. The Y loop source is exported for future consistent inclusion.
- The Matchete parameter input is the saved mass-basis probe point, not the
  running values at the separately computed `continuous_v1` thresholds.
- Renormalized matching requires counterterms, the heavy/light hard-region
  subtraction and coefficient normalization to the numerical Wilson basis;
  these are not supplied by two successful `RunT3Matching` calls.

The X-heavy EOMEFT may fail Hermiticity or canonical normalization. In that
case the runner stops with an informative status: this is a useful failure,
not permission to use `MatchedEFT` or delete quartics.

## Intended project paths (ZIP)

- `Lagrangian/RunT3BCXYNativeFollowup.experimental.wl` — complete independent
  Wolfram runner (not a patch)
- `scripts/RunT3BCXYNativeFollowup.py` — strict CLI, source SHA fingerprints,
  status gates and provenance
- `scripts/CollectT3BCXYNativeEvidence.py` — read-only source-hashed ZIP
  collector of Y/X raw expressions for subsequent finite matching
- `tests/test_t3_bc_xy_native_followup.py` — standalone static/contract tests
- `docs/T3_BC_X_NATIVE_FOLLOWUP_EXPERIMENT.md`
- `docs/T3_BC_X_NATIVE_FOLLOWUP_MILESTONE.md`

Outputs under:
`output/general_t3_mass_eigen_matchete/eft1_scalar_xy_native_followup_v1/`.

The X matching report for model B is
`B_m1_final_xy_native_match_x_x_followup.json`; for C, replace `B_m1` with
`C_m1`. The successful diagnostic status would be
`XFollowupYTreeOnlyXRawMatchingNotCertified`.

## PowerShell instructions

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$env:USERPROFILE\Downloads\t3_bc_xy_native_followup_v1.zip" -DestinationPath . -Force
python -m unittest discover -s tests -p "test_t3_bc_xy_native_followup.py" -v
python -u -X utf8 -m scripts.RunT3BCXYNativeFollowup --all --preflight-only
python -u -X utf8 -m scripts.RunT3BCXYNativeFollowup --all --check-x
# Run only after BOTH check-x statuses pass:
python -u -X utf8 -m scripts.RunT3BCXYNativeFollowup --all --match-x
python -u -X utf8 -m scripts.CollectT3BCXYNativeEvidence --all
```

Evidence path:
`output/audits/t3_bc_xy_native_followup/xy_native_matching_evidence.zip`.
Upload that ZIP for direct comparison and finite-coefficient extraction. The
collector will refuse to produce a ZIP if the expected prior reports are
missing; it never declares physical matching complete.

To avoid overwriting prior runs, the native runner refuses existing output
unless `--overwrite` is specified, and the collector refuses overwriting an
existing evidence archive.

## Verification and limitations

Static tests validate source anchors, gate logic, output contracts, and basic
Wolfram delimiter balance. They **do not** execute Matchete or prove the
physics. The full matched EFT1/EFT2 expressions and numeric Wilson boundary
must be reconciled before promotion to physical results.
