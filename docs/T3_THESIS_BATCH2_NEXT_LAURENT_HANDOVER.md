# T3 A(-4) / D(-2) — Batch 2 next stage: complete-source loop-function and Laurent audit

Source base: **read-only GitHub `main` at `c054b8f6012485cc1be2351ef4d3a903a94cb8fa`**, 10 November 2026. The source commit contains Batch 1 results, but does **not** contain the user's newly generated Batch 2 raw outputs. This stage reads those raw outputs from the **user's local Windows checkout**.

## What is implemented

This package operates on an **already successful**
`output/parallel_t3/full_general_batch2/thesis_full_general_batch2_raw_v1/{A_-4,D_-2}`.
It does not rerun Batch 1 or `RunSequentialT3Matching`.

1. Validates each local raw Batch 2 manifest, native F-only status, complete operator ledger and SHA256 source provenance against the exact Batch 1 files. It rejects changed source and missing operators.
2. Computes **numerical, source-state-dependent heavy flavour diagnostics** using the complete running 3x3 `MF`, `y1`, `y2`, and the `D_-2` `y2Conj`, bilinear and five general quartics. The two candidate flavour contractions use `MF^{-1}` and `(MF^T)^{-1}`, separately, as **nonphysical orientation diagnostics**. Nothing is called a physical Wilson coefficient.
3. Loads the **full, unselected** one-loop Matchete `eft1_after_F_one_loop_inputform.txt` in Wolfram; isolates its exact `hbar^1` factor; rejects higher-loop/tree contamination or remaining heavy fermions; exports the full symbolic hbar coefficient as WXF.
4. Invokes the existing `EvaluateLoopFunctions` on **all** the matching terms; exports the resulting full expression; tries an explicit Laurent decomposition into `epsilon^{-1}` residue and `epsilon^0` regular piece. It checks exact reconstruction of the retained terms and preserves WXF artifacts. If residual loop functions, multiple dimensional regulators or unresolved higher poles remain, it marks the model **Blocked**, without inventing missing physics.
5. Writes checksummed per-model audit reports, logs, and a run summary. No old outputs are overwritten.

**Limitation:** `epsilon^0` of the Matchete full raw one-loop correction is *not* the finite `MSbar` Wilson coefficient at `MF`. This package has no UV-hard minus EFT1-soft subtraction, local operator counterterms, physical native-to-Matchete parameter/CG dictionary, or full EFT1 Wilson RGE. These require further implementation. No actual Mathematica evaluation was possible in the authoring environment; only Python interface logic and offline tests were executed.

## Complete files added (nothing existing replaced)

- `Lagrangian/RunT3Batch2LaurentAudit.wl`
- `Numerical/parallel_t3/wilson/general_batch2_laurent.py`
- `scripts/RunT3ThesisGeneralBatch2Laurent.py`
- `tests/test_t3_thesis_general_batch2_laurent.py`
- `docs/T3_THESIS_BATCH2_NEXT_LAURENT_HANDOVER.md`

## Install and run in Windows PowerShell

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
$zip = 'C:\path\to\T3_FullGeneral_Batch2_Laurent_2026_11_10.zip'
Expand-Archive -LiteralPath $zip -DestinationPath (Get-Location) -Force

python -X utf8 -m unittest discover -s tests -p test_t3_thesis_general_batch2_laurent.py -v
python -X utf8 -m scripts.RunT3ThesisGeneralBatch2Laurent --all --preflight
python -X utf8 -m scripts.RunT3ThesisGeneralBatch2Laurent `
  --all `
  --run-id thesis_full_general_batch2_laurent_v1 `
  --accept-formal-audit-only
```

The preflight requires **both** saved raw matching folders and the complete preexisting Batch 1 result folders, but does not invoke Wolfram. `--accept-formal-audit-only` makes a successful *formal audit* return code 0; without that option, a successful audit returns code 3 because the complete finite numerical Wilson physics is not implemented. If the input or Mathematica audit blocks, the exit code is 2; check the per-model `manifest.json`, `laurent_audit.json`, and `wolfram_laurent.log`.

Default symbolic call timeout: 1200 seconds, with an overall subprocess timeout of 1800 seconds. If Mathematica aborts due to a long expression, use **a new run-id** and increase, for example, `--symbolic-seconds 2400 --timeout-sec 3000`. The audit performs no `Match[]` and therefore should be cheaper than the original fermion matching, but this has not been verified against your real exports.

## Outputs

Per model:

`output/parallel_t3/full_general_batch2/thesis_full_general_batch2_laurent_v1/{MODEL}/`

- `native_flavour_diagnostics.json` — MF singular values, conditioning, diagnostic Yukawa products in both possible inverse orientations, all live D extra operators; *not matched coefficients*.
- `raw_hbar1_coefficient.wxf` — full native operator source with its explicit loop-counting factor removed.
- `evaluated_full_hbar1_coefficient.wxf` — same full expression after Matchete loop-function evaluation.
- `explicit_pole_residue_unsubtracted.wxf` — if one simple explicit dimensional pole is identified.
- `epsilon_zero_regular_unsubtracted.wxf` — if a well-defined formal regular term is identified.
- `laurent_audit.json` — operator-preservation and regulator checks, preview, artifact SHA256, explicit non-certification of MSbar.
- `wolfram_laurent.log`, `manifest.json` — provenance and stage status.

The scalar-threshold Batch 3 remains a **later** step. The required EFT chain remains

`UV → integrate ALL F → EFT1 including BOTH scalar multiplets and Wilsons → integrate ALL scalars → SMEFT`.

## What to inspect after the run

Look first at `A_-4/laurent_audit.json` and `D_-2/laurent_audit.json`. If `status` is `LoopFunctionsEvaluated_ExplicitLaurentSeparated_AuditOnly`, the full native first-stage expression has a formal pole and regular part ready for **subtraction and numerical operator projection**. If `LoopFunctionsEvaluated_NoExplicitRegulator_AuditOnly`, inspect whether Matchete evaluated its loop functions without leaving an explicit epsilon symbol; this does **not** establish finiteness. If the audit is `Blocked`, the recorded reason identifies the precise unsupported symbolic structure. This is the next evidence needed to implement the missing finite Wilson and RGE stage without guessing any of the native gamma/CG conventions.
