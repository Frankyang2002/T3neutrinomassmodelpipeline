# A(-4) and D(-2): native scalar continuation (local-only)

## Source baseline

Read-only GitHub `main` checked on 2026-11-10:
`c7313f048fe29f9f5c110616b1b40c34273eb220`.

This package is a **complete source replacement** of
`Lagrangian/RunT3FullGeneralTwoThresholdNative.wl` and
`tests/test_t3_thesis_final_one_batch.py`, plus a new
`Lagrangian/RunT3DMinus2FullGeneralPhysicalScalar.wl`.
It does not overwrite the successful Batch 1/Batch 2 outputs.

## Concrete changes

1. A(-4): removes the premature parent-kernel `DefineCoupling` call.
   The parent constructs the exact symbolic `W12[p,q]` EFT1 coefficient
   with the full `r=1,2,3` contraction and transports it by WXF to the
   existing **fresh scalar-stage kernel**, which already registers the
   Wilson couplings after `ResetAll`, SM loading, and T3 field setup.
   The F hard one-loop source is preserved separately.
2. D(-2): adds a physical scalar-stage Matchete continuation that reads
   the saved three-flavour `full_general_batch2` F boundary and verifies
   its real self-bilinear mass-basis metadata. It never repeats the
   unrotated F threshold. It reconstructs the CG registry, contracts
   `y1*y2/MF[r]` and `y1*y2Conj/MF[r]` into independent EFT1 Wilson
   flavour tensors, then promotes the three physical scalar masses from
   light to heavy EFT counting. It attempts **one** Match at LoopOrder=1
   removing `NewScalar1`, `T3Scalar2Minus`, `T3Scalar2Plus` together.
3. The tests check the source contract, forbidding the NF=1 fallback and
   preventing the partial scalar-soft coefficient from being marked
   as a finite physical C5.

## Limitations

The stage-1 hard one-loop operator boundary is *preserved*, not
renormalised/mapped into the independent Wilson tensor basis. The full
EFT1 Wilson RG evolution and F hard / scalar soft subtraction are still
missing. Neither scalar-soft output is a physical `C5`; the existing
Python wrapper correctly blocks the downstream numeric C5 and fitter.
The new D continuation is source-implemented but **not executed or
Matchete-verified in this sandbox**, which lacks Wolfram/Matchete.
Existing D legacy scalar candidate is NF=1 and is **not** substituted.

A zero exit code is intentionally **not** returned for a scalar-soft-only
result. The user must read
`output/parallel_t3/full_general_final_one_batch/<RUN_ID>/<MODEL>/native_two_threshold_status.json`
and the corresponding `native_two_threshold.log`.

## Install and run (PowerShell)

From the project root `C:\thesistemp\T3neutrinomassmodelpipeline`, extract
`T3_ThreeFlavour_PhysicalScalar_Continuation.zip` over the root. Keep the
existing outputs. The following is one single two-model execution, not a
re-run of the completed Batch1 or Batch2 stages:

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\T3_ThreeFlavour_PhysicalScalar_Continuation.zip" -DestinationPath . -Force
python -X utf8 -m unittest discover -s tests -p test_t3_thesis_final_one_batch.py -v
python -X utf8 -m scripts.RunT3ThesisFinalOneBatch `
  --all `
  --run-id thesis_physical_scalar_continuation_v1 `
  --debug
```

No `--fit` is appropriate until finite, normalised physical `C5` is
obtained for both models. If the run ID exists, use a fresh one.

## Evidence required for next implementation

Inspect A's `debug/fresh_kernel_stage_2/stdout.log`, D's
`native_two_threshold.log`, and both `native_two_threshold_status.json`
files. Their Matchete diagnostics will determine whether the new matching
stage now runs. Even if these stages pass, finite Wilson renormalisation
and EFT1 operator mixing remain necessary for a full-general result.
