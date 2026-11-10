# T3 A(-4) v8: invalid CG diagnosis and two-heavy-scalar trace

## Verified input (local v7 log)

* `CheckLagrangian` passes for the exploratory UV scalar base + W + `PlusHc` theory.
* `Match[..., EFTOrder -> 5, LoopOrder -> 0]` returns an EFT with no heavy scalar.
* The one-loop `Match` aborts after `DefineCG::notInvariant` (3.05 s).
* Removing the Wilson operator restores successful one-loop matching.
* The selected exploratory EFT is **not proven equivalent** to the first-threshold F EFT1; it is never a physical coefficient.

The Matchete implementation of `DefineCG` validates gauge-tensor invariance and aborts on failure. The existing source does not provide its arguments in the message. This package NEVER changes or disables that validation.

## Workstreams

1. With `T3_A4_TRACE_DEFINECG=1`, wrap the original one-loop `Match` in `TraceScan` matching only `DefineCG[...]`. The trace writes to `A_one_loop_DefineCG_calls_DIAGNOSTIC.txt` and flushes after each call so the last candidate survives `Abort[]`. This uses the official Wolfram `TraceScan`, not a mutation of Matchete definitions. Inspection of the final recorded call is mandatory before any CG correction is proposed.
2. With `T3_A4_PROBE_TWO_HEAVY_TRACE=1`, invoke Matchete's supported `CovariantLoop` on `{NewScalar1, NewScalar2}`. This targets the S1-S2 scalar bubble that can connect the inherited LLSS Wilson vertex with the Higgs-scalar mixing vertex. It is **only a subset of one-loop matching**. If it succeeds, raw WXF, each standard postprocessing stage and, if present, the extracted Weinberg coefficient are saved under `A_two_heavy_scalar_*_DIAGNOSTIC` names.
3. The original one-loop `Match` failure still blocks physical completion; no raw scalar trace is substituted into physical `C5`.

## Install / execute from repo root (PowerShell)

```powershell
Expand-Archive -LiteralPath "$HOME\Downloads\T3_A4_TraceCG_TwoScalarLoop_v8.zip" -DestinationPath . -Force
python -X utf8 -m unittest discover -s tests -p test_t3_a4_definecg_two_heavy.py -v
$env:T3_A4_EXPLORATORY_SOURCE_MATCH = '1'
$env:T3_A4_TRACE_DEFINECG = '1'
$env:T3_A4_PROBE_TWO_HEAVY_TRACE = '1'
python -X utf8 -m scripts.RunT3ThesisFinalOneBatchSafe --model A_-4 --run-id thesis_a4_definecg_scalartrace_v8 --debug
Remove-Item Env:T3_A4_EXPLORATORY_SOURCE_MATCH,Env:T3_A4_TRACE_DEFINECG,Env:T3_A4_PROBE_TWO_HEAVY_TRACE -ErrorAction SilentlyContinue
```

Inspect from `output/parallel_t3/full_general_final_one_batch/thesis_a4_definecg_scalartrace_v8/A_-4/`:

* `A_joint_scalar_match_diagnostics.json`
* `A_one_loop_DefineCG_calls_DIAGNOSTIC.txt` (especially the final recorded call)
* `A_two_heavy_scalar_Weinberg_DIAGNOSTIC.txt` if extraction succeeds
* `native_two_threshold_status.json`

## Scope and limitations

This uses the previously validated three-fermion-flavour W12 boundary; it does not repeat the F threshold, replace the two-real-scalar D(-2) treatment, invent a finite C5, or modify GitHub. In particular F hard one-loop matching, complete EFT1 Wilson RG and renormalisation are outstanding. This file **does not** provide an authoritative neutrino mass or NuFIT result.
