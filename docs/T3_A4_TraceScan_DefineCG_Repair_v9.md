# A(-4) Matchete one-loop CG failure: v9 TraceScan repair

## Verified from the user's v8 execution

- Full three-flavour A(-4) exploratory Hermitian source passed `CheckLagrangian`.
- Joint scalar tree `Match[..., LoopOrder -> 0]` succeeded (v7).
- Joint scalar one-loop `Match[..., LoopOrder -> 1]` in v7 aborted with `DefineCG::notInvariant`.
- The v8 `TraceScan` wrapper did **not** run Matchete: it returned unevaluated `TraceScan` with `TraceScan::optx`, because `TraceOriginal -> True` is not a `TraceScan` option.
- A separate two-heavy-scalar `CovariantLoop` in v8 independently aborted with `DefineCG::notInvariant`.
- Wilson-free scalar one-loop matched successfully in v7. It does **not** yield physical Weinberg matching.

## Change

`TraceScan[f, expr, form]` is documented and `TraceInternal -> True` is a documented option. `TraceOriginal` is **not** in `Options[TraceScan]`; it belongs to `Trace`. v9 removes this option from the **full replacement** `Lagrangian/RunT3AMinus4FullGeneralPhysicalScalar.wl`, so Matchete's `Match` is actually evaluated.

The same valid trace wrapper now surrounds the optional targeted `CovariantLoop` call. Each path has its own flushed file:

- `A_one_loop_DefineCG_calls_DIAGNOSTIC.txt`
- `A_two_heavy_DefineCG_calls_DIAGNOSTIC.txt`

The relevant `...diagnostics.json` entries provide separate return statuses and trace call counts. An unevaluated `TraceScan` is recognized as an instrumentation failure, not an EFT result.

The visitor does **not** redefine `DefineCG`, alter a CG tensor, bypass `TensorInvariantQ`, change the scalar model, or modify any fermion-threshold source. Matchete's invariant assertion remains mandatory.

## Run (PowerShell, project root)

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\T3_A4_TraceScan_Repair_v9.zip" -DestinationPath . -Force
python -X utf8 -m unittest discover -s tests -p test_t3_a4_tracescan_repair_v9.py -v
$env:T3_A4_EXPLORATORY_SOURCE_MATCH = '1'
$env:T3_A4_TRACE_DEFINECG = '1'
$env:T3_A4_PROBE_TWO_HEAVY_TRACE = '1'
python -X utf8 -m scripts.RunT3ThesisFinalOneBatchSafe --model A_-4 --run-id thesis_a4_tracescan_repair_v9 --debug
Remove-Item Env:T3_A4_EXPLORATORY_SOURCE_MATCH,Env:T3_A4_TRACE_DEFINECG,Env:T3_A4_PROBE_TWO_HEAVY_TRACE -ErrorAction SilentlyContinue
```

Inspect output at `output/parallel_t3/full_general_final_one_batch/thesis_a4_tracescan_repair_v9/A_-4/`:

```powershell
$dir = 'output\parallel_t3\full_general_final_one_batch\thesis_a4_tracescan_repair_v9\A_-4'
Get-Content "$dir\A_joint_scalar_match_diagnostics.json"
Get-Content "$dir\A_one_loop_DefineCG_calls_DIAGNOSTIC.txt" -Tail 12
Get-Content "$dir\A_two_heavy_DefineCG_calls_DIAGNOSTIC.txt" -Tail 12
Get-Content "$dir\native_two_threshold_status.json"
```

The latest attempted `DefineCG[symb,indexTypes,tensor]` from the flushed trace should distinguish a representation mismatch from a newly generated invalid contraction. The next *physics* fix must be based on that actual failing call. In particular, a dimension-five CG cannot be replaced by a manually guessed invariant.

## Scientific status

The source-rebuilt A(-4) EFT has **not** been proven equivalent to the original saved F-threshold EFT1. This remains exploratory even if its soft match succeeds. Full UV fermion-hard inheritance, EFT1 Wilson RG, finite renormalisation, and pole cancellation are not yet verified. D(-2)'s previously obtained 68 scalar-soft Weinberg terms are unchanged; these also do not establish a complete physical C5.
