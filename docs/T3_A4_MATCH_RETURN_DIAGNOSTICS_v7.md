# A(-4) Matchete raw `Match` classification — v7

## Source and limitations

This file replaces `Lagrangian/RunT3AMinus4FullGeneralPhysicalScalar.wl`, based on the actual v6 ZIP provided in the conversation and the GitHub reference `bad407a7` (November 10, 2026). The v6 audit reported a Hermitian UV-derived scalar base plus the original holomorphic three-flavour Wilson term and a Matchete-generated Hermitian conjugate. That *exploratory* Lagrangian passed `CheckLagrangian` but gave `Head[Match[...]] === Symbol`.

The original saved EFT1 did **not** pass the equivalence test. The exploratory theory therefore remains diagnostic only; it cannot be used as an authoritative UV-to-SMEFT Wilson coefficient.

## Changes

The rest of the full v6 file is preserved. At the joint scalar matching stage:

1. Match the **same** selected input once at `EFTOrder -> 5, LoopOrder -> 0` as a structural control.
2. Match that same input at `EFTOrder -> 5, LoopOrder -> 1` if the tree stage returned an appropriate expression.
3. Write `A_joint_scalar_match_diagnostics.json` before any failure exit. The report records the exact return identity (`$Failed`, `$Aborted`, `$TimedOut`, another atom, or a nonatomic expression), raw expression head, elapsed time, hbar presence, remaining scalars, and any captured Mathematica message labels.
4. Optionally, after a failed full one-loop Match, run a scalar-base-only one-loop diagnostic *in the same kernel*, controlled by `T3_A4_PROBE_WILSON_FREE_MATCH=1`. **Never** reuse that result as the T3 Wilson or Weinberg coefficient. This optional test is not independent because a prior failed Matchete invocation may have changed internal state.

### Intended project path

- Replace: `Lagrangian/RunT3AMinus4FullGeneralPhysicalScalar.wl` (full file)
- Add: `tests/test_t3_a4_matchete_match_triage.py`
- Add: `docs/T3_A4_MATCH_RETURN_DIAGNOSTICS_v7.md`

## Windows PowerShell commands

From the existing project root (do not overwrite result directories):

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\T3_A4_Matchete_Return_Triage_v7.zip" -DestinationPath . -Force
python -X utf8 -m unittest discover -s tests -p test_t3_a4_matchete_match_triage.py -v

$env:T3_A4_EXPLORATORY_SOURCE_MATCH = "1"
$env:T3_A4_PROBE_WILSON_FREE_MATCH = "1"
python -X utf8 -m scripts.RunT3ThesisFinalOneBatchSafe --model A_-4 --run-id thesis_a4_match_return_triage_v7 --debug
Remove-Item Env:\T3_A4_EXPLORATORY_SOURCE_MATCH -ErrorAction SilentlyContinue
Remove-Item Env:\T3_A4_PROBE_WILSON_FREE_MATCH -ErrorAction SilentlyContinue

$dir = 'output\parallel_t3\full_general_final_one_batch\thesis_a4_match_return_triage_v7\A_-4'
Get-Content "$dir\A_joint_scalar_match_diagnostics.json"
Get-Content "$dir\native_two_threshold_status.json"
Get-Content "$dir\native_two_threshold.log" -Tail 120
```

The wrapper uses a binary-safe capture path, preventing UTF-8 subprocess reader crashes. Use only a **new run ID**.

## Interpretation

- **`tree_control.return_status` not `ExpressionReturned`**: scalar matching input or tree-level Matchete path fails before one-loop machinery; inspect messages in the raw log.
- **Tree passes, one-loop `$Failed`**: explicit Matchete one-loop failure; inspect message names/errors. The optional scalar-base control can help identify whether the dimension-five Wilson insertion triggers it, but its result is diagnostic-only.
- **Tree passes, one-loop `$Aborted`**: likely an explicit Matchete abort or a Wolfram abort; do not infer a numerical problem without logs.
- **One-loop `$TimedOut`**: timing or computation-resource blocker, not a proof of invalid operators.
- **One-loop `ExpressionReturned`, remaining heavy fields non-empty**: matching output has not integrated the intended scalar fields.
- **One-loop `ExpressionReturned`, no explicit `hbar`**: source produces no confirmed one-loop sector; review operator inventory and logs.

No native Matchete result has been verified in this environment. Even if the exploratory scalar match succeeds, the original F boundary equivalence, fermion-hard inheritance, full EFT1 Wilson anomalous dimensions and finite subtraction are still required before physical C5 and fitted neutrino data can be claimed.
