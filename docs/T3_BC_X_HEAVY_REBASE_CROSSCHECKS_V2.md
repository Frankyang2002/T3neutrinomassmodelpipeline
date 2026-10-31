# B/C(-1) X-heavy continuation: exact rebase crosschecks (v2)

**Status:** Experimental native Matchete diagnostic. Physical sequential matching is **not certified**.

## Verified preceding local result

B(-1) `--check-x` returned `XFollowupRebasedEFT2CheckFailed`. Its Y-tree `EOMEFT` checked as all six true with X registered **light**, as did the Y-tree `MatchedEFT` comparator (not a permissible alternative matching input). The user's PowerShell output showed no contents for the requested `rebased_x_heavy_eft2_detailed_check` property, so it is **not possible to identify the failing heavy-X check** from that output. It must not be assumed that the rebase is a Hermiticity failure.

The V1 Wolfram file actually makes its decision by `And@@(TrueQ/@Values[check])` after `Expand[rawTree-freeLight+freeHeavy]//RelabelIndices`. The v2 runner exports *each* check independently and eliminates dependence on indirect association lookups for acceptance.

## Changes

Complete replacements, preserving the original numeric EFT2 output and all native operator sectors:

- `Lagrangian/RunT3BCXYNativeFollowup.experimental.wl`: after heavy-X re-registration, check the unchanged **Y-tree EOMEFT** directly in the heavy-X registry, as well as the existing free-Lagrangian-rebased candidate. Export independent six-flag results, input-form differences, and exact free equality. Reuse all of the successful independent Y-only native matching code. No operator is deleted or symmetrized.
- `scripts/RunT3BCXYNativeFollowup.py`: adds `--diagnose-x`, uses a **new v2 output directory** by default, source-hashes its inputs, and saves explicit pass/fail evidence. `--diagnose-x` runs both models even when an X check fails. A successful diagnostic can serve as the prerequisite for `--match-x` without rerunning `--check-x`.
- `scripts/CollectT3BCXYNativeEvidence.py`: complete replacement; when X matching succeeds, read the v2 X output directory and create an evidence archive with updated v2 filename. Requires Y and X stage outputs and retains physical certification flags false.
- `tests/test_t3_bc_xy_native_followup.py`: updated V1 regression tests for the stronger acceptance gate.
- `tests/test_t3_bc_xy_rebase_v2.py`: tests all six flags and negative controls, preserves Matchete and field provenance checks, checks source syntax/delimiters, and mocks a failing local Matchete result.

## What is a legitimate X-only native input?

Let `E2 = yTree["EOMEFT"]` after Y-only tree matching, and `FL`, `FH` be the X free Lagrangians constructed under the light and heavy registration respectively. Define

\[
E_{2,\mathrm{rebased}}=\mathrm{RelabelIndices}\left[\mathrm{Expand}(E_2-F_L+F_H)\right].
\]

V2 tests both `E2` and `E2_rebased` **after** re-registering X as heavy.

1. If the rebased expression passes all **six** full native Matchete checks, it may be used for the raw matching experiment.
2. Otherwise, the unchanged EOMEFT may be used **only if** all six full checks pass **and** `RelabelIndices[Expand[FH-FL]] === 0` is established in the same kernel. Then the replacement is algebraically trivial, and the source expression is unchanged. No numerical approximate equality is accepted.
3. Otherwise, it fails with `XFollowupRebasedEFT2CheckFailed`, and there is **no** X matching.

This is a fail-closed structural gate, **not** a proof of equivalence of the independent first threshold reconstruction to the original F-first EOMEFT, nor of inherited Y one-loop contributions.

## What to run (PowerShell)

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$env:USERPROFILE\Downloads\t3_bc_xy_rebase_crosscheck_v2.zip" -DestinationPath . -Force
python -m unittest discover -s tests -p "test_t3_bc_xy_native_followup.py" -v
python -m unittest discover -s tests -p "test_t3_bc_xy_rebase_v2.py" -v
python -u -X utf8 -m scripts.RunT3BCXYNativeFollowup --all --preflight-only
python -u -X utf8 -m scripts.RunT3BCXYNativeFollowup --all --diagnose-x
```

The v2 diagnostic is in:

`output/general_t3_mass_eigen_matchete/eft1_scalar_xy_native_followup_v2/B_m1_final_xy_native_diagnose_x_x_followup_rebase_crosschecks.json`

and the equivalent `C_m1...` file. A compact overview for each model is printed directly in PowerShell (direct-heavy check, rebased check, free-difference-zero).

If **both** models report `XFollowupEFT2XHeavyCheckPassedOnly`, then perform native X matching:

```powershell
python -u -X utf8 -m scripts.RunT3BCXYNativeFollowup --all --match-x
python -u -X utf8 -m scripts.CollectT3BCXYNativeEvidence --all
```

The resulting evidence ZIP will be:

`output/audits/t3_bc_xy_native_followup/xy_native_matching_evidence_v2.zip`.

If the X check remains false, do not run `--match-x`. Upload the two `*_rebase_crosschecks.json` files, the two `*_x_followup.json` files and their logs. The v2 comparison will isolate whether the defect is already present in the heavy registry or is introduced by replacing the free Lagrangian.

## Still unresolved

- The independently reconstructed Y matching Lagrangian and the saved F-first EOMEFT are not established equivalent. The former passes Matchete checks and gave raw amplitudes; this does not identify the latter with the physical F-first hard Wilson boundary.
- The X native matching currently acts on **Y-tree EOMEFT only**; inherited finite Y-loop matching is explicitly omitted.
- The original Matchete parameter point is not the `continuous_v1` running threshold point.
- Renormalized matching, counterterms, full three-generation SM Yukawa Wilson effects, and quartic convention certification remain open.

Do **not** replace the existing conditional neutrino masses with these symbolic probes or mark physical matching complete.
