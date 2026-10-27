# B_-1 / C_-1 scalar-threshold EOMEFT continuation fix

## Observed live failure

The scalar-threshold probe failed at `scalar-threshold-EFT1-lagrangian-check-started`, before second-threshold matching, returning `SavedEFT1TreeHeavyScalarCheckFailed` for both B_-1 and C_-1. This does **not** invalidate the successful fermion-first tree/one-loop runs or the 68+68+68 hard-pole audit.

## Verified source-level cause / corrective design

The failed scalar-threshold Python runner supplied `*_eft1_tree.inputform.txt`, namely `matchingResult["MatchedEFT"]` from the first-stage matching. The repository's authoritative `Lagrangian/RunMatching.wl`, `CanonicalTransitionEFT`, explicitly returns `EOMEFT` for the **next threshold**, because `ReplaceEffectiveCouplings` may undo canonical mass/kinetic field redefinitions. This is a verified architectural error in the previous scalar-probe handoff; whether it was the only cause of the live failure is **not yet proven** without rerunning local Matchete.

Separately, the second-stage Wolfram code wrapped `CheckLagrangian` inside `Quiet @ Check[..., $Failed]`. `Check` returns its failure alternative on *any* message, including nonfatal messages, as noted in `Lagrangian/RunMatching.wl`. First-stage Matchete validation already uses `CheckAbort[CheckLagrangian[...],$Aborted]`; the second stage now does likewise and records diagnostics. It does **not** ignore a failed canonicalization or bypass the heavy-mass-basis check.

## Changes

- `Lagrangian/RunGeneralT3FermionFirstEFT1Tree.wl` (full replacement): additionally exports `matchingResult["EOMEFT"]` as `*_eft1_transition_eom.inputform.txt` and certifies its type. Keeps existing expanded `MatchedEFT` and LLSS diagnostics unchanged.
- `scripts/RunGeneralT3FermionFirstEFT1Tree.py` (full replacement): validates the new canonical EOMEFT certificate and artifact; cleans stale transition artifacts before a run.
- `Lagrangian/RunGeneralT3EFT1ScalarThresholdProbe.wl` (full replacement): no longer marks message-emitting successful Matchete checks as `$Failed`. Logs detailed failures and **still fails closed** unless ordinary `CheckLagrangian` and detailed `HeavyMassBasis` checks pass.
- `scripts/RunGeneralT3EFT1ScalarThresholdProbe.py` (full replacement): reads the newly exported canonical EOMEFT from the validated first-stage **tree-only** run; rejects missing/stale transition files. Prior one-loop fermion-threshold outputs remain independently validated and are not rerun.
- Both corresponding `tests/python/test_*.py` are full replacements with canonical handoff regression tests.

The root `Lagrangian/RunMatching.wl` and all Matchete CGs/couplings are untouched. No GitHub changes. The one-loop hard pole/log/finite export and the numerical C5 benchmark code are untouched.

## PowerShell, repository root

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "C:\path\to\T3_General_EFT1_Canonical_EOM_Handoff_Fix.zip" -DestinationPath . -Force

python -m pytest -q tests\python\test_general_t3_fermion_first_eft1_tree.py tests\python\test_general_t3_eft1_scalar_threshold_probe.py

# Preserve existing successful tree-EFT1 outputs before re-running only its cheap tree match:
Copy-Item -Path "output\general_t3_mass_eigen_matchete\eft1_fermion_first_tree" -Destination "output\general_t3_mass_eigen_matchete\eft1_fermion_first_tree_pre_eom_backup" -Recurse -Force

python -m scripts.RunGeneralT3FermionFirstEFT1Tree --all
python -m scripts.RunGeneralT3EFT1ScalarThresholdProbe --model B_-1
python -m scripts.RunGeneralT3EFT1ScalarThresholdProbe --model C_-1
```

The backup command assumes `eft1_fermion_first_tree_pre_eom_backup` does not already exist; give it another name if necessary. **Do not rerun** `RunGeneralT3FermionFirstEFT1OneLoop` unless its prior outputs were removed, since the loop outputs are independent of the EOMEFT export addition.

Successful first-tree re-execution should create for each model:

`output/general_t3_mass_eigen_matchete/eft1_fermion_first_tree/{B,C}_m1_final_eft1_fermion_first_tree_eft1_transition_eom.inputform.txt`

The second-stage scalar-threshold outputs remain in:

`output/general_t3_mass_eigen_matchete/eft1_scalar_threshold_probe/`

## Validation and remaining uncertainty

Offline tests verify the canonical transition contract and reject stale `MatchedEFT` inputs, incomplete EOMEFT, or failed Matchete checks; the actual native scalar-threshold matching and physical sequential C5 remain **unverified** until the user's local Wolfram run. If the second stage fails again, use the new `B_m1_final_eft1_scalar_threshold_probe.log` / `C_m1_...log`; the Matchete messages and detailed check are now visible. The next scientific steps are EFT1 tensor/RGE mapping, hard pole/RGE consistency, scalar matching and physical C5 assembly.
