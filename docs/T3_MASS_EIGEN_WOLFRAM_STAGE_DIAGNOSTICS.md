# Mixed-doublet Matchete live-stage diagnostics

## Reason for replacement

The live `B_-1` eigenbasis run reached Wolfram but returned only
`UncaughtMatcheteError`, exit 5. In the previous Wolfram file that status
is written when `UsingFrontEnd[Catch[T3MEPerform[], ...]]` returns `$Failed`,
`$Aborted`, or a non-association. The result did not preserve the failing
stage. Thus **there is no evidence yet whether `HeavyMassBasis` is valid.**
The provenance mismatch and the numerical mass rotation are already solved.

## Changes (experimental only)

- The existing `Lagrangian/RunGeneralT3MassEigenMatchete.wl` is replaced in
  full. It writes a separate stage-progress JSON before each high-risk step
  and emits `T3ME_STAGE: ...` in stdout. `UncaughtMatcheteError` now includes
  the last stage, a short stage history, the kernel result and message list.
- The existing `scripts/RunGeneralT3MassEigenMatchete.py` is replaced in full.
  It deletes stale progress files and reports the checkpoint and most recent
  output when live Wolfram exits unsuccessfully. SHA-256 provenance checks
  are unchanged and Windows/POSIX key normalization is preserved.
- `tests/python/test_general_t3_mass_eigen_matchete.py` is replaced in full.
  The tests enforce the checkpoint contract and mock an uncaught Matchete
  failure with/without progress data.
- This note is new.

The Lagrangian, scalar rotation, source inventory, CG conventions, and
`CheckLagrangian` pass/fail criteria are **unchanged**. This package does not
fix a yet-unidentified native Matchete issue; it exposes where it fails.

## Install (PowerShell, project root)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "C:\path\to\T3_MassEigen_Wolfram_Stage_Diagnostics.zip" -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_mass_eigen_matchete.py tests\python\test_general_t3_mixed_scalar_mass_basis.py
python -m scripts.RunGeneralT3MassEigenMatchete --model B_-1
```

`--all` is **not necessary** while `B_-1` fails. After a failure, inspect:

```powershell
Get-Content .\output\general_t3_mass_eigen_matchete\B_m1_final_matchete.progress.json
Get-Content .\output\general_t3_mass_eigen_matchete\B_m1_final_matchete.json
Get-Content .\output\general_t3_mass_eigen_matchete\B_m1_final_matchete.log -Tail 120
```

The `last_stage` and nearby Wolfram messages identify the specific operation
that needs correction. `B_-1` and `C_-1` still require **successful** ordinary
`CheckLagrangian` and `HeavyMassBasis` before the mass-eigenstate Lagrangian
can be considered validated.

## Verification

67/67 offline Python tests passed in a workspace containing the previous
mass-rotation source and tests. No Wolfram/Matchete kernel is available in
this packaging environment: the native `CheckLagrangian` remains unverified.
