# Batch 6 B(-1)/C(-1) — canonical EFT1 scalar CheckLagrangian recovery (v5)

## Observed failure

Both models in `mixed_hard_04` pass Python preflight and successfully run the original native **mixed-scalar mass-eigenstate** builder. The native scalar-only continuation then stops with:

```
Canonical EFT1 scalar CheckLagrangian returned False
last_stage=fermion-first-lagrangian-checks-returned
```

The current console excerpt does not include `CheckLagrangian[...,DetailedOutput->True]`, so the precise failed subcheck in that run is **not yet verified**. Prior source `Lagrangian/RunGeneralT3EFT1ScalarThresholdProbe.wl` explicitly documents the **light-to-heavy FreeLag rebase** and a potential `EOMEFT` serialization/Hermiticity representation issue. The repository's `scripts/ResolveT3BCEFT1HermiticityAndYGate.py` documents original-syntax Matchete Hermiticity checks that can be false despite native component equivalence. v5 implements these known, narrow recovery checks rather than disabling `CheckLagrangian`.

## What changes in v5

A full replacement of `scripts/RunT3MixedScalarHardBatch6.py`:

1. Keep the validated initial 11-vertex B/C baseline and seven native extra quartics, native epsilon CG rotations, and the physical scalar basis `Phi=(S1,epsilon S2*)=U (X,Y)`. Never enable original-basis mixed quadratics before rotation.
2. Resume from the **saved canonical EOMEFT** without calling any fermion `RunT3Matching` again. Remove F, preserve the real MF coupling, and capture `FreeLag[X,Y]` in the **light** registry.
3. Register both complex eigen-doublets **heavy**, regenerate the native **heavy** `FreeLag`, and rebase only the free action: `EFT1_heavy = EFT1_saved - FreeLag_light + FreeLag_heavy`. The saved file itself is untouched.
4. Run ordinary and **detailed** `CheckLagrangian` (Hermiticity, canonical normalization, heavy mass basis, gauge invariance and object definitions). Attempt equality-preserving rewrites only when needed.
5. If the checks still fail, reconstruct an alternative representation using **only the already validated native rotated UV action with every F-containing term removed**, plus **the saved expanded GammaCC/LLSS terms**. Require algebraic equality (allowing dummy-index relabelling) to the rebased saved EFT1, and independently require all five detailed Matchete checks. If either proof fails, abort. Never project away quartics, invent conjugates or replace the action by `(L+HC[L])/2`.
6. Only after successful complete checks attempt native one-loop scalar matching and export the provisional scalar-hard C5 candidate.

This is a guarded recovery attempt. It has **not** been run against a live Mathematica/Matchete kernel in the authoring environment; a further native blocker is possible.

## Files in this ZIP — paths relative to project root

- `scripts/RunT3MixedScalarHardBatch6.py` — **complete replacement**
- `tests/python/test_t3_mixed_batch6.py` — **complete replacement**
- `README_T3_SEQUENTIAL_BATCH6.md` — **complete replacement**

Leave `Numerical/running/T3MixedBatch6Numerics.py` and `scripts/RunT3MixedFixedOrderBatch6.py` from v3/v4 unchanged. No first-threshold inputs or any of the eight completed ordinary models are modified.

## PowerShell (repository root)

```powershell
Set-Location 'C:\thesistemp\T3neutrinomassmodelpipeline'
Expand-Archive `
  -LiteralPath "$HOME\Downloads\T3_Batch6_MixedScalar_CanonicalEFT1Repair_v5.zip" `
  -DestinationPath . -Force

python -X utf8 -m pytest -q tests/python/test_t3_mixed_batch6.py
python -X utf8 -m scripts.RunT3MixedScalarHardBatch6 --preflight-only
python -X utf8 -m scripts.RunT3MixedScalarHardBatch6 --run-id mixed_hard_05
```

If both native matches pass:

```powershell
python -X utf8 -m scripts.RunT3MixedFixedOrderBatch6 `
  --hard-run mixed_hard_05 --check-only
python -X utf8 -m scripts.RunT3MixedFixedOrderBatch6 `
  --hard-run mixed_hard_05 --run-id mixed_fixed_05 --conditional-quartic-map
```

Use a **new run ID**: v5 does not overwrite any historical `mixed_hard_01`–`mixed_hard_04` results.

## New diagnostics

For each model, in
`output/general_sequential_pilot/mixed_scalar_hard_batch6/mixed_hard_05/B_m1/`
(or `C_m1/`):

- `hard_scalar_input_audit.json`: preserved F-only source, MF and scalar mass format, heavy/light free-action rebase.
- `hard_scalar_eft1_validation.json`: first detailed Matchete check, equivalence-preserving attempts, accepted method (if any).
- `hard_scalar_reconstruction_audit.json`: saved-vs-reconstructed algebraic equality and independent Matchete checks (if attempted).
- `native_result.progress.json`, `native_result.json`, `stdout.log` and `stderr.log`: native stage and failure context.
- On successful matching only: `hard_c5_candidate.txt`, `hard_scalar_report.json` with status `MixedScalarHardCandidateExported`.

If either native model fails, inspect the JSON checks above before rerunning it. In particular, `full_EFT_validated=false` means the scalar matching was not attempted.

## EFT and physical status

This pilot follows `UV -> integrate out all fermions -> EFT1(SM+X+Y) -> integrate out X and Y -> SMEFT`. The scalar step is at a **common scale**, not separate Y then X thresholds. All three flavour numerical heavy-mass contributions rely on lifting the saved common-MF template to nondegenerate masses, not on native three-mass matching. The Matchete/RGBeta quartic map is conditional. The inherited first-threshold soft pole, full counterterms, pole cancellation, EFT1 Wilson running and certified physical `C5` remain **not established**.

Offline verification performed here: 20/20 Python tests pass, the Wolfram helper's brackets/comments/strings balance, and ZIP integrity is checked. These checks cannot certify the new Wolfram computation.
