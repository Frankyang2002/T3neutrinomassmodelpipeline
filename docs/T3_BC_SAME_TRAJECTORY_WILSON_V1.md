# B/C same-five-mass EFT1 Wilson candidate (v1)

Runs the existing **seven-term** `T3BCCoupledWilsonTensor` 12-real-scalar, six-lepton-pair solver using `eft1_trajectory` **directly** from the saved Workstream-B five-mass UV→EFT1 report. Heavy inputs are decoded from the *same* report's `uv_threshold.state`. No `continuous_v1` scalar RGE is read or re-integrated.

## Matching sequence

`UV(F,S1,S2) @ 1e7→1e5 GeV -> EFT1(S1,S2) @ 1e5→1e3 GeV -> SM+C5 @ 1e3 GeV -> SM @ 100 GeV`

There is one intermediate EFT. Both mixed scalars are grouped at 1 TeV, not separately integrated. The F-hard grouped boundary reuses the established historical grouped direct-source relocation convention, *not* the formal 68-term renormalised hard boundary. The EFT1 Wilson anomalous dimension uses the existing one-generation SM Yukawa approximation with six physical lepton flavour pairs. Stage29 RGBeta↔Matchete quartic map remains conditional.

## Source and certification

Requires Batch 1 audit summary, saved `B_m1`/`C_m1` five-mass reports, existing native source files, and project Wilson seed files. Rehashes five-mass inputs, verifies scalar-mass continuity and re-decodes each saved EFT1 state. Runs project `low_energy` from the **same** EFT1 endpoint.

Missing: independent UV/EFT1 finite counterterms, renormalized F-hard matching, finite scalar-threshold matching, complete three-family Wilson anomalous dimension, physical Matchete/RGBeta operator dictionary. **Not physical C5.** No production files overwritten. Numerical completion must be checked on local Windows repository.

## PowerShell

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
& .\scripts\RunT3BCSameTrajectoryWilson.ps1
```

Output: `output/studies/t3_bc_same_trajectory_wilson_v1/summary.json`, alongside each model's `*_same_trajectory_wilson.json`.
