# A(−4): Hermitian Higgs–triplet portal trial

**Input evidence:** Native Matchete CG data show `T3HiggsPortal2CGInv1` is anti-Hermitian under conjugate-transpose of the Higgs and triplet slots; `Inv2` is Hermitian. The isolated real-coupling `Inv1 + h.c.` expression was eliminated by `GreensSimplify`.

**Trial only:** `Lagrangian/ProbeAminus4HermitianPortal.wl` loads the existing current general Lagrangian builder and temporarily overrides only the construction of `HiggsPortal2Inv1` in that Wolfram process. It uses

```
(1/2) PlusHc[I lambdaH2Inv1[] (Hdag H S2dag S2)_Inv1]
```

instead of the original `(1/2) PlusHc[lambdaH2Inv1[] ...]`, after verifying the exact native CG is anti-Hermitian and checking the `Inv2` tensor is Hermitian. The trial then runs a fresh tree-level F-only Matchete match, tests the isolated corrected portal through `GreensSimplify`, and records the final EFT1 inventory. No existing `.wl` files are overwritten, and existing Stage A/B and audit outputs are untouched.

## PowerShell (repository root)

```powershell
Expand-Archive -LiteralPath "$HOME\Downloads\T3_EFT1_Aminus4_Hermitian_Portal_Trial.zip" -DestinationPath .
$env:T3_EXPERIMENTAL_GENERAL_FERMIONS = '1'
$env:T3_EXPERIMENTAL_GENERAL_QUARTICS = '1'
$env:T3_EXPERIMENTAL_GENERAL_QUADRATICS = '1'
$env:T3_EXPERIMENTAL_FLAVOR_MASSES = '1'
wolframscript -file Lagrangian\ProbeAminus4HermitianPortal.wl `
  output\eft1_general_lagrangian_trace\A_m4_hermitian_trial_01
Get-Content 'output\eft1_general_lagrangian_trace\A_m4_hermitian_trial_01\hermitian_portal_trial.json' -Raw
```

A successful diagnostic status is `TrialBothPortalsSurvive_DoNotTreatAsPhysicalBasisCertificate`. This is a hypothesis to be tested in your Matchete installation, **not a confirmed result**. It does not certify the Matchete-to-RGBeta dictionary, sign, normalization, Wilson physical normalization, or hard matching. The existing general-Lagrangian audit remains blocked until the production source is explicitly corrected, regenerated, and retested.

## Project paths included

- `Lagrangian/ProbeAminus4HermitianPortal.wl` — standalone additive trial.
- `EFT1_A_M4_HERMITIAN_PORTAL_TRIAL.md` — this file.
