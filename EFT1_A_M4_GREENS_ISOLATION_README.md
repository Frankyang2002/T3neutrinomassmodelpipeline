# A(-4) `GreensSimplify` isolated-portal diagnostic

**Project path:** `Lagrangian/TraceAminus4PortalGreens.wl` (new file). This README goes at the repository root. No existing source or numerical outputs are modified.

## Why this diagnostic is necessary

Your fresh trace established:

- `LUV`: both `lambdaH2Inv1` and `lambdaH2Inv2` occur.
- `RawEFT`: both occur.
- `GreenEFT`: only `lambdaH2Inv2` occurs.

`GreensSimplify` performs basis reductions. Counting source coupling names is insufficient to tell whether an independent operator was lost, or whether the source was canceled, contracted, redefined or absorbed into another invariant. Therefore **do not** add `lambdaH2Inv1` by hand or bypass the existing Stage B general-Lagrangian gate.

## Run (PowerShell, from repository root)

```powershell
Set-Location 'C:\thesistemp\T3neutrinomassmodelpipeline'
$env:T3_EXPERIMENTAL_GENERAL_FERMIONS = '1'
$env:T3_EXPERIMENTAL_GENERAL_QUARTICS = '1'
$env:T3_EXPERIMENTAL_GENERAL_QUADRATICS = '1'
$env:T3_EXPERIMENTAL_FLAVOR_MASSES = '1'

wolframscript -file Lagrangian\TraceAminus4PortalGreens.wl `
  output\eft1_general_lagrangian_trace\A_m4_greens_01

Get-Content 'output\eft1_general_lagrangian_trace\A_m4_greens_01\greens_portal_isolation.json' -Raw
```

Use a new output directory if `A_m4_greens_01` already exists.

## What it checks

1. Rebuilds **the actual general** A(-4) `LUV` with all four general-operator switches set.
2. Runs the same `RunT3Matching[... , 5, 0]` as the existing fermion-first matching pipeline.
3. Extracts entire expanded `RawEFT` summands containing each portal coupling, **without stripping scalar fields, CGs or indices**.
4. Calls Matchete's `GreensSimplify` independently on the `lambdaH2Inv1` terms, the `lambdaH2Inv2` terms, and the original `RawEFT` with only `lambdaH2Inv1` terms subtracted.
5. Saves each input/output as full `InputForm` and reports whether `lambdaH2Inv1` contributes a nonzero difference in the fully simplified EFT.
6. If the two Matchete portal CG tensors can be read, reports their raw linear rank. Rank 2 supports **independence of the CG arrays in that frame** but does not establish a physical Matchete/RGBeta coupling dictionary.

**Report:** `output/eft1_general_lagrangian_trace/A_m4_greens_01/greens_portal_isolation.json`.

Key files to inspect if the report is ambiguous:

- `RawEFT_lambdaH2Inv1_only.inputform.txt`
- `GreenEFT_lambdaH2Inv1_only.inputform.txt`
- `GreenEFT_without_lambdaH2Inv1.inputform.txt`
- `GreenEFT_difference_with_without_inv1.inputform.txt`

## Meaning of outcomes

- `Inv1IsolatedGreensZero_CheckSU2CGIdentity`: isolated reduction itself kills the `Inv1` terms. Check whether the source's CG actually yields a nonzero scalar operator after contraction. This is **not proof of a Matchete bug**.
- `FullGreensIndependentOfInv1_CheckRebasisOrCancellation`: output is syntactically unchanged by removing `Inv1`; check cancellation or renormalisable-coupling redefinitions before diagnosing a lost operator.
- `Inv1SurvivesIsolatedReduction_CheckFullRecombination`: the isolated source survives, but in the complete theory it is eliminated or rewritten. Compare `GreenEFT_difference_with_without_inv1.inputform.txt` and inspect its effective couplings and CGs.
- `OneOrMoreIsolatedGreensCallsFailed`: Matchete cannot simplify the isolated partial Lagrangian on its own. The full-source stage is still valid; do not infer a zero contribution.

**Limits:** This performs no Matchete source-code alteration, physical `GammaCC`/Weyl normalization, Matchete/RGBeta scalar basis-map certification, Wilson-ODE running, or scalar-threshold matching. Full general-Lagrangian certification remains blocked until the origin of the missing source coefficient is resolved.
