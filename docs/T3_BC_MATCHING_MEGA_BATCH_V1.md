# T3 B(-1)/C(-1) grouped scalar matching mega batch v1

## Scope

**One intermediate EFT only:** UV (fermions+scalars) -> EFT1 (two scalar doublets) -> SM+O5.

Reconstructs the complete existing *conditional* same-five-mass EFT1 Wilson boundary and integrates the existing seven-term master. Rotates the final 12-real-scalar Wilson and quartic tensors to the running X,Y mass basis at 1 TeV; evaluates **all internal scalar pair contributions with at least one X or Y**, using the project's already implemented zero-momentum MSbar scalar bubble `b0bar_zero`; and projects the resulting **HH** correction on C5. Then evolves baseline and bubble-shifted coefficients down to 100 GeV using the same SM boundary. No extra EFT stage.

The source quartic convention is `V = lambda_abcd phi_a phi_b phi_c phi_d / 4!`. Only diagrams with one LLphi phi insertion and a single quartic coupling are included. The bubble expression is

`delta C_ijab = -(1/(2*16*pi^2)) Sum_(c,d; at least one heavy) lambda_abcd C_ijcd B0bar(0;m_c^2,m_d^2;muS^2)`.

Differentiating the correction in ln(mu) yields `-heavy beta(C)/(16*pi^2)` holding the couplings/masses frozen. The finite part is a conditional diagnostic, not independently matched Matchete result.

## Numerical checks

- Saved five-mass trajectory hash inventory and continuity across F removal.
- Same five-mass EFT1 scalar fields at every Wilson integration node.
- Heavy Majorana basis residual, fermion hard-grouped embedding parity.
- 12-real-field X/Y rotation, scalar-potential and operator symmetries.
- Baseline/adjusted C5 and neutrino masses at 100 GeV.
- Fitted Weinberg projection residual before/after bubble.
- Coefficient log-derivative / heavy-quartic-beta cancellation.
- Matched-scale log sensitivity from muS/2 to 2muS, **frozen background only**.
- Extra native quartic beta sensitivity and second Wilson ODE tolerance.

## Missing physics / interpretation

The uncomputed matching pieces are **not** treated as numerically zero physical contributions: independent UV/EFT1 fermion counterterms, non-bubble one-loop scalar terms, other matching operators, canonical RG dictionary checks, and the full three-family Wilson anomalous dimension. The scalar-threshold bubble includes both masses in a grouped threshold at 1 TeV, not a Y->X intermediate EFT. Its full-mass finite contribution contains higher-order threshold-scheme dependence when scalar masses are nearby. A scan of frozen-loop logarithms is not an honest physical matching-scale independence test.

**Physical C5 / physical neutrino mass predictions remain NOT CERTIFIED.**

## Files

- `Numerical/experimental/T3BCGroupedFiniteScalarBubble.py`: grouped heavy-scalar loop kernel, log check, Weinberg projection.
- `scripts/RunT3BCMatchingMegaBatch.py`: end-to-end same-background source, Wilson, loop, SM continuation, diagnostics.
- `scripts/RunT3BCMatchingMegaBatch.ps1`: test+execution wrapper.
- `tests/test_t3_bc_matching_megabatch_v1.py`: independent kernel tests.
- `docs/T3_BC_MATCHING_MEGA_BATCH_V1.md`: this file.

## Running

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$env:USERPROFILE\Downloads\t3_bc_matching_megabatch_v1.zip" -DestinationPath . -Force
& .\scripts\RunT3BCMatchingMegaBatch.ps1
```

For faster initial execution (one Wilson solve instead of two per model):

```powershell
& .\scripts\RunT3BCMatchingMegaBatch.ps1 -SkipConvergence
```

Outputs: `output/studies/t3_bc_matching_megabatch_v1/Bm1_grouped_finite_bubble.json`, `Cm1_grouped_finite_bubble.json`, and `summary.json`. Existing production files are read-only.
