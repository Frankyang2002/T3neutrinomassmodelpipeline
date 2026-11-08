# A(-4) RGBeta fresh-adjoint-index source correction

## Scope

This batch **replaces the full** `RGE/running/rgbeta/T3RGBetaModel.wl`
(1133 original lines). It changes exactly two A(-4)-active quartic-invariant
registrations: `T3RGBetaPortalAdjInvariant[3]` (for `lambdaH2Adj`) and
`lambdaS2Adj`. The contracted adjoint index is localised **within each pure
function call**, matching the independently evaluated `both` intervention.
No installed RGBeta package or GitHub source is modified.

The source is pinned to original SHA-256
`a58c79f18b2abdf9d476588552dfa52e2f131f322094e5b72d1198162aa20abc`,
and corrected SHA-256
`147b25f172ebcc3644e3ae1be8a7686fe66e1b7c1673ca17dd10e444ce0fff0d`.

## Physical and numerical significance

Previous native intervention (already run on Windows):

| Adjoint dummy-index treatment | Failed scalar-only probes |
| --- | ---: |
| Neither local | 34/55 |
| Triplet local | 27/55 |
| Portal local | 19/55 |
| Both local | 0/55 |

This batch goes beyond that isolated intervention by regenerating the **full
17-parameter** native A(-4) EFT1 beta export using the established
`RGE/running/rgbeta/RunGeneralT3RGBeta.wl`, validating all nine complete
quartic beta expressions at three nonzero gauge/SM-Yukawa benchmark points
against the prior successful native `both` intervention, validating unchanged
nonquartic beta channels, checking all 55 independent scalar-only Hessian
probes, and integrating the **renormalisable** EFT1 ODE to 1 TeV.

The newly solved trajectory is compared to BOTH previously saved seven-scale
results: archived original RGBeta and the previously generated five-monomial
scalar-only reference diagnostic. The corrected source need not agree with
the old scalar-only subtraction away from the gaugeless scalar subspace;
those differences are recorded, not suppressed.

## Required pre-existing data

These were produced in preceding successful user batches:

- `output/a_minus4_next_batch/fresh_index_causal_v2/native/both.json`
- `output/a_minus4_next_batch/conditional_quartic_ode/native/A_m4_eft1_tensor_basis.json`
- `output/a_minus4_next_batch/conditional_quartic_ode/reference_diagnostic/A_m4_eft1_tensor_basis.json`
- `output/a_minus4_eft1_batch2/A_m4_projected_raw.json`
- `output/a_minus4_eft1_batch2/A_m4_eft1_tensor_basis.json`

The already installed
`Numerical/running/Aminus4QuarticHessianIdentities.py` and
`scripts/RunAminus4FreshIndexCausalBatch.py` are imported strictly to reuse
the independent exact Hessian and verified 55-probe evaluator, respectively.
**Do not rerun the old four-mode experiment after replacing the source:**
that historical script deliberately pins the old source SHA-256.

## Outputs

`output/a_minus4_next_batch/fresh_index_production/`

- `A_m4_fixed_native_rgbeta.json`: newly exported native RGBeta results.
- `wolfram_export.log`: complete WolframScript command and stdout/stderr.
- `A_m4_fixed_17beta_compatible.json`: same new 17 beta functions with F-
  integrated EFT1 metadata required by the existing numerical solver.
- `running/A_m4_eft1_tensors.npz`: 45-Weyl/12-real-scalar sampled tensor trajectory.
- `running/A_m4_eft1_tensor_basis.json`: renormalisable coupling trajectory,
  current source provenance and explicit physical-certification blockers.
- `coupling_comparison.csv`: seven-scale quartic differences against the old
  five-monomial scalar-only diagnostic.
- `summary.json`: exact checks, SHA hashes, endpoint comparisons and blockers.

## Certification boundary

What this CAN support: native RGBeta source-level resolution of the A(-4)
scalar-only discrepancy and an internally consistent **renormalisable** EFT1
running trajectory.

What it does NOT support: absolute Matchete/RGBeta coupling normalization and
triplet CG phase mapping, Matchete `GammaCC` to master-Weyl signed factor,
physical `C_ijab(mu_F)`, full Wilson evolution, or simultaneous threshold
matching of both scalars to SMEFT. The previous physical-Wilson fail-closed
guard remains unchanged.
