# Section 5.3: Zee-style compensated sensitivity, general T3 candidate

## Install paths

Copy the ZIP contents into the **local** repository root. This adds only:

- `scripts/RunT3GeneralZeeSensitivity.py` (new numerical driver)
- `tests/test_t3_general_zee_sensitivity.py` (pure unit tests)
- `docs/T3_GENERAL_ZEE_SENSITIVITY.md` (this guide)

No existing files are modified and no GitHub push is performed.

## Correct source and treatment

The older `scripts.RunT3FrozenFullScalar64 --plan` asks for missing historical
`output/studies/fermion_start_lambda_comparison/{5e-9,0p5}/fitted/*.json`.
**This new driver does not read those files.** It instead reads already fitted
*general-Lagrangian* candidates from:

`output/studies/general_scalar_extended/thesis_general_all16/thesis16_general_c5_full_v1/<MODEL>/`

- `fitted_config.json` and `result.json` for each model;
- `general_uv_overrides.json` for nonmixed models (including the saved extra quartics,
  allowed complex self-bilinears and extra Yukawas).

For the 14 nonmixed cases, the driver reuses:

- `scripts.RunT3GeneralAll16FromReference.nonmixed_inputs`,
- `scripts.RunT3GeneralAll16FromReference.NonMixedEvaluator`,
- `output/studies/general_scalar_extended/complete_scalar_symbolic/full16_general_v1/<MODEL>/{uv_rgbeta,eft1_rgbeta}.json`,
- `output/general_matchete_matching/<MODEL>/c5_coefficient.txt` and its matching status.

For B(-1)/C(-1), it reuses the separate **conditional** general-quartic + mixed-bilinear UV/EFT1 evaluator from
`RunT3GeneralScalarBCExploratoryFit`, with the saved one-scale mass-basis coefficient and
`output/general_t3_quadratic_eft1/*_live_rgbeta.json`.

The saved *general-fit* input for B/C is substituted into the existing reference-source
inventory; historical fitting directories are not read.

## UV deformations

For all three directions, `gamma` is positive. `lambda` denotes `lambdaT3`.

- T1: `lambda -> gamma*lambda`, `y1,y2 -> gamma**(-1/2)*(y1,y2)`.
- T2: `lambda -> gamma*lambda`, `y2 -> gamma**(-1)*y2`, `y1` unchanged.
- T3: `lambda` unchanged, `y1 -> gamma*y1`, `y2 -> gamma**(-1)*y2`.

This is the established T3 T1/T2/T3 convention, analogous to the *method* of
Zee section 4.4; T3 is **not** the exact Zee S3 transformation.

The nominal UV tensor `lambdaT3*(y1@y2.T + y2@y1.T)` is verified unchanged
at every gamma point. This is *not* a proof that the **complete general-C5**
source is unchanged: extra Yukawa channels and extra scalar interactions may
have additional dependence. The extra `y1Conj`/`y2Conj` Yukawas and independent
additional quartics/bilinears are held fixed at the UV boundary.

In the B/C adapter, the existing general beta evaluator sets extra quartics
as `quartic_fraction*abs(lambdaT3)`. To avoid accidentally moving independent
scalar couplings, the driver counter-scales `quartic_fraction` by
`abs(lambda_original)/abs(lambda_gamma)` for T1/T2. If this requires a fraction
above 1, it **fails** that point rather than using a changed scalar boundary.

Every point reruns the UV beta evolution (1e7 -> 1e5 GeV), fermion boundary,
scalar EFT1 beta evolution (1e5 -> 1e3 GeV), available one-scale scalar
matching, and SM+Weinberg evolution to 100 GeV.

## PowerShell from local repository root

```powershell
# Verify newly installed numerical driver's pure scaling contracts.
python -X utf8 -m pytest -q tests\test_t3_general_zee_sensitivity.py

# Check the GENERAL fitted/canonical symbolic inputs, with no physics run or outputs.
python -X utf8 -m scripts.RunT3GeneralZeeSensitivity `
  --models A_-2,B_-1,C_-1 --preflight

# One-model pilot, Zee paper ±1% window (33 full integrations total):
python -u -X utf8 -m scripts.RunT3GeneralZeeSensitivity `
  --models A_-2 --run-id zee_general_a2_v1 `
  --gamma-min 0.99 --gamma-max 1.01 --points 11

# Later, if the pilot and source checks succeed, scan A/B/C conditional models:
python -u -X utf8 -m scripts.RunT3GeneralZeeSensitivity `
  --models A_-2,B_-1,C_-1 --run-id zee_general_abc_v1 `
  --gamma-min 0.99 --gamma-max 1.01 --points 11

# Optional: all 16 model candidates, each on its own matching scheme:
python -u -X utf8 -m scripts.RunT3GeneralZeeSensitivity `
  --models all --run-id zee_general_all16_v1
```

Never reuse an output `--run-id`: a prior output folder causes an immediate
error. Failed points are recorded, not interpolated. A failed gamma=1 baseline
blocks the model from being labelled complete. The saved baseline predictions
are cross-checked against the new gamma=1 calculations.

## Outputs

`output/studies/general_scalar_extended/zee_style_general_sensitivity/<run-id>/`

- `summary.json`: per-model status and saved baseline chi-square;
- `points.csv`: model, T1/T2/T3, gamma, five low-energy observables, chi-square, errors;
- `<MODEL>/scan.json`: gamma points plus input SHA256 fingerprints and matching caveats;
- `<MODEL>/solar_sensitivity.png`: solar mixing versus solar mass-squared splitting, if matplotlib is installed.

The five computed observables are `delta_m21_sq_ev2`, `delta_m3l_sq_ev2`,
`sin2_theta12`, `sin2_theta13`, and `sin2_theta23`. Real-Yukawa inputs
only; CP-phase sensitivity is not computed.

## Interpretation constraints

1. **Not the historical reduced scalar potential.** All 14 nonmixed systems
   use canonical *general* UV and scalar-EFT RGBeta exports and the general
   Matchete source. B/C use a separate conditional general beta/C5 adapter.
2. **Not fully sequential general-potential matching.** For nonmixed systems,
   the saved scalar C5 is one-scale; for B/C it uses a conditional signed
   quartic map. Full fermion-first higher-dimensional Wilson matching and
   two-threshold physical certification remain unavailable in these evaluators.
3. **Not a fresh fit.** It rescans saved fitted inputs, some of which already
   have chi-square >5.86. The scan records those failed baseline screens.
4. **Not an isolated RG-only effect.** Each point uses a changed UV input.
   Changes can come from general matching as well as UV/EFT running. A dedicated
   running-vs-frozen-parameter comparison would be required to isolate RG effects.
5. **Method differences across representations.** Do not compare B/C and
   nonmixed absolute predictions as if they used a single certified matching
   prescription.

Development verification: Python syntax compiled and pure Python scaling
contract tests passed. **No full RG integration or Matchete evaluation was
possible in the artifact container**, which does not have the user's full local
repository or installed Wolfram environment. The local preflight and gamma=1
checks are therefore essential.
