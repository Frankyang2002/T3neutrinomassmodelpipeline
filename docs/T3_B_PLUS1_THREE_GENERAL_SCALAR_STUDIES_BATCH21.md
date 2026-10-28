# Batch 21 — B(+1) three-scenario conditional general-scalar numerical fits

## Scope

Runs exactly the original three study prescriptions for **B(+1)**, using the
original `lambdaT3` values and each scenario's historical fitted real Yukawa
matrices as seeds:

| `--case` | UV interval | Original lambdaT3 |
|---|---|---:|
| `small_uv_frozen_5e-9` | **frozen**, boundary imposed at 100 TeV | 5e-9 |
| `large_uv_frozen_0p5` | **frozen**, boundary imposed at 100 TeV | 0.5 |
| `small_uv_running_5e-9` | **coupled UV RGE**, 10^7 GeV to 100 TeV | 5e-9 |

After the fermion boundary the **native Batch20 RGBeta 18+1** system is
integrated from 100 TeV to 1 TeV. All three fit the 18 *real* Yukawa entries;
`lambdaT3` remains fixed. The new extra complex quartic is assigned a **real,
illustrative RGBeta boundary of -0.01 by default**, independently of the
small/large topology quartic. That is a numerical choice, not a derived physical
coupling map. Alternative: `--extra-real <value>` or `--extra-fraction 0.1`
(the latter uses `0.1 * abs(lambdaT3)` as the extra boundary instead).

The C5 evaluator uses **the model-specific existing**
`output/full/hypercharge/T3_B_alpha_p1/data/final_weinberg_coefficient.json`,
containing the hierarchical hard and direct logarithmic term already used by
Stage36. Batch19 showed no direct extra-quartic dependence in the one-loop
Matchete coefficient; **the new quartic affects C5 through RG running** in
this candidate approximation. A fast one-generation kernel is derived from
that existing symbolic source and checked against Stage36's full symbolic
`evaluate_final_weinberg_json` at the first mass-safe trial (relative error
must be < 3e-8). It is not a newly matched general-scalar Wilson result.

All three candidates evolve SM+Weinberg from 1 TeV to 100 GeV and use the
same five-observable NuFIT target and charged-lepton mass-basis fit as the
existing project. Refit cost is reduced in frozen UV by evaluating the scalar
EFT1 trajectory once; UV-running points independently rerun UV and EFT1.

## Install

Extract the archive into the project root. It adds:

- `scripts/FitT3BPlus1GeneralThreeStudies.py`
- `tests/python/test_t3_b_plus1_general_three_studies.py`
- `docs/T3_B_PLUS1_THREE_GENERAL_SCALAR_STUDIES_BATCH21.md`

It depends on **Batch20** being installed and on a successful, unchanged live
Batch20 output at
`output/studies/general_scalar_extended/b_p1_eft1_running/b_p1_eft1_general_v1/`.
The native beta SHA256 and original zero-extra control status are rechecked.
The Stage31 B(+1) coupled UV beta export and original three fitted source
configs must also exist. No Mathematica call occurs during the fit itself.

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
python -m pytest -q tests\python\test_t3_b_plus1_general_three_studies.py

# First run a cheap 3-scenario pilot (historical + 4 Sobol points, no local fit).
python -m scripts.FitT3BPlus1GeneralThreeStudies `
  --all --run-id b_p1_three_pilot_v1 --sobol-points 4 --max-nfev 0

# Full pilot optimisation. Reuse the same Batch20 betas, but a NEW run ID.
python -m scripts.FitT3BPlus1GeneralThreeStudies `
  --all --run-id b_p1_three_fit_v1 --sobol-points 8 `
  --max-nfev 25 --seed-count 1
```

Read outputs for every case:

```powershell
Get-Content -Encoding UTF8 output\studies\general_scalar_extended\small_uv_frozen_5e-9\b_p1_three_fit_v1\summary.md
Get-Content -Encoding UTF8 output\studies\general_scalar_extended\large_uv_frozen_0p5\b_p1_three_fit_v1\summary.md
Get-Content -Encoding UTF8 output\studies\general_scalar_extended\small_uv_running_5e-9\b_p1_three_fit_v1\summary.md
```

Each case writes a new immutable `fit.json`, `summary.md`, and (if a mass-safe
candidate exists) `fitted_config.json`. No historical input/output is changed.
The fitter prints the status and chi2 for each trial as it evaluates it.
When a case has no accepted fit, the report distinguishes `PoorFit` from
`NoMassSafeNeutrinoFit` and records up to 15 numerical or mass-screen failures.
Do not label these cases successful without checking their live summaries.

## Numerical screens

- `rhoF = max(sqrt(mS1Sq(muF)),sqrt(mS2Sq(muF)))/min(MF_singular(muF)) <= 0.5`.
- Positive scalar *running quadratic* masses at 1 TeV within 0.5–5 TeV.
- Yukawa maximum singular value at the fermion threshold <= sqrt(4 pi).
- Scalar coordinate magnitudes at the 1 TeV endpoint and the extra-quartic
  boundary <= 4 pi. This **does not prove** vacuum stability or a global
  perturbativity bound over every scale.
- Recomputed five-observable chi2 <= 5.86 constitutes a *candidate accepted*
  fit. The fit has 18 Yukawa entries for five observables, so it is not a
  prediction of the angles.

## Limitations

- B(+1) only; this does not add all remaining 42 results.
- Matchete–RGBeta general quartic coefficient normalization remains conditional.
- The reference Stage36 B(+1) hard-plus-direct C5 is reused; the extra quartic
  is not directly rematched at either threshold.
- No matched `C12`/Wilson running from the fermion threshold; grouped scalar
  matching is an approximation.
- The displayed masses are unbroken-phase running quadratic masses, not pole
  masses; vacuum stability and unitarity are not tested.
