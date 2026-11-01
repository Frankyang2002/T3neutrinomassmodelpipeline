# Full-general Matchete C5 flavour channels: local 16-model exploratory fit

**What changes:** Replaces the old *restricted/original-topology* C5 kernel in the
14 non-mixed T3 model fits with C5 matched from the **full symmetry-allowed**
Matchete Lagrangian. The existing RGBeta full-general UV/EFT1 evolution,
Sobol sampling, least-squares refinement, low-energy Weinberg evolution,
NuFIT observable extraction and B(-1)/C(-1) mixed-scalar branch are reused.

The Matchete coefficient is decomposed into every distinct quadratic Yukawa
channel using `Lagrangian/EvaluateGeneralC5Channels.wl`. Each resulting
coefficient is combined with **the actual evolved UV heavy-mass-basis Yukawa
matrices**, including `y1Conj`/`y2Conj` when these occur in the full-general
coefficient, before symmetrising the Majorana flavour matrix. No missing
channel is replaced with a restricted C5 kernel.

This is still provisional **one-scale** matching, not a complete sequential
fermion-first EFT1 Wilson calculation; the Matchete/RGBeta quartic dictionary
and same-heavy-leg rotations of additional Yukawas are assumptions recorded
in diagnostics. B(-1)/C(-1) retain their existing direct mixed-scalar matching.
The 18 real entries of y1/y2 are fitted; the extra Yukawa boundary is fixed
relative to y1 or y2 and evolved independently through the general RGEs.

## Full replacement files and project destinations

- `Numerical/matching/GeneralT3FlavourC5.py` **NEW**: channel-based complex
  flavour C5 evaluation, using the existing Wolfram kernel evaluator.
- `scripts/RunT3GeneralAll16FromReference.py` **REPLACE**: route non-mixed
  models through full-general Matchete instead of original topology C5.
- `tests/python/test_thesis_general_c5.py` **NEW**: test flavour pairings,
  relative complex phases, channel completeness and missing-Yukawa detection.
- `tests/python/test_thesis_general16.py` **REPLACE**: updated campaign tests.

Keep `scripts/RunT3BCGeneralFromReference.py` from the previously installed
16-model ZIP unchanged. No original physics modules or GitHub files are edited.

## Installation

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive "$HOME\Downloads\T3_Full_General_C5_Flavour_and_16_Fits.zip" -DestinationPath . -Force
python -X utf8 -m pytest -q tests/python/test_thesis_general_c5.py tests/python/test_thesis_general16.py
```

## 1. Build the full-general Matchete C5 for one model

```powershell
python -X utf8 -m scripts.RunT3GeneralMatcheteMatching --model A_-2
```

This uses `RunModel.wl` with `T3_EXPERIMENTAL_GENERAL_FERMIONS=1`,
`T3_EXPERIMENTAL_GENERAL_QUADRATICS=1`, and
`T3_EXPERIMENTAL_GENERAL_QUARTICS=1`. The coefficient should be at
`output/general_matchete_matching/A_m2/c5_coefficient.txt` with
`general_matching_status.json` reporting `Matched`. The general Matchete
symbolic matching can fail for individual models; inspect
`general_matching.log` in that model directory rather than using a restricted
coefficient. The matching generator already excludes B(-1)/C(-1).

## 2. First end-to-end flavour-C5 numerical smoke fit

```powershell
python -X utf8 -m scripts.RunT3GeneralAll16FromReference `
    --run-id general_c5_A_m2_test --models A_-2 `
    --generate-general-c5 --sobol-points 2 `
    --max-nfev 8 --no-local-refinement
```

This needs the existing canonical full-general UV/EFT1 symbolic outputs under
`output/studies/general_scalar_extended/complete_scalar_symbolic/full16_general_v1/`
(or supply `--canonical-run-id <your ID>`).
The calculation uses the normal ordering NuFIT target and writes the predicted
mixing angles and splittings to `observables.csv`.

## 3. Full 16-model fit

Use **a new run ID**, since `thesis16_v1` used the previous restricted-C5
matching scheme. B(-1) and C(-1) can reuse their saved direct mixed-scalar runs:

```powershell
python -X utf8 -m scripts.RunT3GeneralAll16FromReference `
    --run-id thesis16_general_c5_v1 `
    --reuse-b-run thesis_b_m1_128 `
    --reuse-c-run thesis16_v1 `
    --generate-general-c5 --export-bc-betas `
    --sobol-points 128 --max-nfev 60
```

`--reuse-c-run thesis16_v1` works only if the previous C(-1) direct fit exists
at `output/studies/general_scalar_extended/fresh_reference_bc/thesis16_v1/C_m1/`.
If that saved direct result is absent, remove `--reuse-c-run` to reoptimise C(-1).

An initial `--preflight` without `--generate-general-c5` lists absent source
coefficients but does not generate them. Fitting with
`--generate-general-c5` calls the existing general Matchete matching generator
when a coefficient is missing. Each model is independently marked as fitted
or blocked, with errors in `<model>/error.json`.

Resume without changing the physics or Sobol settings:

```powershell
python -X utf8 -m scripts.RunT3GeneralAll16FromReference `
    --run-id thesis16_general_c5_v1 --resume `
    --reuse-b-run thesis_b_m1_128 --reuse-c-run thesis16_v1 `
    --generate-general-c5 --export-bc-betas `
    --sobol-points 128 --max-nfev 60
```

## Results

`output/studies/general_scalar_extended/thesis_general_all16/<run_id>/`
contains `summary.csv`, `summary.json`, `settings.json` and one folder per
model containing `result.json`, `fitted_config.json`,
`general_uv_overrides.json`, and `observables.csv` on successful non-mixed
fits. `result.json` records `matching_channels`, heavy-mass-basis residual,
quartic mapping assumptions, fitted Yukawas and conditional C5 limitations.

The interactive RGE and sensitivity oval HTML renderers have **not** been
connected in this package; this first step produces the general-C5 flavour
observable values and optimal input configs.

### Practical caveats

- `EvaluateGeneralC5Channels.wl` is invoked per numerical evaluation; for
  128+ Sobol points per model this may be costly. Run `--sobol-points 2` for one
  model first, then 8 or 16, before a 128-point all-model study.
- `hbar` is passed as `1/(16*pi^2)` to the numerical coefficient evaluator;
  check the Matchete convention separately when comparing physical masses.
- Additional Yukawas `y1Conj`/`y2Conj` are independent matrices in the RGE
  state, but this pilot **fixes** their initial values as a fraction of their
  ordinary counterpart, and fits only y1/y2.
- Self-bilinear mass splitting is inserted only if the *actual* symbolic C5
  requests split-mass variables. If an allowed extra Yukawa is absent from a
  full-general model's exported C5, that model is blocked, not silently fitted.
- No certified full sequential EFT matching is claimed. This package does not
  prove all general matching diagrams or group factors are captured.
