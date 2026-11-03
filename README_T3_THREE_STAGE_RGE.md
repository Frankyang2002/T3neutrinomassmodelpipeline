# Three-stage RGE display: real SMEFT C5 and neutrino mass running

This additive package provides **one new project script**. It reuses the
already implemented one-loop full-flavour SM+Weinberg beta functions and the
saved general 16-model shared-A(-2) gamma=1 centre data. It does not change
UV/EFT1 running, Matchete matching, symbolic beta functions, fitted Yukawas,
existing sensitivity reports, or any GitHub resources.

## Install / file locations

Extract the ZIP in your local repository root (preserve paths):

- `scripts/BuildT3SharedA2ThreeStageRGE.py` **new** report adapter.
- `tests/python/test_t3_three_stage_rge.py` **new** offline tests.
- `README_T3_THREE_STAGE_RGE.md` **new** documentation.

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive "$HOME\Downloads\T3_Three_Stage_RGE_Display.zip" -DestinationPath . -Force
python -X utf8 -m pytest -q tests/python/test_t3_three_stage_rge.py
```

## Run with the saved general-Lagrangian gamma=1 centres

```powershell
python -X utf8 -m scripts.BuildT3SharedA2ThreeStageRGE `
  --run-id shared_A2_general16_gamma08_12_v1 `
  --points 41

Start-Process "output\studies\general_scalar_extended\shared_a2_general16_reports\shared_A2_general16_gamma08_12_v1\interactive_RGE_with_SMEFT.html"
```

If the shared run has not yet been generated locally, first use the **original**
21-point general sensitivity and UV/EFT1 RGE pipeline, with existing data and
Matchete source present:

```powershell
python -X utf8 -m scripts.RunT3SharedA2General16Reports `
  --donor-run-id thesis16_general_c5_full_v1 `
  --run-id shared_A2_general16_gamma08_12_v1 `
  --gamma-min 0.8 --gamma-max 1.2 --gamma-points 21 `
  --rge-points 41
```

If you already have gamma diagnostics and centres but need the 14 genuine UV/EFT1
curves, run the existing script with `--rge-only --rge-points 41` using the
same donor/run IDs instead of repeating the full sensitivity scan.

## What it actually generates

Inside the saved run directory:

- `interactive_RGE_with_SMEFT.html`: offline interactive selectors for the
  **actual existing UV and EFT1 solver samples** and **new numerically
  integrated SMEFT samples** for the 14 nonmixed models.
- `<model>/smeft_trajectories.json`: 41 full-flavour C5/neutrino mass samples,
  plus solver and saved-endpoint consistency diagnostics.
- `<model>/smeft_observables.csv`: scale, C5 norm and selected elements,
  `m1,m2,m3`, mass-squared differences, `sin2_theta12/13/23`, running SM
  gauge couplings and Higgs quartic.
- `smeft_status.json`: explicit success/failure provenance.

The original `interactive_RGE.html`, `sensitivity_oval.html`, `centre.json`,
`rge_trajectories.json`, gamma scan and Matchete source outputs remain intact.
If any individual SMEFT run disagrees with its saved `C5_low` or low-energy
oscillation prediction, it is **excluded** and an error JSON is written.

## Precise physics limitations

1. The saved **general** C5 boundary is from **conditional one-scale Matchete
   matching**, not a certified full sequence of integrating out all fermions,
   then all scalars. The report labels this explicitly.
2. `C5` in EFT1 is **not** produced by this adapter. The physical LLSS
   operator closure, normalisation, and general antilinear `C*` anomalous
   dimensions are still not certified in the current general solver. The
   dashboard displays this missing stage instead of inserting unverified data.
3. Mixed B(-1)/C(-1) do not have compatible saved dense UV/EFT1 curves and
   are omitted, as in the current interactive-RGE report.
4. In SMEFT, `mν(μ) = -v² C5(μ)/2` uses constant `v=246.22 GeV` and running
   charged-lepton basis, matching the existing general fit. These running
   mass parameters are not physical pole masses.
5. The **Higgs quadratic mass-squared** `mH²` is not in the existing
   `SMWeinbergEvolution` state. This report shows **neutrino** mass parameters
   and the Higgs quartic `lambdaH`, not a fictitious `mH²` running.

To obtain the missing EFT1 C5 curves legitimately, the next source task is to
connect the matched LLSS Wilson tensor to the complete physical EFT1 beta
system and evolve every sourced and zero-boundary component, including the
complex-conjugate terms. Only then should this stage enter the numerical
comparison or thesis claims.
