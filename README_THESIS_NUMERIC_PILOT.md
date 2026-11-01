# T3 thesis numerical pilot — original physics/plotting routines

Destination paths in `C:\thesistemp\T3neutrinomassmodelpipeline`:

- `scripts/PrepareThesisNumericPilot.py` — NEW input-config builder and optional runner.
- `scripts/BuildThesisNumericHTML.py` — NEW presentation adapter using the ORIGINAL `InteractiveModelComparison` and `BuildFittedModelEnvelopes` renderers.
- `tests/python/test_thesis_numeric_pilot.py` — offline formatting/control-flow regression tests.

**No existing repository file is replaced. No GitHub modification.** No Matchete/RGBeta/ODE/matching/fit equations are rewritten.

## Scope

The original numerical pipeline uses the **restricted** T3 threshold-matching convention. The 16 fitted-reference configs under `Numerical/benchmarks/fitted_reference` were fitted to NuFIT, so reproducing their five observables is calibration, **not** independent prediction. Do not describe these as full general-interaction neutrino-mass predictions.

## Commands (PowerShell, from repository root)

```powershell
$env:PYTHONUTF8 = '1'
python -X utf8 -m pytest -q tests/python/test_thesis_numeric_pilot.py

# Copy the committed A(-2) fitted-reference input, keep the actual fitted Yukawa matrices,
# and enable a local 0.8-to-1.2 gamma sensitivity scan (11 points).
python -X utf8 -m scripts.PrepareThesisNumericPilot --model A_-2

# See/read and optionally modify this JSON before running:
# configs/thesis_numerics/T3_dS1_1_dS2_3_dF_2_alpha_m2_reference.json

# Execute existing restricted matching, RGE evolution, observables, sensitivity scan.
# Reuses the config created in the preceding step. The pipeline resets the model-specific generated config; save a copy of that existing file first if you want to retain it. The wrapper `--run` mode automatically saves such a backup before resetting.
python -X utf8 pipeline.py --dims 1 3 2 --alpha -2 --restricted-production `
  --threshold F --threshold S1 S2 --study thesis_A_m2 `
  --numerical configs/thesis_numerics/T3_dS1_1_dS2_3_dF_2_alpha_m2_reference.json `
  --reset-numerical-configs

# Generate ORIGINAL interactive plot style and save real computed observables.
python -X utf8 -m scripts.BuildThesisNumericHTML --study thesis_A_m2
```

Saved outputs expected (when numeric integration succeeded):

- `Reports/output/thesis_A_m2/Numerical/interactive_RGE.html` — original interactive RGE chart.
- `Reports/output/thesis_A_m2/Numerical/sensitivity_oval.html` — original empirical sensitivity-oval chart, **not** an experimental confidence region; absent if no successful scan.
- `Reports/output/thesis_A_m2/Numerical/observables.csv` — actual computed low-scale masses and five observables, with baseline chi2 if the sensitivity scan succeeds.
- `Reports/output/thesis_A_m2/Numerical/observables.json` — same data.
- `Reports/output/thesis_A_m2/Numerical/README.md` — links and scope.

Open:
```powershell
Start-Process Reports/output/thesis_A_m2/Numerical/interactive_RGE.html
Start-Process Reports/output/thesis_A_m2/Numerical/sensitivity_oval.html
Import-Csv Reports/output/thesis_A_m2/Numerical/observables.csv | Format-List
```

## Optional exploratory inputs/refit

Use a new output path (preparation refuses to overwrite an existing config):

```powershell
python -X utf8 -m scripts.PrepareThesisNumericPilot --model A_-2 --lambda-t3 0.05 `
  --yukawa-factor 1.0 --output configs/thesis_numerics/A_m2_lambda005.json
```

To **refit** rather than reuse an already-optimised reference:

```powershell
python -X utf8 -m scripts.PrepareThesisNumericPilot --model A_-2 --refit `
  --sobol-points 64 --output configs/thesis_numerics/A_m2_refit.json
```

Then run `pipeline.py` with that new config and a distinct `--study`, e.g. `thesis_A_m2_refit`. The existing Sobol/local refiner may require many numerical evaluations. The refit config searches a *local neighbourhood* of the reference; it is not a guarantee of a global optimum.

## Inputs and proposed values

The starter copies values from `Numerical/benchmarks/fitted_reference/T3_dS1_1_dS2_3_dF_2_alpha_m2.json`. They use the thesis hierarchy 1e7 -> 1e5 -> 1e3 -> 1e2 GeV, MF=(100,110,120) TeV, lambdaT3=0.1 (A-D), 0.05 (E), fixed SM gauge values, and model-dependent fitted y1/y2 with typical amplitudes O(1e-4). The previously-used generic default configuration starts from diagonal O(0.05) Yukawas instead and should not replace the fitted seed for this pilot.

**Reminder:** thesis Appendix C's lambdaT3=0.5 frozen-UV 16-model exercise is a different comparison scenario, not the starting fit used here.

The `--run` launcher automatically backs up an existing `configs/generated_models/<model>.json` to `configs/thesis_numerics` before invoking the original pipeline.
