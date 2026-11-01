# Direct B(-1)/C(-1) exploratory general-scalar numerical pilot

This **new launcher** avoids the missing historical `small_lambda_1e-7` and
`lambda_5e-9_refit` outputs. It loads the **tracked** `Numerical/benchmarks/fitted_reference`
configuration for B(-1) or C(-1), keeps the full candidate B/C general-scalar
quartic and mixed-bilinear beta-function state, and reuses the repository's
existing Matchete ordered one-scale C5 evaluator, NuFIT numerical evaluation,
Sobol sampler, and staged least squares.

It does **not** edit any existing physics script or GitHub. The previous
`RunT3GeneralScalarBCExploratoryFit.py` is left untouched.

## Install

Extract the ZIP at the project root, so this new file is:

`C:\thesistemp\T3neutrinomassmodelpipeline\scripts\RunT3BCGeneralFromReference.py`

Activate the existing venv, then from the project root:

```powershell
Get-Command wolframscript
python -X utf8 -m scripts.RunT3BCGeneralFromReference --help
```

## First quick numerical input -> output check

```powershell
python -X utf8 -m scripts.RunT3BCGeneralFromReference `
  --model B_-1 --run-id thesis_b_m1_smoke `
  --export-betas --sobol-points 2 --no-local-refinement
```

`--export-betas` runs the existing `RunGeneralT3QuadraticEFT1.export_live` helper
for **both UV and EFT1** using locally installed Wolfram/RGBeta. It does NOT
require the old Stage33 trajectories. Once betas are exported successfully,
subsequent experiments can omit `--export-betas`.

## Proper Sobol + local least squares

```powershell
python -X utf8 -m scripts.RunT3BCGeneralFromReference `
  --model B_-1 --run-id thesis_b_m1_128 `
  --sobol-points 128 --max-nfev 60
```

Repeat with `--model C_-1` and a new `--run-id` (the C(-1) beta exports must
also be generated once using `--export-betas`).

To change the fixed topology quartic, use `--lambda-t3 0.05` or similar.
The tracked initial values have `lambdaT3=0.1`. If changed, the Yukawa matrices
are rescaled by `sqrt(0.1/lambdaT3)` **as an initial guess only**, then refitted.
Other additional quartic boundaries are chosen by `--quartic-fraction 0.1`;
the allowed B/C scalar mixing input by `--mu2-fraction 0.01`.

## Outputs

`output/studies/general_scalar_extended/fresh_reference_bc/<run-id>/B_m1/`

- `fitted_config.json`: optimised numerical input config
- `result.json`: chi2, five observables, C5 matrices, parameters, metadata
- `observables.csv`: NuFIT target versus fitted predictions

No historical fitted files are synthesized. No old outputs are overwritten.
This pilot intentionally does not generate interactive RGE or oval HTML; the
existing HTML generators require time-sampled trajectories/sensitivity scans,
which are a separate downstream integration.

## Explicit scope

This is an **exploratory B/C-only** workflow: existing candidate full scalar
quartic + mixed quadratic RGBeta UV/EFT1 running, and an existing Matchete
**one-scale** C5 insertion with the conditional `lambda_M = -lambda_RGBeta`
map. It is not certified general sequential EFT matching. The B(-1)/C(-1)
models do not have additional conjugate Yukawas in this catalogue. This does
not yet provide 16-model complete-general fits.

If the run reports missing Wolfram/RGBeta libraries or a failed beta export,
that is a local dependency or unsupported source issue. Do not substitute
restricted RGEs for those general beta files.
