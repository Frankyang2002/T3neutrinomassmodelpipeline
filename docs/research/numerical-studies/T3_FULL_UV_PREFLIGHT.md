# Full UV running: scalar-screen preflight

This tool reuses `scripts/RunT3GeneralScalar48Campaign.py` and saved RGBeta exports; it makes **no changes** to mass restrictions, beta functions or C5 matching.

## Why this is necessary

Under full UV running the scalar masses and Yukawas evolve together. Previous searches reported scalar threshold masses ~38–42 TeV, above the screen 0.5–5 TeV. Mass normalisation alone does not fix this. This script evaluates only the full UV and EFT1 running and reports the physical screen, saving one record per model. It performs **no neutrino fit**.

The `--frozen-run` seed is an optional exploratory cross-scenario input; it should not be interpreted as an accepted running-UV configuration.

## Run

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
python -X utf8 -m pytest -q tests\test_t3_fulluv_running_preflight.py
python -u -X utf8 -m scripts.ProbeT3FullUVRunning `
  --run-id fulluv_preflight_v1 `
  --source-campaign overnight_six_v1_bulk48 `
  --models A_+0,A_-2,D_+0,E_-2 `
  --frozen-run frozen_massnorm_a0_v1
```

To resume after an interrupted run, add `--resume`. Report: `output/studies/general_scalar_extended/bulk_48/fulluv_preflight_v1/fulluv_preflight.jsonl`.

If no seeds pass, the next step is a **scalar UV boundary search** (mass inputs/quartics), not more mass-normalized Yukawa Sobol. If any seed passes, run the authoritative C5/mixing evaluator on it before fitting its Yukawas. The normalization law used in frozen mode needs separate verification because full UV running modifies scalar masses as Yukawas change.
