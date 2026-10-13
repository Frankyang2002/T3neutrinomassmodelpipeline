# Numerical comparison scenario correction

Cause: the earlier resume code prioritized `resolved` configurations, which
can contain globally overwritten fitted values. The 64 numerical jobs completed,
but 15/16 model combinations showed scenario-independent results. For one
verified class-A model, all four `resolved` files had `Y1_11=-0.20606645548`
and `lambdaT3=0.001`, whereas `input` files contained the intended benchmark
scenarios.

This version prioritizes `configs/generated_model_sets/comparison/<scenario>/input/`
for every model and writes a separate `numerical_resume` file with sensitivity
disabled. Neither original inputs nor symbolic artifacts are modified.

Run from repository root (PowerShell):

```powershell
python -m pytest -q tests/test_resume_pipeline_numerical_results.py
python -m Numerical.orchestration.ResumePipelineNumericalResults --run --force --limit 1
```

Inspect that first resumed model's `running_diagnostics.json`, especially its
`uv_running.y1_frobenius_norm` initial value, then run:

```powershell
python -m Numerical.orchestration.ResumePipelineNumericalResults --run --force
python -m Numerical.diagnostics.ExtractFullComparisonResults
```

CAUTION: `--force` overwrites numerical `running_diagnostics.json` and figures,
not the Matchete/RGBeta intermediate artifacts. The exported observables
are fixed unfit benchmarks, not guaranteed realistic oscillation fits.
