# T3 fixed-comparison numerical-only recovery

## Verified observations (local run, October 2026)

The 64-case `python pipeline.py --full` evaluation encountered two separate orchestration failures, **not evidence of an incorrect CG or RGE result**:

1. The default single threshold `(F,S1,S2)` produced no intermediate scalar EFT; the integrated numerical continuation requires `eft1_rgbeta_rge.json`. For symbolic runs requiring this intermediate EFT, pass `--threshold F --threshold S1 S2`.
2. The numerical figure-output implementation originally used `Path(__file__).resolve().parents[1]` in `Numerical/orchestration/PipelineNumericalResults.py`, pointing to `Numerical/` instead of the repository root. The local fix uses `parents[2]` and was verified to output under `Reports/output/...`.

One existing benchmark (`largeY_largeL`, class E, dimensions `(3,3,2)`, alpha `+2`) subsequently ran through `run_pipeline_numerical_results` successfully, writing `data/running_diagnostics.json` and figures, **without Matchete/RGBeta reruns**. This is one verified numerical case; the other cases have not been counted as numerical successes.

## Resume command

`Numerical/orchestration/ResumePipelineNumericalResults.py` discovers pre-existing directories beneath `output/full/comparison/<scenario>/<model-key>/<physical-model>/`. It requires three files in each physical-model `data/` directory:

- `uv_rgbeta_rge.json`
- `eft1_rgbeta_rge.json`
- `final_weinberg_coefficient.json`

For configuration, it prefers that scenario's `configs/generated_model_sets/comparison/<scenario>/resolved/<key>.json`, then its `input/<key>.json`. It does **not** fall back to the global `Numerical/configs/generated_models/<key>.json`: that file is shared across scenarios and is overwritten by successive benchmark runs.

The resume runner refuses configs that enable benchmark search, skips existing successful `running_diagnostics.json` unless `--force` is specified, and records individual failures without stopping other cases. It does not update the original `output/full/full_numerical_summary.json`, which describes the prior full pipeline run.

```powershell
python -m Numerical.orchestration.ResumePipelineNumericalResults
python -m Numerical.orchestration.ResumePipelineNumericalResults --run --limit 1
python -m Numerical.orchestration.ResumePipelineNumericalResults --run
```

The first command is dry-run inventory; the second runs at most one unfinished ready case. The third processes all ready cases. Summaries are written to `Reports/output/full/numerical_resume_summary.json` and `.csv`.

## Limitations

A successful numerical resume verifies this numerical continuation under its fixed scenario configuration. It does not repeat symbolic matching, recompute RGEs, prove phenomenological agreement, or turn previously failed symbolic cases into successes. Broken-phase mass splitting remains outside the multiplicity-only benchmark approximation.

## October 2026 resume retry: NuFIT sensitivity path

A batch resume attempted 63 previously unfinished cases, and all 63 stopped because
an **optional sensitivity scan** tried to read
`Numerical/configs/nufit_6_1_ic24_no.json` rather than the project's
`configs/nufit_6_1_ic24_no.json`. This is an incorrect relative-root
resolution in `PipelineNumericalResults.py`'s sensitivity branch
(`parents[1]` rather than repository-root `parents[2]`).

The four fixed comparison scenarios are not intended to conduct sensitivity
scans. The resume utility therefore creates a separate, persistent effective
configuration for each attempted case at:

`configs/generated_model_sets/comparison/<scenario>/numerical_resume/<model-key>.json`

It copies the selected input/resolved benchmark and forces
`"sensitivity": {"enabled": false, ...}` while keeping all other numerical
parameters, model representation, benchmark-search settings, and fixed Yukawa
and lambdaT3 values. It does not mutate the original scenario configs and
does not re-run Wolfram/Matchete/RGBeta. The effective config path is recorded
in the batch summary under `resume_config` and in running diagnostics.

If a separate sensitivity study is requested in the future, fix the standalone
sensitivity path in `PipelineNumericalResults.py` rather than depending on this
comparison-specific bypass. This package does not change that independent
sensitivity implementation.
