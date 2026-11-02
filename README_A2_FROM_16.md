# A(-2) fitted-Yukawa benchmark from the existing 16-model run

This package provides a full replacement of `scripts/RunT3Aminus2YukawaBenchmark.py` plus an offline regression test. It does **not** change the main 16-model fitter, the Lagrangian, or the general Matchete coefficient.

## Source files

Expected at project root `C:\thesistemp\T3neutrinomassmodelpipeline`:

- `output/studies/general_scalar_extended/thesis_general_all16/thesis16_general_c5_full_v1/A_m2/result.json`
- `output/studies/general_scalar_extended/thesis_general_all16/thesis16_general_c5_full_v1/A_m2/fitted_config.json`
- `output/studies/general_scalar_extended/thesis_general_all16/thesis16_general_c5_full_v1/settings.json`
- `output/general_matchete_matching/A_m2/c5_coefficient.txt`

The benchmark reads 18 fitted **real** UV Yukawa entries (`y1`,`y2`) from the A(-2) result and checks that the entries in `fitted_config.json` match the optimizer output. It also checks the `c5_coefficient.txt` SHA256 recorded at fit time and checks the source campaign used general Matchete matching. **If A(-2) is unfinished or files have changed, it aborts rather than use the original reference seed.**

The independent Yukawa `y2Conj` was **not fitted** by the campaign. Its nominal boundary remains `0.05*y2_reference`, and the benchmark varies it over `0`, `0.05*y2_reference`, `0.10*y2_reference`, holding the fitted `y1`,`y2` fixed and the self-bilinear nonzero. It keeps the fitter's original convention, **not** `y2Conj=f*y2_fitted`.

## Install and run in PowerShell

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive "$HOME\Downloads\T3_Aminus2_From_16_Fit.zip" -DestinationPath . -Force
python -X utf8 -m pytest -q tests/python/test_thesis_a2_yukawa_benchmark.py tests/python/test_thesis_a2_yukawa_from_16.py

$base = "output\studies\general_scalar_extended\thesis_general_all16\thesis16_general_c5_full_v1\A_m2"
Test-Path "$base\result.json"
Test-Path "$base\fitted_config.json"
if (Test-Path "$base\result.json") {
    $r = Get-Content "$base\result.json" -Raw | ConvertFrom-Json
    $r | Select-Object model,status,chi2_final,best_source
}

python -X utf8 -m scripts.RunT3Aminus2YukawaBenchmark `
    --fit-run-id thesis16_general_c5_full_v1 `
    --benchmark-id a_m2_from_16_general_v1

Import-Csv "output\studies\general_scalar_extended\a_m2_yukawa_benchmarks\a_m2_from_16_general_v1\summary.csv" | Format-Table -AutoSize
```

Benchmark outputs: `output/studies/general_scalar_extended/a_m2_yukawa_benchmarks/a_m2_from_16_general_v1/` including `summary.csv`, `yukawa_matrices.csv`, `manifest.json`, and per-case JSON with C5 and neutrino observables. The manifest includes the exact fit result and fitted-config paths and SHA256 values.

The benchmark is a one-scale conditional general Matchete-C5 evaluation, not a validated sequential-EFT prediction. It does not re-optimise Yukawas for each `y2Conj` fraction.
