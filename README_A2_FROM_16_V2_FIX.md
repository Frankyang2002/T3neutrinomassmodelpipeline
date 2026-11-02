# A(-2) benchmark: fix 16-model campaign provenance check

The running 16-model campaign writes `matching_scheme = "full-general-Matchete-Yukawa-channel-massbasis-v2"` in its `settings.json`. The preceding benchmark script mistakenly required `"full-general-Matchete-Yukawa-channel-v1"`, rejecting the correct campaign before any numerical evaluation.

This package replaces only:

- `scripts/RunT3Aminus2YukawaBenchmark.py` — now requires the exact **v2 mass-basis** matching scheme and continues to validate the saved best-fit Yukawa values against `fitted_config.json` and the SHA256 of the general Matchete C5 coefficient.
- `tests/python/test_thesis_a2_yukawa_from_16.py` — updates the simulated 16-model settings to v2 and adds regression coverage rejecting v1.

No general Lagrangian, Matchete, RGEs, C5 formula or 16-model optimisation files are modified. The benchmark keeps the nonzero self-bilinear and varies the separate `y2Conj` boundary while reusing fitted A(-2) `Y1` and `Y2`.

From `C:\thesistemp\T3neutrinomassmodelpipeline`:

```powershell
Expand-Archive "$HOME\Downloads\T3_Aminus2_From_16_V2_Fix.zip" -DestinationPath . -Force
python -X utf8 -m pytest -q tests/python/test_thesis_a2_yukawa_benchmark.py tests/python/test_thesis_a2_yukawa_from_16.py
python -X utf8 -m scripts.RunT3Aminus2YukawaBenchmark `
  --fit-run-id thesis16_general_c5_full_v1 `
  --benchmark-id a_m2_from_16_general_v1
Import-Csv "output\studies\general_scalar_extended\a_m2_yukawa_benchmarks\a_m2_from_16_general_v1\summary.csv" |
  Format-Table -AutoSize
```

The script requires the A(-2) result and fitted configuration to have been saved by the completed model step; the *overall* 16-model campaign may still be running. The three benchmark cases have fixed fitted Y1/Y2 and `y2Conj` equal to 0, 0.05, 0.10 times the **reference** Y2. It does not refit independently and remains a conditional general one-scale matching calculation, not certified sequential EFT matching.
