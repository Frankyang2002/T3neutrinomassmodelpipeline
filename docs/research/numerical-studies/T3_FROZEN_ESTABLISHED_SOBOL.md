# Frozen-UV established Sobol fitting

New local runner: `scripts/RunT3FrozenEstablishedSobol.py`.

This uses the actual repository `Numerical.fitting.SobolBenchmarkSearch.generate_sobol_parameter_points` and `refine_sobol_candidates`, including the original reduced-dimension staged optimisation (5, 8, 12 active parameters), instead of the general-scalar campaign's single-best local refinement. The original conditional general-quartic `RunT3GeneralScalar48Campaign.Evaluator` still provides RG running, C5 and NuFIT residuals. No physics code is changed.

Frozen only: small and large lambdaT3; original scenario lambdaT3 values remain fixed (including E-class half convention). The old source campaign's accepted fits are preserved unless `--refit-accepted` is passed. B(-1)/C(-1) are referenced through their dedicated special matching evaluator and saved validation. Legacy C5 parity failure cases will be explicitly marked `FailedOrUnavailable`; they are not physically excluded. All 32 model/scenario cells appear in summary if required numerical source inputs exist; failures and unavailable entries are separately labelled.

## Execution from repo root (PowerShell)

```powershell
python -X utf8 -m pytest -q tests\test_t3_frozen_established_sobol.py
python -u -X utf8 -m scripts.RunT3FrozenEstablishedSobol `
    --source-campaign overnight_six_v1_bulk48 `
    --run-id frozen_established_v1 `
    --sobol-points 128 `
    --rounds 3 `
    --seed-count 5 `
    *> .\frozen_established_v1.log
```

Rerun the same command after interruption: completed round JSON checkpoints are reused. The search does not overwrite source campaign data. For a focused pilot add `--models A_+0,A_-2`; increase `--rounds` once the pilot is verified. Output: `output/studies/general_scalar_extended/bulk_48/frozen_established_v1/summary.csv`, `summary.json`, and per-case `round_001.json`, `result.json`, `fitted_config.json` when viable. The saved configuration is exploratory, not certified sequential Wilson matching.

Testing done on this package is static/syntax only. Full physics search must run on local repo.
