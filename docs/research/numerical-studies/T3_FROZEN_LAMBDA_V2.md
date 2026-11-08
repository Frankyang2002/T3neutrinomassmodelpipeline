# Frozen Lambda report v2

This update preserves the original repository HTML/canvas generators (`Numerical.plotting.InteractiveModelComparison` and `Numerical.plotting.BuildFittedModelEnvelopes`). No original files are changed.

**Implemented:** NuFIT Gaussian overlays and JUNO six-year projection are restored using source values extracted from the user-supplied `envelopes.html`; accepted frozen fit chi2 labels are shown, and a 32-combination coverage audit is generated.

**Not implemented:** historical fitted-lambda and shared A(-2) new numerical evaluation. They are identified separately rather than claiming their saved data are frozen general-scalar results. Full C5(mu)/PMNS running cannot be inferred from the saved scalar EFT1 trajectories.

**Benchmarks:** small and large frozen lambda values are included when their accepted state files exist. The historical Lambda benchmark is a separate comparison, not a third frozen fitting case.

## PowerShell

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
python -X utf8 -m pytest -q tests\test_frozen_display_pipeline.py
python -u -X utf8 -m scripts.GenerateT3FrozenRGEData --run-id rge_v1 --run-ids frozen_massnorm_a0_v1,frozen_massnorm_remaining_v1 --source-campaign overnight_six_v1_bulk48
python -u -X utf8 -m scripts.GenerateT3FrozenSensitivityScans --run-id sensitivity_v1 --run-ids frozen_massnorm_a0_v1,frozen_massnorm_remaining_v1 --source-campaign overnight_six_v1_bulk48
python -X utf8 -m scripts.BuildT3FrozenLambdaNative --run-ids frozen_massnorm_a0_v1,frozen_massnorm_remaining_v1 --output-dir frozen_lambda_native_v2
Start-Process .\Reports\output\display_frozen_general_scalar\frozen_lambda_native_v2\index.html
```

`coverage.csv` lists all 16 models for each frozen case and whether fit/RGE/sensitivity data are actually present. Historical shared A(-2) is not relabelled as new frozen data. For full comparison, identify and evaluate the saved historical shared A(-2) benchmarks against the frozen general-scalar evaluator.
