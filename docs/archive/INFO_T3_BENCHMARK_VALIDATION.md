# T3 multiplicity-only benchmark validation

This is the runbook for the charge-multiplicity-normalized 16-model comparison implemented by `studies/FullT3Study.py`.

## Exact benchmark rule

The four fixed comparison scenarios use `Y=0.005` or `0.5`, and `lambdaT3_reference=0.1` or `1.0`. The revised factors are:

```python
LAMBDA_T3_NORMALISATION = {"A": 1.0, "B": 1.0, "C": 1.0, "D": 1.0, "E": 0.5}
```

Thus `lambdaT3_effective` equals 0.1 or 1 for A–D, and 0.05 or 0.5 for E. This preserves relative representation CG weights and only compensates two equal-weight, equal-loop-integral class-E charge channels in the unbroken-SU(2) benchmark approximation.

## Commands (PowerShell, project root)

```powershell
python -m pytest -q tests/test_t3_multiplicity_benchmark.py
python -m pytest -q tests/test_full_study_architecture.py
python -m pytest -q
python pipeline.py --full
```

The exact full-study option should be confirmed by `python pipeline.py --help` if the local CLI differs. The full study invokes external symbolic/RGE stages and can be long-running; a passing local unit test does not establish 64-case end-to-end success.

## Output and ownership

The comparison generator writes per-model input files below `configs/generated_model_sets/comparison/<scenario>/input/` and the summary to `output/full/full_numerical_summary.json`. The summary includes the `lambdaT3_normalisation` basis/factors, and each configuration records `comparison_point.lambdaT3_reference`, `.lambdaT3_effective`, `.normalisation_factor`, and `.normalisation_basis`.

Only `studies/FullT3Study.py` benchmark initialization changes. Do not modify `RGE/`, `Lagrangian/`, `Numerical/running/`, or the physical group coefficients as part of this update. Existing generated output/configuration files are historical results until regenerated.

## Status

Code-level check: syntax and factor tests on the replacement files can be run without Matchete. Full numerical outputs for the revised 64 benchmark points remain **pending user-side execution**.
