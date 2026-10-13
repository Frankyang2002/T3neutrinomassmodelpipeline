# Fitted sensitivity comparison

Add `Numerical/orchestration/RunFittedSensitivityComparison.py` to the project.

Prerequisites: existing local `output/full/comparison/smallY_smallL` RGE/Weinberg artifacts and all three fit result directories in `Reports/output/full/` (`fit_feasibility`, `fit_feasibility_retry`, `fit_feasibility_focused`). It chooses the smallest reported chi2 for each model, independently evaluates gamma=1 through the existing model-specific numerical evaluator, then scans T1/T2/T3. No multiplicity factor is added to Class E.

PowerShell from repository root:

```
python -m Numerical.orchestration.RunFittedSensitivityComparison
python -m Numerical.orchestration.RunFittedSensitivityComparison --run --limit 2
python -m Numerical.orchestration.RunFittedSensitivityComparison --run
```

Output: `Reports/output/full/fitted_sensitivity/smallY_smallL/comparison_summary.csv`, per-model `sensitivity_diagnostics.json`, and 10 PNG observable-pair plots under `plots/`. `--force` recomputes the fitted sensitivity outputs only.

The scan uses 11 gamma points in [0.95,1.05] by default, ensuring gamma=1 is sampled. Baseline verification threshold is chi2 <=0.01. This is not a new optimisation and does not certify the coupling-regime labels.
