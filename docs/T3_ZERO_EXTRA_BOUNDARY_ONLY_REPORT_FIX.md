# Zero-additional-scalar: boundary_zero-only RGE and sensitivity

This package corrects the previous *report adapter* and is **local-only**.
It does not modify the T3 physics kernels, reference HTML canvas, sensitivity
renderer, or the saved 16-model zero-extra control results.

## Files (full replacements)

- `scripts/GenerateT3ZeroExtraGamma.py` — adds `--mode` (default `boundary_zero`); fixes B/C endpoint lookup by passing the individual model directory `.../boundary_zero/B_m1` or `.../boundary_zero/C_m1` to the existing gamma calculator.
- `Numerical/plotting/BuildT3ZeroExtraReferenceReports.py` — adds `--mode` (default `boundary_zero`); the selected mode alone is rendered and validated. Incomplete `held_zero` gamma scans are not read.
- `scripts/RunT3ZeroExtraReferenceReports.ps1` — defaults to `boundary_zero`; opens `boundary_zero/interactive_rge.html`.
- `tests/test_t3_zero_extra_reference_reports.py` — complete test file retaining the five original tests and adding two mode/path-regression tests.

The source study remains:

`output/studies/general_scalar_extended/zero_extra_all16/zero_extra_all16_low_v1/`

The gamma grid remains `0.95,0.96,...,1.05`, in directions T1, T2, T3.
For boundary-only the required new scan is **16 × 3 × 11 = 528 calculations**.
The 32 validated 1 TeV–100 GeV Weinberg trajectories reported by the user are reused.

## PowerShell (from repository root)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\T3_zero_extra_boundary_only_fix.zip" -DestinationPath . -Force
python -X utf8 -m pytest -q tests\test_t3_zero_extra_reference_reports.py
Unblock-File .\scripts\RunT3ZeroExtraReferenceReports.ps1
.\scripts\RunT3ZeroExtraReferenceReports.ps1 -Mode boundary_zero
```

If your original ZIP is elsewhere, substitute its real path.

Alternatively, because the Weinberg running finished in the user's log,
run only the two remaining commands:

```powershell
python -u -X utf8 -m scripts.GenerateT3ZeroExtraGamma --run-id zero_extra_all16_low_v1 --scan-id gamma_095_105_v1 --mode boundary_zero --resume
python -X utf8 -m Numerical.plotting.BuildT3ZeroExtraReferenceReports --run-id zero_extra_all16_low_v1 --weinberg-output-id weinberg_running_v1 --scan-id gamma_095_105_v1 --mode boundary_zero
Start-Process .\Reports\output\display_zero_extra_all16\zero_extra_all16_low_v1\boundary_zero\interactive_rge.html
```

**Do not delete** `gamma_095_105_v1/held_zero/B_m1.csv` or rerun the `held_zero`
mode. Its previous failed scan is left untouched and is ignored by this run.

The only written report mode is:

`Reports/output/display_zero_extra_all16/zero_extra_all16_low_v1/boundary_zero/`

with `interactive_rge.html`, `sensitivity_ovals.html`, and `all_results.html`.
An index and source-provenance report are also written in the enclosing directory.
The report uses the exact `optimised_rge.html`-derived canvas from the existing
`BuildT3FullScalarFrozen64HTML` adapter. Full sequential Wilson matching and
intermediate `C12(mu)` numerical integration are not certified/saved.

## Verified vs pending

- Verified by source inspection: old B/C matching expects model-directory-specific `*_eft1_quadratic_running.json`, whereas the previous gamma adapter passed its parent mode directory.
- Verified locally: seven isolated tests pass, and the replacement Python files compile.
- Not yet verified: the 528 real numerical gamma evaluations in the user's Windows checkout, including B(-1)/C(-1).
