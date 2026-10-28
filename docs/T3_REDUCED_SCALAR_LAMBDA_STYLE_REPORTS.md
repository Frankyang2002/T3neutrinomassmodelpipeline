# Reduced T3 scalar-sector reports using the original lambda-study renderer

## Correct report family

This rebuild deliberately uses the established **non-general-scalar** presentation interfaces:

- `Numerical/plotting/BuildLambda5e9Displays.py`: `page`, `esc`, `images_as_gallery` for HTML layout and gallery copying. This is the same interface used by the original `lambda_5e-9` presentation pages. It does not itself generate PDFs.
- `Reports/ReportGeneration.py`: `compile_latex_document` for the PDF, following the pipeline's established LaTeX report compilation procedure.

**Not used:** `scripts/BuildGeneralScalar48FinalReports.py` and its ReportLab atlas. `Numerical/plotting/InteractiveModelComparison.py` is also not called, because the frozen-UV saved dataset does not supply the same physical `running_diagnostics.json` curves that its selectors require; in particular, no fictitious Weinberg evolution is generated above a matching threshold.

Historical non-general-scalar run instructions can be found in `docs/REGENERATE_CURRENT_REPORTS.md`, `docs/FERMION_START_REUSE_EXISTING_REPORTS.md`, and `Numerical/plotting/BuildLambda5e9Displays.py`.

## Install (local PowerShell, from repository root)

Place `T3_ReducedScalar_Lambda_Style_Reports.zip` in `C:\thesistemp\T3neutrinomassmodelpipeline`.

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath .\T3_ReducedScalar_Lambda_Style_Reports.zip -DestinationPath . -Force
python -m pytest -q tests\python\test_reduced_scalar_lambda_style_reports.py
python -m scripts.BuildReducedScalarLambdaStyleReports --check
python -m scripts.BuildReducedScalarLambdaStyleReports
```

PDF compilation requires the same `pdflatex` or `latexmk` used by `Reports.ReportGeneration.compile_latex_document`. On Windows, MiKTeX or TeX Live supplies this executable. Check it with:

```powershell
Get-Command pdflatex,latexmk -ErrorAction SilentlyContinue
```

If no LaTeX compiler is installed, use `--html-only` (the `.tex` source is still created). No additional physics dependencies or solver runs are needed.

## Output

`Reports/output/reduced_scalar_lambda_style_v1/`:

- `index.html`: lambda-study-style cross-model overview.
- `optimised.html`: two frozen-UV scenarios, 16 individual fits each.
- `shared_a_m2.html`: same 32 model/scenario evaluations at unchanged A(-2) Yukawas.
- `running_uv.html`: exploratory historical UV mass/hierarchy diagnostics with saved figures.
- `model_A_m2.html` etc.: 16 separate model pages, including the five saved observable pulls in both frozen scenarios.
- `benchmark_report.tex` / `benchmark_report.pdf`: printable tables, interpretation, existing saved PNG figures. PDF is absent in `--html-only` mode.
- `assets/`: copied original images; no figure values modified.
- `provenance.json`: exact file SHA-256 hashes plus renderer/compiler identity.

Open:

```powershell
Start-Process .\Reports\output\reduced_scalar_lambda_style_v1\index.html
Start-Process .\Reports\output\reduced_scalar_lambda_style_v1\benchmark_report.pdf
```

The destination intentionally differs from `Reports/output/display_lambda_5e9/`, `Reports/output/display_lambda_5e10/`, `Reports/output/analytical/`, and the previous `Reports/output/reduced_scalar_final_v2/` adapter. No original report is overwritten. Re-runs should use a **new** output path:

```powershell
python -m scripts.BuildReducedScalarLambdaStyleReports --output Reports/output/reduced_scalar_lambda_style_v2
```

## Input and caveats

Input is exactly `output/studies/reduced_scalar_final_reports/final_reduced_v2/`, especially `frozen_comparison.csv`, `oscillation_observable_pulls.csv`, `running_uv_formal_endpoints.csv`, `manifest.json`, and any saved PNGs. The adapter validates 32 individual fits, 10 acceptable shared A(-2) evaluations, all 320 pull rows, 16 historical running endpoints and the topology-coupling convention (class E is half the A-D input). Invalid/incomplete input fails without rendering.

Primary conclusions: these 18-real-parameter Yukawa fits reproduce five fitted neutrino observables, so a small chi-square is not independent predictive evidence. The shared A(-2) Yukawa transfer frequently fails. UV-running scalar masses beyond the 0.5-5 TeV window represent failure of the *chosen boundary*, not model exclusion; formal nominal 1 TeV EFT1 endpoints are not physical rematched thresholds for heavy scalars. The reduced zero-extra-scalar boundary and existing C5 matching-normalisation caveats remain unchanged.

The script does not change Yukawas, scalar potential, matching artifacts, source datasets, GitHub, or the earlier HTML/PDF reports.
