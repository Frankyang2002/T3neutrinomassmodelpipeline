# Installation — Chat A Batch 3

**Additive** files at project-relative paths:

- `Numerical/lo_patterns_16/analytic_explanation.py`
- `Numerical/lo_patterns_16/BATCH03_METHOD.md`
- `Numerical/lo_patterns_16/BATCH03_INSTALL.md`
- `scripts/ExplainLOPatterns16.py`
- `tests/test_lo_analytic_explanation_16.py`

Requires the previous `Numerical/lo_patterns_16/result_analysis.py` from Batch 2, the installed `Numerical/lo_certification/direct_uv` and `Numerical/direct_one_loop` packages, and complete results at `provisional_01` and `analysis_batch02_01`.

From `C:\thesistemp\T3neutrinomassmodelpipeline`:

```powershell
$repo='C:\thesistemp\T3neutrinomassmodelpipeline'
$zip=Join-Path $HOME 'Downloads\Chat_A_Batch03_Analytical_Explanation.zip'
Expand-Archive -LiteralPath $zip -DestinationPath $repo
Set-Location $repo

git rev-parse HEAD
.\.venv\Scripts\python.exe -m pytest -q tests\test_lo_analytic_explanation_16.py
.\.venv\Scripts\python.exe scripts\ExplainLOPatterns16.py
```

A complete run produces a NEW `output/general_observable/lo_pattern_study16/analysis_batch03_01` directory containing:

- `BATCH03_ANALYTICAL_REPORT.md` — all-16 equations/numerical interpretation.
- `batch03_analytic_explanation.json` — all-16 conditional model summaries.
- `batch03_model_patterns.csv` — three masses, hierarchy, loop-kernel residuals.
- `batch03_charge_channels.csv` — each source charge, scalar mixing and exact matrix norm.
- `batch03_matched_shapecases.csv` — comparisons with one best complex scale factor.
- `models/<model>.json` — original matrix, independently reconstructed exact and small-mixing matrices, full charge matrices and heavy-mass group data.
- `SOURCE_MANIFEST.json`, `HANDOFF.md`, `HANDOFF.json`.

All original inputs and previous outputs are read-only and SHA-checked. The script fails closed on any unsupported source drift, additional operators, failed identity or incorrect certification claim. **These are provisional numerical patterns, not a certified physical 4-spinor-to-Weyl mapping.**
