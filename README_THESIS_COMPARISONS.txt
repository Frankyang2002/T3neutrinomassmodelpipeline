T3 VERIFIED EXISTING COMPARISONS — OCTOBER 24, 2026

Source: GitHub Frankyang2002/T3neutrinomassmodelpipeline main, commit
648a3cb564bd6b0d1dc91ca355e5f542f20491f9

The repository already contains UV-on/fermion-start sensitivity data,
a same-seed A(-2) staged ablation, and the full 16-model 0.5 scalar fit.
Four figures and 3 CSV tables were generated from the saved results.
NO new RGE integration, symbolic calculation, or fitting has been run.

Files in project:
- Numerical/analysis/BuildExistingRGComparisonThesisFigures.py
- docs/THESIS_EXISTING_COMPARISONS_PLACEMENT.md
- docs/THESIS_EXISTING_COMPARISONS_CAPTIONS.tex
- Reports/output/thesis_assets/data/verified_existing_comparison_snapshot.json
- Reports/output/thesis_assets/figures/fig11_* through fig14_* (PNG + PDF)
- Reports/output/thesis_assets/tables/table_existing_scalar_comparison.csv
- Reports/output/thesis_assets/tables/table_uv_sensitivity_maxima.csv
- Reports/output/thesis_assets/tables/table_staged_ablation.csv

Recommended:
Fig11 chapter 5.3 (main text, descriptive UV/fermion comparison)
Fig12 chapter 5.3 or appendix (same-seed, but poor chi2)
Fig13 chapter 5.4 (low vs high scalar fits)
Fig14 chapter 5.5 (representation dependence with shared Yukawas)

Use PowerShell at local repository root:
Expand-Archive -LiteralPath .\T3_VERIFIED_EXISTING_COMPARISONS_FIGURES.zip -DestinationPath . -Force
python -m Numerical.analysis.BuildExistingRGComparisonThesisFigures --source local

For same saved values even if local original outputs unavailable:
python -m Numerical.analysis.BuildExistingRGComparisonThesisFigures --source snapshot

The old thesis guide incorrectly labels both the UV and high-scalar studies
as pending. This package's placement guide supersedes those statements.
