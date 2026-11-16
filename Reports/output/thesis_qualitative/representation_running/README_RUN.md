# Representation versus running — local read-only analysis

The bundled `scripts/AnalyzeT3RepresentationRunning.py` is a complete replacement/new analysis file, not a diff. Install it at exactly `scripts/AnalyzeT3RepresentationRunning.py` in the repository root. It reads existing data and writes *only* into `Reports/output/thesis_qualitative/representation_running/`.

Required inputs are existing repository files under `output/studies/general_scalar_extended/` plus `Numerical/configs/generated_models/`. **Do not copy the included presentation CSVs into any physics-input directory.**

PowerShell commands (with current directory at the local repository root):

```powershell
# Optional: confirm the local version matches the GitHub main SHA reviewed.
git rev-parse HEAD
# Expected: 6ebe9937a72880c6f180bcd026261779800ac71f

# Extract the downloaded ZIP from the folder in which it was saved.
Expand-Archive -LiteralPath '.\T3_representation_running_analysis.zip' -DestinationPath '.' -Force

# Only plotting/numerical packages; no matching/RGE/fitting dependencies.
python -m pip install numpy matplotlib
python .\scripts\AnalyzeT3RepresentationRunning.py

# Review generated files.
Get-ChildItem '.\Reports\output\thesis_qualitative\representation_running' -Recurse
```

Reproduced from existing local source files, the script will generate:

- `csv/all_256_saved_cases.csv` — U/E/W/S switches, full masses, mixings (sin² and degrees), splittings and C5-log diagnostic.
- `csv/full_1111_baselines.csv` — fitted 1111 neutrino spectra.
- `csv/frozen_stage_observable_changes.csv` — signed absolute changes and percentages when reference values are meaningful; explicitly tags Wilson as conditional.
- `csv/representation_stage_heatmap_values.csv` — quantitative model × running-stage response.
- `csv/model_family_stage_summary.csv` — within-family min/median/max and sign counts.
- `csv/conditional_wilson_log_diagnostics.csv` — distinct from independently evolved stages.
- `csv/compensated_coupling_sensitivity_slopes.csv` — centered differences of **already-saved** sensitivity points, not refits.
- `plots/independent_stage_mass_heatmap.png` — 16 × 3 independently evolved effects.
- `plots/conditional_wilson_toggle.png` — separate conditional Wilson sensitivity.
- `plots/observable_stage_changes.png` — selected oscillation-observable responses.
- `physical_interpretation.md` and `data_validation.json` — refreshed from full local data.

The ZIP already contains three *preview* plots, two preview CSVs, and an interpretation calculated from the checked GitHub saved outputs, rounded for presentation. The script produces additional, full-precision outputs after it reads the existing **local** source data. It does not contact GitHub and requires no Mathematica or matching/RGE reruns.

**Important:** W-off is `C5_one_scale - LLSS_log` bookkeeping at a fixed one-scale anchor, not an independently hard-matched physical prediction. All accepted fit targets are calibrated to the same oscillation data. The conditional B(-1)/C(-1) results use channel normalization; the independent audit validates selected mass-log coefficients, not complete hard matching.
