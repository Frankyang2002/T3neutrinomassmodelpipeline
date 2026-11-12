# Replot A(-4) and D(-2) EFT1 LLSS ablation results — default Matplotlib font

## Project paths

- `scripts/RunT3EFT1LLSSBothDonors.py` — full replacement two-donor report wrapper.
- `tests/test_t3_eft1_llss_both_donors.py` — validation and style tests.
- `README_EFT1_LLSS_BOTH_DONORS_LARGE_FONT.md` — this document.

Requires the existing `scripts/RunT3EFT1LLSSIsolation.py` from the previous package (its `--plot-only` mode). No physics files or previously generated data are altered.

## Data provenance (verified on GitHub)

The published results are located at:

- `output/parallel_t3/eft1_llss_both_donors_end/A_m4/all16_three_way.csv`
- `output/parallel_t3/eft1_llss_both_donors_end/D_m2/all16_three_way.csv`

Each has its own `provenance.json`. The saved D(-2) summary reports 16/16 completed with no blocked models; the A(-4) summary reports 16 valid input models. The wrapper checks all 16 model identities, source donor, and all 3 case observables before plotting.

**No numerical physics is recalculated.** If saved results exist locally in those paths, they are used directly. Otherwise, the wrapper downloads only the two public CSVs and their provenance JSONs from GitHub `main` and checks them before plotting. You can force download with `--github-only`; this does not write to GitHub.

## Run from the project root (PowerShell)

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline

Expand-Archive -LiteralPath "$HOME\Downloads\t3_eft1_llss_both_donors_readable.zip" -DestinationPath . -Force

.\.venv\Scripts\python.exe scripts\RunT3EFT1LLSSBothDonors.py `
  --output output\parallel_t3\eft1_llss_readable_v1
```

The default uses Matplotlib's standard text and math fonts. It produces **title-free** figures at 21 inches wide and 300 DPI with larger labels and annotation text.

## Control font sizes

```powershell
.\.venv\Scripts\python.exe scripts\RunT3EFT1LLSSBothDonors.py `
  --output output\parallel_t3\eft1_llss_readable_v2 `
  --font-size 16 `
  --cell-fontsize 12 `
  --top-tick-fontsize 19 `
  --bottom-tick-fontsize 16 `
  --y-tick-fontsize 17 `
  --axis-label-fontsize 18 `
  --colorbar-fontsize 15 `
  --figure-width 23 `
  --row-height 0.72
```

Use a **new `--output` directory** for each replot. No existing results are overwritten.

### Default sizes

- `--font-size 14` (base for all text)
- axis label 16.1 pt
- bottom `Running` / `Mixing off` tick labels 14 pt
- top observable labels 17.5 pt
- model row labels 15.4 pt
- cell annotations 12.6 pt
- colour bar 14 pt
- width 21 in, row height 0.67 in, DPI 300

Individual flags override their corresponding base sizes; no `font.family` override is made.

## Outputs

- `output/parallel_t3/eft1_llss_readable_v1/A_m4/plots/all16_heatmap_percent.png`
- `output/parallel_t3/eft1_llss_readable_v1/D_m2/plots/all16_heatmap_percent.png`
- Same PDF paths, plus individual mass-squared and mixing-angle PNG/PDF heatmaps.
- Both donor-specific Markdown/LaTeX reports in each `reports/` folder, plus CSV tables.
- `output/parallel_t3/eft1_llss_readable_v1/both_donors_summary.json` documenting sources and file hashes.

The two columns per observable retain the same frozen reference:

\[
\delta_{\mathrm{Running}}=100\frac{O_{\mathrm{full}}-O_{\mathrm{frozen}}}{O_{\mathrm{frozen}}},\qquad
\delta_{\mathrm{Mixing\ off}}=100\frac{O_{\mathrm{LLSS-off}}-O_{\mathrm{frozen}}}{O_{\mathrm{frozen}}}.
\]

The underlying numerical results are conditional exploratory computations, not certified one-loop predictions.
