# Original simultaneous-EFT1 LLSS ablation — split thesis heatmaps

Intended project paths:

- `scripts/RunT3OriginalSimultaneousLLSSHeatmap.py`
- `tests/test_t3_original_simultaneous_llss.py`
- `README_ORIGINAL_SIMULTANEOUS_LLSS.md`

## What changed

This keeps the original simultaneous-EFT1 LLSS calculation unchanged, but extends the existing report generation.

Outputs now include:

- combined 16 x 10 heatmap: `plots/all16_heatmap_percent.png`
- mass-squared-only heatmap: `plots/all16_masssq_heatmap_percent.png`
- mixing-angle-only heatmap: `plots/all16_mixing_heatmap_percent.png`
- corresponding PDF versions
- CSV heatmap tables under `csv/`

The font remains the Matplotlib default. Only the font sizes are user-controlled.

## Physics definition

These comparisons are unchanged:

- `Running = 100 * (1111 - 1001) / 1001`
- `Mixing off = 100 * (LLSS_OFF - 1001) / 1001`

For `--plot-only`, the script simply replots the existing `all16_three_way.csv`.

## PowerShell

From the project root:

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline

Expand-Archive `
  -LiteralPath "$HOME\Downloads\t3_original_llss_split_plots.zip" `
  -DestinationPath . -Force

# Replot from an existing A(-4) LLSS run without rerunning physics
.\.venv\Scripts\python.exe scripts\RunT3OriginalSimultaneousLLSSHeatmap.py `
  --plot-only output\studies\general_scalar_extended\reports\llss_original_all16_sourcefix_v2 `
  --output output\studies\general_scalar_extended\reports\llss_original_all16_splitplots `
  --font-size 16 `
  --cell-fontsize 12 `
  --top-tick-fontsize 18 `
  --bottom-tick-fontsize 15 `
  --y-tick-fontsize 15 `
  --figure-width 23
```

If you need to recompute the LLSS-off data first, use the same script without `--plot-only`, exactly as in the previous README.
