# Original EFT1 running-column reconciliation

Project replacement path: `scripts/RunT3OriginalSimultaneousLLSSHeatmap.py`.

This full replacement retains the original LLSS-off calculation and all three plots. New `--stagefreeze-csv` option (with `--plot-only`) reads **EFT1** rows from the original three-stage `frozen_stage_observable_changes.csv` via the repository's `PlotT3OscillationHeatmaps.read_rows`. It uses these numbers for the displayed Running columns. Mixing-off numbers are still from the supplied `all16_three_way.csv`.

**Scientific caveat:** When the Running percentages differ, the LLSS-off CSV and original stage-freeze CSV have not been shown to share the same 1001 denominator. The output is only a visual juxtaposition; do not interpret it as a consistent numerical ablation until the source runs are reconciled. `running_source_audit.csv` and `.json` report every difference and baseline-consistency warning. No RGE rerun.

First locate original CSV:

```powershell
Get-ChildItem output\studies\general_scalar_extended -Recurse -Filter frozen_stage_observable_changes.csv | Select-Object FullName
```

Then run (adjust `--stagefreeze-csv` to the CSV actually used for your original figure):

```powershell
.\.venv\Scripts\python.exe scripts\RunT3OriginalSimultaneousLLSSHeatmap.py `
  --plot-only output\studies\general_scalar_extended\reports\llss_original_all16_sourcefix_v4 `
  --stagefreeze-csv output\studies\general_scalar_extended\reports\simultaneous_16_stagefreeze_original_heatmaps_v1\frozen_stage_observable_changes.csv `
  --output output\studies\general_scalar_extended\reports\llss_original_running_reconciled_v1 `
  --font-size 16 --cell-fontsize 18 --top-tick-fontsize 18 `
  --bottom-tick-fontsize 15 --y-tick-fontsize 15 --figure-width 12
```

Four-column figure: `plots/all16_masssq_heatmap_percent.png` (also PDF). Other figures remain available. Inspect `running_source_audit.csv` before using results in thesis.
