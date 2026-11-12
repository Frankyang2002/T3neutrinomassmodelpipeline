# Thesis-consistent EFT1 LLSS heatmaps

Project file: `scripts/RunT3EFT1LLSSThesisHeatmap.py` (new/additive).
Test: `tests/test_t3_eft1_llss_thesis_heatmap.py` (new/additive).

## Physics definitions

Both columns are percentages relative to the SAME original thesis `1001` baseline:

- **Running:** `100*(O_FULL_1111 - O_FROZEN_1001)/O_FROZEN_1001`.
- **Mixing off:** `100*(O_LLSS_OFF - O_FROZEN_1001)/O_FROZEN_1001`.

`1001` freezes **both** renormalisable EFT1 running and dimension-five Wilson running, while retaining the UV and SMEFT running stages and both thresholds. `1111` runs everything. `LLSS_OFF` runs everything **except** the quartic `LLSS -> LLHH` Wilson RGE source term, subtracted at every Wilson RHS evaluation. Both scalars are integrated out simultaneously at MS.

### Numerical integrity

Original thesis FULL and FROZEN values are read directly from `output/parallel_t3/thesis_two_benchmarks_all16_stagefreeze_v2/running_comparisons.csv` (only `stage_frozen=EFT1`). Its stored percentages are verified against the raw observables. FULL values are checked against the original 32-case result in `output/parallel_t3/thesis_generalised_tables_v1/results.json`. The same A(-4)/D(-2) donor inputs and fixed Yukawas from `A4_Dm2_targeted_fit01` are verified against the original thesis table provenance. The only newly computed physics is `LLSS_OFF`.

This deliberate reuse means the **Running** column is exactly the thesis stage-freeze result, not a recomputation with a different frozen definition. The code refuses missing or mismatched input, does not substitute missing data and does not edit GitHub.

The known local source hash mismatch in `RGE/running/rgbeta/RunGeneralT3QuadraticRGBeta.wl` is permitted by default, while cached UV/EFT1 beta JSONs are still hash-verified using the existing contract, and all other required source mismatches remain blocked. Pass `--strict-batch1-sources` to disable this limited exception.

### Run (PowerShell, repository root)

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
.\.venv\Scripts\python.exe scripts\RunT3EFT1LLSSThesisHeatmap.py `
  --reference both `
  --output output\parallel_t3\eft1_llss_thesis_consistent_v1 `
  --font-size 16 `
  --cell-fontsize 12 `
  --top-tick-fontsize 19 `
  --bottom-tick-fontsize 16 `
  --y-tick-fontsize 17 `
  --axis-label-fontsize 18 `
  --colorbar-fontsize 15 `
  --figure-width 23
```

To smoke-test only A(-4), use `--reference A_-4 --models A_-4` and a new output directory. To redraw an already completed run with different fonts, rerun with `--plot-only` and the **same** output directory (only images and report are overwritten).

Matplotlib's native default text and math fonts are used. No figure title is present. The axes show five observable headings and paired `Running`/`Mixing off` labels.

The runner depends on the existing production `scripts/RunT3EFT1LLSSIsolation.py` for the tensor-source ablation, cached-beta verifier and plotting renderer. That file is **not** replaced by this package.

### Outputs

Per donor (`A_m4`, `D_m2`) under `--output`:

- `thesis_consistent_three_way.csv` (FULL, FROZEN, LLSS_OFF raw observables and percentages)
- `plots/all16_heatmap_percent.png` and `.pdf`
- `plots/all16_masssq_heatmap_percent.png` and `.pdf`
- `plots/all16_mixing_heatmap_percent.png` and `.pdf`
- `csv/all16_heatmap_table.csv`
- `reports/thesis_consistent_LLSS_report.md`
- `reports/thesis_consistent_LLSS_table.tex`
- `provenance.json`, per-model diagnostics

**Scientific limitation:** Candidate seven-slot Wilson ADM, scalar basis map, spinor-normalisation conversion and full independent finite F-hard matching are not certified. The figures are exploratory and should be qualified as such in the thesis.
