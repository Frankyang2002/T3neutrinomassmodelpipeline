# EFT1 LLSS isolation: 16-model heatmaps and donor reports

Project files in this ZIP (all full files; extract over the older ZIP):

- `scripts/RunT3EFT1LLSSIsolation.py`
- `tests/test_t3_eft1_llss_isolation.py`
- `README_EFT1_LLSS_ISOLATION.md`

## Physics definitions (unchanged)

Each model is evaluated in three numerical conditions using the SAME fixed
UV Yukawas, UV RG flow, all-fermion threshold, EFT1 renormalisable background,
joint two-scalar threshold, and SMEFT evolution.

- `W_FROZEN`: dimension-five Wilson beta set to zero (all EFT1 Wilson running off).
- `LLSS_OFF`: only the quartic `LLSS -> LLHH` source in the full tensor Wilson
  beta is removed *during* the numerical integration.
- `W_FULL`: all candidate Wilson beta terms retained.

The heatmap displays **two columns per observable**, both relative to the same
frozen-Wilson denominator:

```
Running   = 100 * (O_W_FULL   - O_W_FROZEN) / O_W_FROZEN
Mixing off = 100 * (O_LLSS_OFF - O_W_FROZEN) / O_W_FROZEN
```

Five observable groups: `dm21_sq_eV2`, `dm31_sq_eV2`, `sin2_theta12`,
`sin2_theta13`, `sin2_theta23`. The row inventory has 16 models; blocked
models produce blank cells, not invented values.

**Important distinction:** This is NOT the historical `EFT1 off` stage freeze,
which also disabled renormalisable scalar running. In the present diagnostic,
renormalisable scalar RGEs remain ON in all three cases. The purpose is to
isolate the Wilson LLSS quartic source.

## Plot styling

- No figure title. The thesis caption should supply the description.
- Top x-axis: five mathematical observable names.
- Bottom x-axis, repeated under each observable: **Running**, **Mixing off**.
- Matplotlib default font and math fonts are used, without any font-family
  override. Font sizes remain configurable for readability.
- The colourbar has a reserved column so it cannot overlap data.
- Both PNG and vector PDF figures are written (no LaTeX installation required).

Options include:

| Option | Default | Purpose |
| --- | --- | --- |
| `--axis-label-fontsize` | `12` | Y-axis label |
| `--bottom-tick-fontsize` | `10` | `Running`, `Mixing off` labels |
| `--top-tick-fontsize` | `13` | Five observable headings |
| `--y-tick-fontsize` | `11` | Model names |
| `--cell-fontsize` | `8` | Percentages inside cells |
| `--colorbar-fontsize` | `10` | Colourbar ticks and label |
| `--figure-width` | `16` | Figure width (inches) |
| `--row-height` | `0.48` | Height per model (inches), plus margin |
| `--dpi` | `300` | PNG export quality |

## Recommended PowerShell workflow

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline

Expand-Archive -LiteralPath "$HOME\Downloads\t3_eft1_llss_matplotlib_default.zip" `
  -DestinationPath . -Force

.\.venv\Scripts\python.exe -m unittest discover `
  -s tests -p test_t3_eft1_llss_isolation.py -v
```

### A(-4) donor: replot existing successful 16-model results (NO physics rerun)

The input directory is the one containing the successfully completed
`all16_three_way.csv` from the uploaded heatmap. The example below uses
the output path suggested in the previous run; change it if yours differs.

```powershell
.\.venv\Scripts\python.exe scripts\RunT3EFT1LLSSIsolation.py `
  --plot-only output\parallel_t3\eft1_llss_cached_all16 `
  --output output\parallel_t3\eft1_llss_A4_thesis_report `
  --bottom-tick-fontsize 10 --top-tick-fontsize 14 `
  --y-tick-fontsize 11 --cell-fontsize 9 `
  --colorbar-fontsize 11 --figure-width 16 --row-height 0.5
```

### D(-2) donor: run all 16 models with the same three Wilson switches

This is a **new donor calculation**, not merely relabelling the A(-4) data.
All original scenario physics is retained except the chosen Wilson beta
switch. The cached Batch 1 exception is retained from the previous runner
and still verifies the saved UV/EFT1 beta export hashes.

```powershell
.\.venv\Scripts\python.exe scripts\RunT3EFT1LLSSIsolation.py `
  --reference D_-2 `
  --reuse-cached-batch1 `
  --output output\parallel_t3\eft1_llss_Dm2_thesis_report `
  --bottom-tick-fontsize 10 --top-tick-fontsize 14 `
  --y-tick-fontsize 11 --cell-fontsize 9 `
  --colorbar-fontsize 11 --figure-width 16 --row-height 0.5
```

Use new output directories for these commands. The script deliberately
refuses to overwrite existing results. A failed/blocked model is retained as
a blank row and is identified in the report.

If D(-2) donor numerical results already exist from this runner, use
`--plot-only <existing_Dm2_results>` with a new `--output` instead of
recomputing.

## Outputs under each donor's chosen `--output` directory

Figures:

- `plots/all16_heatmap_percent.png` and `.pdf` (16 x 10)
- `plots/all16_masssq_heatmap_percent.png` and `.pdf` (16 x 4)
- `plots/all16_mixing_heatmap_percent.png` and `.pdf` (16 x 6)
- `plots/selected6_*.png` and `.pdf`

Numerical tables:

- `all16_three_way.csv` — three raw observable sets and C5 diagnostics
- `selected6_three_way.csv`
- `csv/all16_heatmap_table.csv` — 16 x 10 signed percentage matrix
- `csv/all16_heatmap_long.csv` — long-format representation
- split CSVs for mass-squared differences and mixing angles

Donor-specific **reports** generated from the computed data:

- A(-4): `reports/eft1_llss_A_m4_report.md`, `.tex`, `_summary.json`
- D(-2): `reports/eft1_llss_D_m2_report.md`, `.tex`, `_summary.json`

Each report contains definitions, full 16-model comparison tables for all
five observables, median absolute percentages, completion status, and
methodological qualifications. The `.tex` fragment can be inserted into the
thesis; it does not overwrite the thesis source.

## Numerical scope

Results remain exploratory/conditional: source-to-Weyl normalisation,
independent hard matching, complete scalar maps, and the full candidate ADM
are not fully certified. Large mass-squared percent changes with very small
angle changes can arise if the neutrino mass matrix is approximately rescaled.
The B(-1)/C(-1) self-/mixed-channel exceptions are not suppressed in the
reported table. Never interpret results from the A(-4) Yukawa donor as
predictions for a D(-2) donor without a new numerical run.

## Default Matplotlib font (no override)

The plotter uses the default Matplotlib font settings. There is no
`--font-family` flag, and no explicit text or math font selection in the code.
Existing saved A(-4) or D(-2) results can be replotted without rerunning RGEs.

```powershell
.\.venv\Scripts\python.exe scripts\RunT3EFT1LLSSIsolation.py `
  --plot-only output\parallel_t3\eft1_llss_cached_all16 `
  --output output\parallel_t3\eft1_llss_A4_default_font
```

The `--output` path must be new. Use a completed D(-2) results directory
as the `--plot-only` input for its corresponding report.

## Combined A(-4) + D(-2) thesis reports (Matplotlib default font, no titles)

The additional file `scripts/RunT3EFT1LLSSBothDonors.py` creates two donor
subdirectories in a **new** output directory. It uses the existing three-way
runner and does not change physics. It always removes figure titles and uses
the default Matplotlib font. Adjust `--font-size` (all text) or any individual
size flag.

If A(-4) was already computed, but D(-2) still needs its own 16-model run:

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline

.\.venv\Scripts\python.exe scripts\RunT3EFT1LLSSBothDonors.py `
  --a4-source output\parallel_t3\eft1_llss_cached_all16 `
  --reuse-cached-batch1 `
  --output output\parallel_t3\eft1_llss_both_donors_v1 `
  --font-size 10 `
  --top-tick-fontsize 13 `
  --cell-fontsize 8
```

If both sets of numerical data have already been computed, replot only:

```powershell
.\.venv\Scripts\python.exe scripts\RunT3EFT1LLSSBothDonors.py `
  --a4-source output\parallel_t3\eft1_llss_cached_all16 `
  --d2-source output\parallel_t3\eft1_llss_Dm2_thesis_report `
  --output output\parallel_t3\eft1_llss_both_donors_v2 `
  --font-size 10
```

If `--a4-source` or `--d2-source` is omitted, that donor is computed afresh;
`--reuse-cached-batch1` is applicable to freshly calculated runs only.
Sources must contain 16 distinct models and match the chosen donor reference.
Each donor retains the same five-observable, two-column-per-observable heatmap:
**Running versus frozen** and **Mixing off versus frozen**.

Main files (relative to chosen `--output` directory):

- `A_m4/plots/all16_heatmap_percent.png` and `.pdf`
- `A_m4/reports/eft1_llss_A_m4_report.md` and `.tex`
- `D_m2/plots/all16_heatmap_percent.png` and `.pdf`
- `D_m2/reports/eft1_llss_D_m2_report.md` and `.tex`
- `both_donors_summary.json` (which datasets were reused/recomputed)

Font controls:

- `--font-size`: all figure text sizes together, unless individually overridden
- `--axis-label-fontsize`, `--bottom-tick-fontsize`, `--top-tick-fontsize`
- `--y-tick-fontsize`, `--cell-fontsize`, `--colorbar-fontsize`
- `--figure-width`, `--row-height`, `--dpi`

No `--font-family` or title flag is needed.
