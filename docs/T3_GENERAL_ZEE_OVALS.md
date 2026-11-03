# Section 5.3 — General T3 Zee-style sensitivity figures (Windows-safe v2)

## Purpose and scope

This replaces the original ellipse-grid rendering in
`scripts/PlotT3GeneralZeeOvals.py` with a plotting style based on **Fig. 9** of
Schmidt & Vandeleur, *Running of neutrino mass parameters in the Zee model*.

- **All 16 T3 models** are plotted, not only B(-1) and C(-1).
- The saved **merged-fit** UV sensitivity data are read without changing fitted
  Yukawas or running RGEs/Matchete again.
- Three coloured gamma trajectories T1/T2/T3 and a black star at gamma=1 are
  drawn on each panel. Endpoints are marked.
- Filled **marginal** NuFIT 6.1 NO ±1σ/±3σ (blue) and measured JUNO 2026 ±1σ
  (green) bands replace the previous outline-only Gaussian ellipses.
- These bands are **one-dimensional marginal intervals**, with vertical and
  horizontal strips; their intersection is *not* a joint 2D confidence contour.
- An optional `--model-oval` flag superposes the descriptive empirical response
  oval. That oval is **not** a statistical confidence region.
- The matching input remains the existing conditional general UV/EFT1 RGE and
  one-scale C5 result, **not** a certified full sequential EFT Wilson match.

The script does not need the historical restricted scalar potential results.

## Why the old plot command failed on Windows

The user's traceback ended in `PIL.Image.save` when trying to open an existing
PNG destination under `Reports/output/thesis_general_5_3/zee_general_merged_ovals_v1`.
The log alone does **not** establish whether an existing file lock, access
control, or a Windows/Pillow pathname interaction caused `OSError [Errno 22]`.

This replacement:

1. Writes by default to a **new** directory, `zee_style_all16_v2`.
2. Saves each image to a `BytesIO` buffer first; Pillow never opens the Windows
   destination pathname directly.
3. Writes **individual models before** the large 4×4 grid.
4. Saves **PNG and vector PDF** for every figure.
5. Supports `--no-grid` if only the large combined figure is problematic.

## Install paths

The ZIP contains complete replacements, preserving the repository structure:

- `scripts/PlotT3GeneralZeeOvals.py`
- `tests/test_t3_general_zee_ovals.py`
- `docs/T3_GENERAL_ZEE_OVALS.md`

No GitHub changes. No physics code or saved scan modifications.

## PowerShell commands

Run from `C:\thesistemp\T3neutrinomassmodelpipeline`:

```powershell
Expand-Archive -LiteralPath "$HOME\Downloads\T3_General_Zee_Style_Windows_Fix_v2.zip" `
  -DestinationPath . -Force

python -X utf8 -m pytest -q tests\test_t3_general_zee_ovals.py

python -X utf8 -m scripts.PlotT3GeneralZeeOvals
```

Open the new folder:

```powershell
Start-Process "Reports\output\thesis_general_5_3\zee_style_all16_v2"
```

Output structure:

```
Reports/output/thesis_general_5_3/zee_style_all16_v2/
    model_grid_zee_solar.png
    model_grid_zee_solar.pdf
    models/A_m4_zee_solar.png
    models/A_m4_zee_solar.pdf
    ...
    models/E_p2_zee_solar.png
    models/E_p2_zee_solar.pdf
    local_sensitivity_slopes.csv
    solar_sensitivity_ranking.csv
    baseline_nufit_juno_comparison.csv
    provenance.json
```

The 4×4 grid and all individual figures use the saved data at:

`output/studies/general_scalar_extended/zee_style_general_sensitivity/zee_general_all16_merged_v1/`

If Windows still rejects a large grid PNG (or you want only the 16 figures):

```powershell
python -X utf8 -m scripts.PlotT3GeneralZeeOvals `
  --no-grid `
  --output-dir "Reports\output\thesis_general_5_3\zee_style_models_only_v2"
```

If an existing output PNG is open/locked, choose another new output directory:

```powershell
python -X utf8 -m scripts.PlotT3GeneralZeeOvals `
  --output-dir "Reports\output\thesis_general_5_3\zee_style_all16_v2_new"
```

To overlay the optional *descriptive* model covariance oval:

```powershell
python -X utf8 -m scripts.PlotT3GeneralZeeOvals `
  --model-oval `
  --output-dir "Reports\output\thesis_general_5_3\zee_style_with_model_ovals_v2"
```

## Scientific interpretation

The coloured T1/T2/T3 traces represent full numerical response to the chosen
UV compensated-coupling perturbation over gamma in [0.99, 1.01]. They are not
oscillation-fit credible regions. Near-perfect fits may be almost invisible
at the NuFIT/JUNO uncertainty scale, while B(-1)/C(-1) can have long solar
mass-splitting excursions. Never interpret the shaded marginal bands as a
combined NuFIT+JUNO likelihood.

Data references (in the source file): `configs/nufit_6_1_ic24_no.json` and
JUNO 207.2 live-day solar measurement, arXiv:2609.37895.
