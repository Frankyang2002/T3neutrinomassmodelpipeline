# T3 general-Lagrangian sensitivity: dual-view figures (v5)

## Scope

This is a plotting-only update for the completed 16-model, 528-point
Zee-style T1/T2/T3 general-T3 scan. It **does not** refit the Yukawas,
rerun the UV/EFT1 RGEs, integrate out particles, or modify any scan outputs.

Input (unchanged):

`output/studies/general_scalar_extended/zee_style_general_sensitivity/zee_general_all16_merged_v1/`

The fitted source uses the individually optimised merged 16-model candidate
Yukawa fits. The underlying calculations retain **conditional one-scale**
C5 matching, not certified complete sequential matching. The intended physical
EFT sequence is UV -> integrate out all fermions -> integrate out all scalars
-> SMEFT.

## Changes from v4

- **No arrows** on any T1/T2/T3 trajectory.
- **Larger plotted point markers** (4.9 pt for detail; 4.2 pt for overview).
- **No black centre stars**; the gamma=1 point is simply one of the saved
  trajectory points.
- Individual panel titles are just `A(-4)`, `B(-1)`, etc.; no
  `zoom to trajectories` suffix. The combined close-up title is now
  `Detail: compensated-coupling sensitivity for all 16 models`.
- **No bottom explanatory footnote** on either combined grid.
- **Readable close-up axis numbers:** for each model and observable, plot
  `delta O = O(gamma) - O(gamma=1)` and divide by a labelled power of ten.
  For example, the x label might read `delta(sin^2 theta_12) [10^-8]` and
  the y label `delta(Delta m^2_21) [10^-11 eV^2]`. Tick labels are therefore
  small ordinary numbers; Matplotlib's cramped `1e-7+3.088e-1` offset
  notation is never used on the close-up grid. Scaling varies *independently*
  by model and observable and is stated in every axis label.
- The **overview uses the original, absolute observable values**, so the
  experimental regions remain directly interpretable.

The coordinates and covariance ellipses are subject to the *same affine
coordinate changes*, so no slopes, relative orientation, gamma points or
response oval shapes are altered. The descriptive model ovals can be
extremely thin when all three sensitivity directions are nearly collinear.
No artificial transverse variance is added.

## Statistical interpretation

- NuFIT 6.1 (NO, IC24): independent-error 68.27% and 99.73% two-dimensional
  Gaussian *approximations* using marginal errors, not official contours.
- JUNO 2026 (207.2 live days): independent-error 68.27% Gaussian ellipse
  from the measured two solar marginal errors, not the older six-year forecast.
- Purple dashed model response ovals: 2 times the sample covariance standard
  deviation of the 30 off-centre saved gamma perturbations, centred on the
  gamma=1 point. **These are descriptive, not confidence regions.**

## Installation (PowerShell, local repository root)

The ZIP contains full replacements for:

- `scripts/PlotT3GeneralZeeOvals.py`
- `tests/test_t3_general_zee_ovals.py`
- `docs/T3_GENERAL_ZEE_OVALS.md`

```powershell
Expand-Archive -LiteralPath "$HOME\Downloads\T3_General_Zee_Ovals_Clean_DualView_v5.zip" `
    -DestinationPath . -Force
python -X utf8 -m pytest -q tests\test_t3_general_zee_ovals.py
python -u -X utf8 -m scripts.PlotT3GeneralZeeOvals
```

## Output

Default destination (new, distinct from v4):

`Reports/output/thesis_general_5_3/zee_ovals_all16_dualview_v5/`

Two combined 4x4 panels (both PNG and vector PDF):

- `model_grid_zee_zoom_out.png` — experimental ellipses + model trajectories
  using **absolute** solar observables.
- `model_grid_zee_zoom_in.png` — close-up of the same T1/T2/T3 trajectories
  and descriptive ovals, using **baseline-subtracted and power-of-ten** axes.

The same two views are generated individually for every model in `models/`,
with filenames such as `B_m1_zee_zoom_in.pdf`. The diagnostic CSV files
(`solar_sensitivity_ranking.csv`, `local_sensitivity_slopes.csv`, and
`baseline_nufit_juno_comparison.csv`) are still generated unchanged.

To open the combined PNGs in PowerShell:

```powershell
$figs = "Reports\output\thesis_general_5_3\zee_ovals_all16_dualview_v5"
Start-Process "$figs\model_grid_zee_zoom_out.png"
Start-Process "$figs\model_grid_zee_zoom_in.png"
```

You can select another folder using `--output-dir <directory>`. Optional
`--no-grid`, `--no-per-model`, `--no-model-oval`, and
`--juno-mode measured|future|both` are retained from v4. The default is
**measured** JUNO 2026, not a forecast. Neither input scans nor existing
figures are overwritten by default.
