# T3 Section 5.3: merged-fit Zee sensitivity against NuFIT and JUNO

**New file, no physics-code modification:** `scripts/PlotT3GeneralZeeOvals.py`.

Read-only postprocessor for the completed 16-model campaign:

`output/studies/general_scalar_extended/zee_style_general_sensitivity/zee_general_all16_merged_v1/`

It reads `points.csv` and `summary.json`; it does **not** refit, run UV RGEs,
change the matched coefficients, or modify the existing scans.

## Experimental inputs

- **NuFIT 6.1**, November 2025, IC24 + SK atmospheric, normal ordering: repository
  `configs/nufit_6_1_ic24_no.json`. The five uncertainties are Gaussianised
  symmetric *marginal* 1σ errors; the correlation matrix is not supplied.
- **JUNO**, 207.2-live-day measurement (JUNO Collaboration, arXiv:2609.37895,
  29 September 2026):
  `sin²θ12 = 0.3036 ± 0.0064`,
  `Δm²21 = (7.388 ± 0.078) × 10⁻⁵ eV²`.
  These are **measured**, not the older JUNO six-year projection.

For display, both experimental 2D regions assume *zero correlation* and use
`r² = 2.30` (68.27% joint) and `r² = 11.83` (99.73% joint).
These are visual approximations, **not official experimental joint likelihood contours**.
JUNO data are drawn only for the solar (`sin2_theta12`, `delta_m21_sq_ev2`) pair.

The **model-sensitivity oval** is computed from all 30 off-centre scan points
for each model:

`Σ_model = Σ_p (O_p - O_gamma=1)(O_p - O_gamma=1)^T/(30-1)`.

The drawn semi-axes are twice the square roots of its two eigenvalues. It is
an empirical descriptive variation envelope, **not a model confidence region**.
The central sensitivity coefficient is the finite difference

`dO/dlnγ ≈ [O(γ+) - O(γ-)]/[ln(γ+) - ln(γ-)]`

with nearest saved gamma points on either side of 1. CSV summaries report
an estimated one-percent UV variation `0.01*dO/dlnγ` in NuFIT/JUNO 1σ units.
The JUNO solar χ² in the output is only the sum of two *diagonal Gaussian*
residual squares, **not the official JUNO likelihood** and not the same as
saved five-observable NuFIT fit χ².

## Installation and PowerShell commands

Extract the provided ZIP **into your local repository root**:

```powershell
Expand-Archive -LiteralPath "$HOME\Downloads\T3_General_Zee_Ovals_2026.zip" `
    -DestinationPath . -Force
python -X utf8 -m pytest -q tests\test_t3_general_zee_ovals.py

python -X utf8 -m scripts.PlotT3GeneralZeeOvals
```

Outputs:

`Reports/output/thesis_general_5_3/zee_general_merged_ovals_v1/`

- `model_grid_sin2_theta12_vs_delta_m21_sq_ev2.png`: 16-panel model comparison.
- `models/*.png`: individual model plots with NuFIT/JUNO overlays.
- `local_sensitivity_slopes.csv`: five-observable central finite-difference slopes.
- `solar_sensitivity_ranking.csv`: solar sensitivity measured in JUNO/NuFIT σ units.
- `baseline_nufit_juno_comparison.csv`: gamma=1 fit and independent solar Gaussian χ².
- `provenance.json`: sources, definitions, restrictions, input hashes.

Other observable pairs:

```powershell
python -X utf8 -m scripts.PlotT3GeneralZeeOvals `
    --x sin2_theta13 --y sin2_theta23 `
    --output-dir Reports\output\thesis_general_5_3\atmospheric_ovals_v1
```

## EFT scope

The saved general-Lagrangian sensitivity scans use the existing **conditional**
renormalisable UV/effective-scalar RGE and one-scale C5 framework.
They are **not** a certified physical sequential Wilson calculation
`UV → integrate out ALL fermions → integrate out ALL scalars → SMEFT`.
Do not describe these candidate ovals as a full-general EFT prediction.
