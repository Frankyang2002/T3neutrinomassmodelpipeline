# General-scalar `continuous_v1`: original Display Lambda HTML presentation

## Scope

This replaces the **standalone SVG-based** `scripts/BuildGeneralScalarInteractiveDisplay.py` experience with the **existing Display Lambda HTML/canvas presentation**. The implementation reads the original HTML templates directly from:

- `Numerical/plotting/InteractiveModelComparison.py`: RGE plot layout, zoom, model selectors, log-scale axis and threshold markers.
- `Numerical/plotting/BuildFittedModelEnvelopes.py`: observable-pair plotting, T1/T2/T3 directions, gamma filters and sample-covariance ovals.

The original files are **not edited**. `Numerical/plotting/BuildLambda5e9Envelopes.py` is also unchanged; its numerical scan driver is **not called**.

The only new reporting Python entry point is `Numerical/plotting/BuildGeneralScalarLambdaReports.py`. It builds four linked offline HTML files and CSV/provenance outputs, using the saved `continuous_v1` results only. It never reruns an RGE, a fit, a threshold matching calculation or a sensitivity scan.

## Install / remove previous generator (PowerShell)

Extract `T3_general_scalar_lambda_reports.zip` directly into the local repository root, keeping the directory structure. Then execute:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
powershell -ExecutionPolicy Bypass -File .\scripts\RemoveOldGeneralScalarDisplay.ps1
python -X utf8 -m pytest -q tests\test_general_scalar_lambda_reports.py
python -X utf8 -m Numerical.plotting.BuildGeneralScalarLambdaReports --run-id continuous_v1
Start-Process .\Reports\output\display_general_scalar\continuous_v1\index.html
```

No `git push`, `git add`, `git commit` or other remote write is involved. The removal script deletes **only** these two superseded local files if present:

- `scripts/BuildGeneralScalarInteractiveDisplay.py`
- `tests/test_general_scalar_interactive_display.py`

The documentation file `docs/T3_GENERAL_SCALAR_INTERACTIVE_DISPLAY.md` is replaced by this guide.

## Outputs

All generated output stays under:

`Reports/output/display_general_scalar/continuous_v1/`

| File | Content |
|---|---|
| `index.html` | Home page with links to the three reports, CSVs and provenance |
| `interactive_rge.html` | Original Lambda-style interactive RGE plot for all 16 saved scalar trajectories |
| `sensitivity_report.html` | Original Lambda-style T1/T2/T3 gamma/oval page; reports unavailable scans honestly |
| `all_results.html` | Sortable/filterable 16-model preferred table and separate matching-scheme view |
| `all_results.csv` | Exactly 16 selected model rows, including their matching source |
| `matching_variants.csv` | All available independent saved scheme records |
| `provenance.json` | SHA-256 of consumed input files and reused templates; counts and limitations |
| `parameter_sensitivity_points.csv` | Only when a real optional general-scalar scan CSV was supplied |

## Data actually plotted

The required input directory is:

`output/studies/general_scalar_extended/continuous_full_scalar/continuous_v1/`

The 16 `MODEL/continuous_trajectory.json` records contain 16 sampled points in each of `uv_trajectory` and `eft1_trajectory` for the inspected dataset. `interactive_rge.html` plots the saved:

- `mS1Sq_gev2`, `mS2Sq_gev2` in GeV²;
- `mu2_gev2` real, imaginary and absolute components where present;
- every actual `extra_quartics` component, with real, imaginary and magnitude views;
- UV and EFT1 scale intervals with fermion threshold and scalar endpoint markers.

The scalar trajectory files do **not** provide a physical sampled `C5(μ)`, PMNS/mixing angle evolution, or the full post-scalar SM+Weinberg running. The replacement does **not** infer those curves from endpoint neutrino masses. As in the original Lambda plot, browser selectors and y-axis zoom do not recalculate the physics.

## Matching scheme table

The reporting script reads distinct candidates from:

1. `neutrino_v1/summary.csv` (all 16 baseline benchmark results).
2. `self_anomalous_massdep_v1/combined_16_summary.csv` (**six** conditional mass-dependent self-bilinear corrections; excludes ten pass-through rows).
3. `direct_wilson_sequential_v1/summary.csv` (saved direct Wilson results).
4. `bc_extra_wilson_sequential_v1/summary.csv` (B(-1)/C(-1) selected extra-quartic results).
5. `bc_coupled_wilson_v1/summary.csv` (B(-1)/C(-1) coupled L/eC candidate).

The preferred table uses the available B/C coupled source for B(-1)/C(-1), the mass-dependent self-bilinear candidate for the six self cases, and direct Wilson transport for the remaining models where available, with baseline fallback. This is only a **presentation selection**, not a combined or certified best-order EFT calculation.

Neutrino masses are reported largest/middle/smallest rather than assigning neutrino mass-eigenstate labels without ordering information. The displayed algebraic squared-mass gaps are explicitly distinguished from saved oscillation `Δm²`. If a matching source saved no mixing angles, those columns remain blank; no baseline values are silently copied across matching approximations.

## Genuine T1/T2/T3 sensitivity

**No real `continuous_v1` general-scalar gamma scan was verified in the inspected repository.** Consequently, the sensitivity page intentionally contains **no fabricated ovals** by default. The original fitted-Lambda envelopes from historical reduced-potential runs are never loaded.

If a genuine general-scalar T1/T2/T3 parameter perturbation CSV is later **already saved** under `continuous_v1`, enable it by adding:

```powershell
python -X utf8 -m Numerical.plotting.BuildGeneralScalarLambdaReports `
    --run-id continuous_v1 `
    --sensitivity-csv .\output\studies\general_scalar_extended\continuous_full_scalar\continuous_v1\general_scalar_gamma_points.csv
```

The optional CSV must be **inside the selected `continuous_v1` run directory**, and must use one row per genuinely evaluated model/direction/gamma sample:

`model,direction,gamma,mass_hi_ev,mass_mid_ev,mass_lo_ev,delta_m21_sq_ev2,delta_m3l_sq_ev2,sin2_theta12,sin2_theta13,sin2_theta23`

Requirements: directions `T1`, `T2`, `T3`; at least three distinct positive gamma values per direction, including `gamma=1`; at least two fully sampled common observables; consistent `gamma=1` centres across directions. The optional columns may be omitted except for at least two real numeric observables per point. Ellipses are **descriptive 2× empirical sample-covariance envelopes**, not experimental confidence contours. No NuFIT/JUNO contours are shown unless separately verified in the original template; their controls are disabled here.

Generating a genuinely new scan would require **new physics calculations**, so it is deliberately outside the scope of this reporting-only replacement.

## Tests and repeatability

`tests/test_general_scalar_lambda_reports.py` uses synthetic fixtures to verify the schema contract, 16-model coverage, no cross-source mixing contamination, exclusion of historical sensitivity input, and failure before writing outputs on malformed scalar trajectories. It does **not** establish numerical validity of the original physics integration.

`provenance.json` records each saved input and each original Lambda HTML generator SHA-256, making the exported reports reproducible against the exact local files used. Generated HTML and tables are self-contained and need no network connection to open.
