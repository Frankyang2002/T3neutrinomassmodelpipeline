# General-scalar 16-model interactive display

## Purpose and scope

This batch assembles **three offline HTML reports**, a home page, and two CSV tables from the completed candidate `continuous_v1` run. **No physics is rerun**: the generator reads only saved numerical JSON/CSV results. It never modifies the numerical inputs or repository on GitHub.

The output is a general-scalar-specific display and does **not** reuse the earlier reduced-potential `output/full/comparison` dashboard as if it were a new general-scalar result.

## Added project files

1. `scripts/BuildGeneralScalarInteractiveDisplay.py`: complete reporting CLI (Python standard library only).
2. `tests/test_general_scalar_interactive_display.py`: synthetic-fixture regressions.
3. `docs/T3_GENERAL_SCALAR_INTERACTIVE_DISPLAY.md`: this guide.

All are **new files**. No existing source is replaced.

## Install and run (PowerShell)

Run from `C:\thesistemp\T3neutrinomassmodelpipeline` after extracting the provided ZIP directly into the repository root:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
python -X utf8 -m pytest -q tests/test_general_scalar_interactive_display.py
python -m scripts.BuildGeneralScalarInteractiveDisplay --run-id continuous_v1
Start-Process .\Reports\output\display_general_scalar\continuous_v1\index.html
```

Open individual reports:

```powershell
Start-Process .\Reports\output\display_general_scalar\continuous_v1\interactive_rge.html
Start-Process .\Reports\output\display_general_scalar\continuous_v1\sensitivity_report.html
Start-Process .\Reports\output\display_general_scalar\continuous_v1\all_results.html
```

Input root (default):

`output/studies/general_scalar_extended/continuous_full_scalar/continuous_v1/`

Output root (default):

`Reports/output/display_general_scalar/continuous_v1/`

The generator requires all sixteen successful `MODEL/continuous_trajectory.json` records and the saved `neutrino_v1/summary.csv`. It reads additional matching summaries if present; it explicitly reports missing optional summaries. Missing/invalid required inputs stop **before** any report is written.

## Report contents

### `interactive_rge.html`

* All sixteen **actual saved** `uv_trajectory` and `eft1_trajectory` sample sequences, with independent UV/EFT1 selection.
* Running real scalar mass-squared parameters `mS1Sq_gev2` and `mS2Sq_gev2`.
* New general-scalar quartics in each saved sample's `extra_quartics`, with real, imaginary, or absolute-value display.
* Holomorphic quadratic `mu2_gev2`, including the zero/no-operator cases.
* Logarithmic energy-scale axis and fermion threshold marker; model and quantity selectors; numeric hover values.

**Not provided by these files:** a time-resolved physical `C5(mu)`, neutrino mixing-angle trajectory, or full SM+Weinberg running. The continuous scalar trajectory JSON only preserves the scalar subset above. These curves are **not** fabricated from endpoint results.

### `sensitivity_report.html`

Two explicitly distinct concepts are supported:

1. **Matching-approximation sensitivity**, immediately available from saved data: for each model and saved observable `X`, display

   `100 * [X(scheme) - X(neutrino_v1)] / abs[X(neutrino_v1])`.

   This describes a *change of theoretical approximation*, **not** the derivative of a fitted prediction with respect to a UV parameter. There is no statistical interpretation.

2. **Actual T1/T2/T3 gamma sensitivity**, **only** if you separately supply a saved *general-scalar* response CSV generated using the matching convention of the corresponding benchmark. Sample paths and optional 1-standard-deviation sample-covariance ellipses are then drawn; these ellipses are descriptive, **not experimental or confidence contours**. Without such a file, the panel states that sensitivity ovals are unavailable. It never copies old reduced-potential data or synthesizes gamma points.

Optional (already computed scan data only):

```powershell
python -m scripts.BuildGeneralScalarInteractiveDisplay --run-id continuous_v1 `
  --sensitivity-csv .\output\studies\general_scalar_extended\your_general_scalar_gamma_scan.csv
```

The scan CSV schema is **one row per saved gamma sample**, with mandatory `model`, `direction` (`T1`, `T2`, `T3`), and `gamma` (positive). At least 3 distinct gamma values including `gamma=1` are required **per model/direction**. Supply at least 2 numeric observables from:

- `mass_hi_ev`, `mass_mid_ev`, `mass_lo_ev` (eV, descending),
- `delta_m21_sq_ev2`, `delta_m3l_sq_ev2`,
- `sin2_theta12`, `sin2_theta13`, `sin2_theta23`.

This display script **does not generate the scan CSV**. Generating one correctly requires rerunning the general-scalar physics calculation under the chosen UV perturbations. It should not be inferred from central values.

### `all_results.html`, `all_results.csv`, `matching_variants.csv`

* One default candidate row for each of 16 models, with model, source status, three masses, and only the measured/saved observables available for that **same** source.
* A separate, sortable/filterable list containing **every saved distinct matching scheme**.
* Additional algebraic squared-mass gaps, **explicitly distinct** from the source's oscillation-ordering observables: `(middle² - smallest²)` and `(largest² - smallest²)`. These are always derived from the three saved nonnegative masses, without assuming their physical ordering.
* No mixing angles are copied into a B/C selected-source or coupled-Wilson row, because those corresponding CSV files save only masses.

The **default model-specific selection** is a presentation convenience, **not** a physically certified best prediction:

| Cases | Preferred available candidate |
|---|---|
| B(-1), C(-1) | Coupled `L/eC` Wilson candidate |
| A(+0), A(-2), D(+0), D(-2), E(+0), E(-2) | Mass-dependent self-bilinear one-scale candidate |
| B(+1), C(+1), D(+2), E(+2) | Direct sequential Wilson candidate |
| A(-4), B(-3), C(-3), E(-4) | Legacy benchmark control (no saved hierarchical hard boundary) |

If an optional source is absent, selection falls back in a documented order without inventing values. The rows across **different** matching schemes are **not** combined or interpreted as a unified finite-order EFT result. The generator specifically excludes the ten unrelated pass-through rows of `self_anomalous_massdep_v1/combined_16_summary.csv` from the six self-bilinear corrections.

## Expected outputs

```
Reports/output/display_general_scalar/continuous_v1/
    index.html
    interactive_rge.html
    sensitivity_report.html
    all_results.html
    all_results.csv
    matching_variants.csv
    provenance.json
    parameter_sensitivity_points.csv   # only when a real CSV was provided
```

HTML is self-contained (SVG + vanilla JavaScript + inline CSS); there are **no internet assets, CDN scripts, Matplotlib, or Plotly dependencies**.

The source manifest contains SHA-256 fingerprints of each input read, record counts, selected scheme counts, and limitations. Re-running the generator replaces only its *dedicated display folder files*; numerical source files remain untouched.

## Validation and limitations

The included tests use synthetic fixtures to check the actual saved schema, 16-model coverage, source selection, missing data handling, no fabricated sensitivity points, and withheld B/C mixing angles. They do **not** independently validate numerical RGEs or model-specific matching physics.

Known physics limitations: candidate scalar/RGBeta conventions, finite scalar threshold matching, self-bilinear mass basis, and missing full general-scalar sequential matching remain uncertified. A display does not remove these limitations.
