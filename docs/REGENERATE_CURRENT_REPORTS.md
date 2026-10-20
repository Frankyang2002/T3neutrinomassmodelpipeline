# Current T3 reports: rebuild and retirement guide

Checked against public GitHub `main` commit `b75a1408307fe738b0dca25412011788b836eec0` (2026-10-20). This ZIP provides **local replacement files**, not changes to GitHub. Existing fitted and A-m2 source data is preserved.

## 1. Install

Place the accompanying ZIP at the project root and execute in PowerShell:

```powershell
Expand-Archive -LiteralPath .\T3_FERMION_THRESHOLD_REPORT_RENAME.zip -DestinationPath . -Force
```

Project paths installed:

- `Numerical/plotting/FermionThresholdVariationFigures.py` — **new** read-only plot renderer
- `tools/RebuildCurrentReports.py` — complete replacement for the prior provided local helper
- `tools/RebuildReports.py` — complete replacement for the GitHub and previously supplied helper
- `tools/RebuildSavedStudyPlots.py` — complete replacement for the GitHub helper
- `tools/LabelSensitivityDirections.py` — complete replacement from the previous sensitivity-label package (bundled for standalone installation)
- `scripts/RetireOldStudyDisplays.ps1` — **new**, safe preview-then-apply removal tool
- `docs/REGENERATE_CURRENT_REPORTS.md` — this guide

## 2. Regenerate the renamed plots and other displays

```powershell
# New presentation: read saved scan.json, no numerical reruns
python -m tools.RebuildReports --study fermion-threshold --study weinberg

# (Optional) Rebuild the separately fitted heavy-fermion mass scan too
python -m tools.RebuildReports --study fitted-mass

# Rebuild all supported HTML dashboards and static figures
python -m tools.RebuildCurrentReports
```

You may alternatively use:

```powershell
python -m tools.RebuildCurrentReports --only displays
python -m tools.RebuildCurrentReports --only plots
python -m tools.RebuildSavedStudyPlots --study fermion-threshold
python -m tools.RebuildSavedStudyPlots --study weinberg
```

**Note:** `--only displays` does not rebuild static PNGs; run the plot command as well after replacing the files.

## 3. Output locations

| New presentation output | Source saved numerical data (unchanged) |
|---|---|
| `Reports/output/display/fermion_threshold_variation/smallY_smallL/fermion_threshold_running_effect_*.png` | `Reports/output/full/staged_rg_comparison/smallY_smallL/scan.json` |
| `Reports/output/display/weinberg_running/smallY_smallL/weinberg_running_*.png` | `Reports/output/full/weinberg_running_comparison/smallY_smallL/scan.json` |
| `Reports/output/full/fermion_scale/smallY_smallL/fermion_scale_*.png` | `Reports/output/full/fermion_scale/smallY_smallL/scan.json` |

The historical *numerical source* folder name `staged_rg_comparison` is retained intentionally; moving it would break the older numerical runner. The **user-facing report title, PNG filename, plot folder, and supported CLI name** are all renamed.

Every model legend now uses `T3A-m2`, `T3B-p0`, etc. The underlying canonical project short-name format is `T3-A-m2` (with a hyphen after T3); removing that extra hyphen is **presentation only** to match the requested style. All physical model identifiers remain unchanged.

### Correct physics labels

The three panels correspond to:

- UV: `000 → 100`, percentage `(O100/O000 - 1) × 100`;
- EFT1: `100 → 110`, percentage `(O110/O100 - 1) × 100`;
- Weinberg: `110 → 111`, percentage `(O111/O110 - 1) × 100`.

Each pair compares the **same final low-energy observable** with the stated segment switched on versus frozen. In all panels, matching is retained. In the scan, both the heavy fermion mass matrix `M_F` and the threshold `mu_F` vary, while `mu_S` remains fixed; therefore the horizontal axis `mu_F/mu_S` is a **threshold variation**, not an RGE time coordinate. The percentages are not a rate of change or a per-energy-step derivative.

The separate Weinberg figure compares evolved versus frozen `C5` in the final EFT. The fitted mass-scan panel instead compares to the same model's baseline at mass factor one; **it is not** an RGE switch-off test.

## 4. Retire old generated outputs after successful rebuild

```powershell
# Preview what will be removed
.\scripts\RetireOldStudyDisplays.ps1

# Apply after both sets of five replacement figures have been created
.\scripts\RetireOldStudyDisplays.ps1 -Apply
```

This retires:

- `Reports/output/display/study_comparisons/` (Study A/B interactive HTMLs);
- `Reports/output/display/rebuilt_studies/` (misleading old folder name, including old staged/Weinberg figures);
- `Reports/output/display_stage5/rebuilt_studies/` if present;
- `Reports/output/full/staged_rg_comparison/smallY_smallL/staged_rg_*.png` (old titled/labelled images **only**);
- `Reports/output/full/weinberg_running_comparison/smallY_smallL/weinberg_running_*.png` (old titled/labelled images **only**);
- retired `Reports/output/display/interactive_correlations.html` if still present.

**It preserves** `Reports/output/full/staged_rg_comparison/smallY_smallL/scan.json`, `scan.csv`, and all underlying fitted benchmark / symbolic / numerical artifacts.

The retirement script defaults to preview mode, refuses `-Apply` until five files exist in **each** renamed plot folder, and does not delete arbitrary report directories.

## 5. What is retired, and what is not

The supported rebuild wrappers **no longer call** `Numerical.plotting.BuildStudyABInteractive` and no longer produce the Study A/B interactive HTML files. The original standalone historical generator remains in the repository because removing research source code is not needed to retire generated reports. Do not run it if you want the HTML retired.

Similarly, the old physics study runners `RunFittedStagedRGComparison.py` and `RunFittedWeinbergRunningComparison.py` are **not modified** by this plotting-only package. They still write source studies in the historical `Reports/output/full/` directories and may emit old-style PNGs if explicitly rerun with `--run`. To refresh the renamed figure presentation, use the supported `tools.RebuildReports` or `tools.RebuildCurrentReports` commands after any numerical rerun. No numerical study is re-executed by the new plot utilities.

The separate `fitted_model_envelopes_all.html` and `common_input_model_envelopes_all.html` retain their T1/T2/T3 parameter-deformation labels. They are unrelated to this threshold-variation renaming.
