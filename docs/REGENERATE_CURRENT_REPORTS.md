# T3 saved-data report rebuild — 20 October 2026

Source inspected: GitHub `main` at `b3ae6d1ca812554a16a1fd450f01a0c895e9b90e`. This package is for **local installation only**; it does not change GitHub. Inspect local edits with `git status`/`git diff` before overwriting local files.

## Installation from PowerShell

From `C:\thesistemp\T3neutrinomassmodelpipeline`:

```powershell
# Refuse to overwrite a newer checkout or changed local versions.
$expected = 'b3ae6d1'
if ((git rev-parse --short=7 HEAD).Trim() -ne $expected) {
    throw 'Local HEAD differs from b3ae6d1. Inspect latest local files before installation.'
}
$targets = @(
    'Numerical/plotting/FermionThresholdVariationFigures.py',
    'Numerical/plotting/DisplayModelLabels.py',
    'tools/RebuildCurrentReports.py',
    'tools/PlotDisplayWrapper.py',
    'scripts/RetireOldStudyDisplays.ps1',
    'docs/REGENERATE_CURRENT_REPORTS.md'
)
if (git status --porcelain -- $targets) {
    throw 'Some target files have local modifications. Review or back them up first.'
}
Expand-Archive -LiteralPath .\T3_REPORT_DISPLAY_FIX_b3ae6d1.zip -DestinationPath . -Force
python -m compileall -q Numerical\plotting\FermionThresholdVariationFigures.py Numerical\plotting\DisplayModelLabels.py tools\RebuildCurrentReports.py tools\PlotDisplayWrapper.py
```

Only install over a matching source: any modifications to these files in your **local** checkout should be reviewed first. The archive contains no numerical data and no symbolic artifacts.

## Rebuild using existing JSON (no physics reruns)

```powershell
python -m tools.RebuildReports --study fermion-threshold --study weinberg
python -m tools.RebuildReports --study fitted-mass
# Optional: rebuild all available current static figures and interactive reports
python -m tools.RebuildCurrentReports
```

`tools.RebuildCurrentReports --only plots` regenerates static PNGs only; `--only displays` regenerates HTML only. Missing optional saved datasets may cause `[SKIP]` messages; use `--strict` if missing datasets should be errors.

The two primary figure sets are stored at:

- `Reports/output/display/fermion_threshold_variation/smallY_smallL/fermion_threshold_running_effect_*.png`
- `Reports/output/display/weinberg_running/smallY_smallL/weinberg_running_*.png`

The existing fitted mass figures continue to be written to `Reports/output/full/fermion_scale/smallY_smallL/fermion_scale_*.png` to preserve the current project output convention. These are generated PNGs, not inputs.

The unchanged numerical sources are, respectively:

- `Reports/output/full/staged_rg_comparison/smallY_smallL/scan.json`
- `Reports/output/full/weinberg_running_comparison/smallY_smallL/scan.json`
- `Reports/output/full/fermion_scale/smallY_smallL/scan.json`

Their directory names remain historical because the numerical runners write there. They are **not** renamed, modified or deleted. The source rows have 60 successful points and four model-level `MissingArtifacts` rows without `mass_ratio` for each study in GitHub commit `b3ae6d1`. The new renderer excludes such model-level records from the ratio sort; numerical failure points with ratios still break plotted curves.

## What the figures mean

**Fermion-threshold variation — RGE running versus frozen segments** uses the horizontal variable `μ_F/μ_S` with both `M_F` and `μ_F` varied and `μ_S` fixed. It is **not** an energy-running trajectory.

For each observable at the final low-energy endpoint, the plotted incremental percentages are exactly the existing saved definitions:

- UV running enabled vs frozen: `000 → 100`, `100 × (O100/O000 − 1)`.
- EFT1 running enabled vs frozen: `100 → 110`, `100 × (O110/O100 − 1)`.
- Weinberg running enabled vs frozen: `110 → 111`, `100 × (O111/O110 − 1)`.

Matching remains in both compared calculations; no percentages are energy derivatives. The separate Weinberg scan compares evolved with frozen `C5` in the final EFT; the separately fitted fermion-scale scan compares against the same model's mass-factor-one benchmark, **not** against frozen RGEs.

Legend labels follow `Numerical/orchestration/ModelNames.py`: e.g. `T3A-m2`, `T3B-m1`, `T3E-p0`. The original keys, scan data, file/directory names and numerical array values are unchanged. Other static figures rebuilt through `tools.RebuildCurrentReports` have plot **text artists only** relabelled on save (also in the separate fitted sensitivity thesis plotting subprocess). Fitted and common-input envelope HTML plot selectors are also normalised to the same short spelling without modifying their numerical payload keys. Other HTML dashboards already have readable model names.

## Preview, apply and verify retirement

```powershell
.\scripts\RetireOldStudyDisplays.ps1
.\scripts\RetireOldStudyDisplays.ps1 -Apply -WhatIf
.\scripts\RetireOldStudyDisplays.ps1 -Apply

# Confirm new outputs and check that removed files stay retired
Get-ChildItem .\Reports\output\display\fermion_threshold_variation\smallY_smallL\*.png
Get-ChildItem .\Reports\output\display\weinberg_running\smallY_smallL\*.png
Test-Path .\Reports\output\display\study_comparisons
Test-Path .\Reports\output\display\rebuilt_studies
Test-Path .\Reports\output\display\interactive_correlations.html

# Confirm numerical sources are still intact
Get-Item .\Reports\output\full\staged_rg_comparison\smallY_smallL\scan.json
Get-Item .\Reports\output\full\weinberg_running_comparison\smallY_smallL\scan.csv
Get-Item .\Reports\output\full\fermion_scale\smallY_smallL\scan.json
```

Deletion is confined to old generated `study_comparisons/`, `rebuilt_studies/`, `display_stage5/rebuilt_studies/`, old `staged_rg_*.png` and duplicate source-folder `weinberg_running_*.png`, plus `interactive_correlations.html` in the retired display roots. The script verifies all 10 replacement PNGs exist and are nontrivial before `-Apply`, and refuses to delete old display directories if numerical scans unexpectedly appear inside them. Default invocation is **preview only**. Standalone archived physics study generators might still generate their own historical PNGs if explicitly run with `--run`; they are deliberately **not** invoked by the report rebuild and should not be run for plot-only changes.

Study A/B HTML generation is not called by the current rebuild. The historical source generator is retained, and numerical inputs under `Reports/output/full/common_scalar_rg_controls/` and `heavy_fermion_scale_comparison/` are untouched.
