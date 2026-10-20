# Reports/output/full — retention audit

**Repository basis:** read-only inspection of GitHub `main` at commit `166ed4c57b8c32eb0bb55e5b0cc9e30b67e7bb64` on 2026-10-20. This is **not** a scan of the user's Windows filesystem. Local untracked files/empty folders may differ.

## Supported physics input families

1. **Common-input:** UV ordinary couplings, heavy parameters and scales inherited from fitted `T3-A-m2` by `RunCommonInputRGE.py`, including the explicitly coded class-E half-`lambdaT3` convention.
2. **Individual fits:** independently optimised Yukawa matrices/parameters selected from feasibility, retry or focused fit candidates. Some studies share a scanned scalar coupling across these fits; that does **not** make them obsolete E-p0/common-Yukawa studies.

## Retain all current tracked numerical study families

| Under `Reports/output/full/` | Tracked files | Basis / why keep |
|---|---:|---|
| `common_input_rge/` | 145 | Common A-m2 reference: diagnostics and running figures |
| `common_input_sensitivity/` | 17 | Common-reference gamma scans used by current envelope HTML |
| `common_input_coupling_scan/` | 28 | Common-reference coupling variations, RGE/frozen results |
| `common_input_fermion_scale/` | 23 | Common-reference fermion mass/threshold scan |
| `fit_feasibility/` | 34 | Fitted configurations and fit selection inputs |
| `fit_feasibility_retry/` | 18 | Best-fit selection inputs (not obsolete solely because called retry) |
| `fit_feasibility_focused/` | 2 | Best-fit selection inputs (not obsolete solely because called focused) |
| `fitted_rge/` | 131 | Individually fitted running diagnostics and current figures |
| `fitted_sensitivity/` | 56 | Active wide-range fitted gamma scan and thesis figures |
| `fermion_scale/` | 7 | Individually fitted mass-ratio results; current renderer reads `scan.json` |
| `staged_rg_comparison/` | 2 | **Retain** saved 000,100,110,111 numerical endpoint comparisons; historical internal folder name |
| `weinberg_running_comparison/` | 2 | **Retain** saved Weinberg evolved/frozen endpoint comparison data |
| `common_scalar_comparison/` | 28 | **Individual fitted Yukawas**; a common lambda grid, not four historical benchmarks |
| `common_scalar_rg_controls/` | 12 | Running vs frozen control points corresponding to common scalar scans |
| `common_scalar_beta_dominance/` | 5 | Derived UV beta decomposition linked to current scalar study |
| `heavy_fermion_scale_comparison/` | 12 | **Individual fitted Yukawas** plus common-scalar points; retain numerical Study B data |
| `representation_comparison/` | 2 | Compact aggregation of fitted UV/EFT1/Weinberg endpoint contributions |
| `common_scalar_diagnostics/` | 2 | 80-point scalar study audit/summary; not used by current rebuild, but real data |
| `fitted_sensitivity_all/` | 18 | Earlier narrow `gamma=0.8...1.2`, 21-point fitted scan; not used by current rebuild, but distinct saved numerical samples; archive only if intentionally discarding that scan |

For the 19 tracked groups: **544 files total; zero 0-byte files and zero empty tracked directories**. These statistics do not say whether a Windows checkout has additional untracked/ignored files or empty directories.

## Retire (if present locally)

The four fixed benchmark and old shared-parameter/E-p0 report lineage:

- `four_benchmark_correlations/`
- `common_parameter_comparison/`, `common_parameter_comparison_v2/`
- `common_yukawa_optimization/`, `common_yukawa_least_squares/`, `common_yukawa_final_summary/`
- `common_yukawa_correlation_scans/`, `common_yukawa_three_panel_scans/`
- `scalar_coupling_correlations/`
- `benchmark_summary/`, `yukawa_grouping/`
- `comparison/`, `within_model_comparison/` (old umbrella/four-scenario lineage)

These 13 directory names are absent from the audited GitHub full tree; they may still exist in a local checkout. The accompanying script will **quarantine**, not irreversibly delete, those directories when `-Apply` is explicitly requested.

Also remove old `staged_rg_*.png` from `staged_rg_comparison/<scenario>/` and `weinberg_running_*.png` from `weinberg_running_comparison/<scenario>/` **only** when the corresponding five current display PNGs exist and are nonempty. Do not delete these studies' `scan.json` and `scan.csv`.

## Manual review (not automatic removal)

If found locally: `beta_dominance/`, `observable_sensitivity/`, `common_scalar_pilot/`, `hypercharge/`, `dimensions/`, `smoke/` and unknown study roots. The current main rebuild doesn't use them, but absence from the main workflow alone does not prove all contained data/artifacts are disposable. Check for unique saved numerical or symbolic content first.

Do **not** automatically remove `normalized_fitted_sensitivity/`: `tools/RebuildCurrentReports.py` currently generates it from `fitted_sensitivity/`, even though another historical cleanup script names it for retirement. Keeping or moving that derived output is a separate report-layout decision.

**Avoid running `scripts/CleanReportOutput.ps1 -Apply` for this selective job.** It intentionally deletes all PNGs and HTML under `Reports/output`, including actively retained figures, and its scope is much broader.

## Install and run (PowerShell at repository root)

```powershell
Expand-Archive -LiteralPath .\T3_REPORTS_FULL_AUDIT.zip -DestinationPath . -Force
.\scripts\AuditAndPruneReportsFull.ps1
.\scripts\AuditAndPruneReportsFull.ps1 -Apply -WhatIf
.\scripts\AuditAndPruneReportsFull.ps1 -Apply
```

The script is preview-only by default; deletion requires explicit `-Apply` (add `-Confirm` for individual confirmations). Retired trees are moved to `.cleanup_quarantine/reports_full/<timestamp>/` outside the report tree. Empty files are deleted and empty directories are removed bottom-up. The script does **not** delete unrecognized nonempty directories, current numerical studies, fitted benchmark JSON/CSV, or files outside `Reports/output/full` (aside from creating the quarantine for recoverability).

For an independent exact local inventory:

```powershell
Get-ChildItem .\Reports\output\full -Directory | Select-Object Name
Get-ChildItem .\Reports\output\full -Recurse -File -Force |
    Where-Object Length -eq 0 | Select-Object FullName
Get-ChildItem .\Reports\output\full -Recurse -Directory -Force |
    Where-Object { -not (Get-ChildItem -LiteralPath $_.FullName -Force | Select-Object -First 1) } |
    Select-Object FullName
```
