# Frozen T3 RGE: user's `optimised_rge.html` layout and narrow sensitivities

Local-repository implementation only. Based on GitHub main SHA `519a0b3fa9965f161b9b044384bb41243cff8103` and the attached user file `optimised_rge.html`.

## Installed paths (exact)

- `Numerical/plotting/BuildT3FullScalarFrozen64HTML.py` **full replacement** of the previously generated, mismatched HTML adapter. Uses the bundled user-supplied template verbatim, not the new custom dashboard or the scalar reporter.
- `Numerical/plotting/templates/T3OptimisedReferenceRGE.html` **new file**, an exact structural extraction of uploaded `optimised_rge.html`: CSS, DOM, and canvas/JS controls are retained; the ~4 MB *historical* DATA/META arrays, menu options and title are replaced with placeholders for live frozen-64 results. No old numerics leak into the new runs.
- `scripts/GenerateT3Frozen64WeinbergRunning.py` **full file**, unchanged from the previous package; uses existing solver and validates against saved final SM+Weinberg endpoint.
- `scripts/RunT3Frozen64NarrowGammaReport.ps1` **new file**, runs/resumes the actual narrow sensitivity numerical scan, generates final Weinberg trajectories and rebuilds the HTML.
- `tests/test_t3_frozen64_weinberg_report.py` **full replacement** contract tests for exact template controls, original CSS, real-data adapter, strict gamma scope and outlier classification.
- `docs/T3_FROZEN64_REFERENCE_RGE_GAMMA095_105.md` **new file** (this document).

No changes are required to `InteractiveModelComparison.py`, `BuildFittedModelEnvelopes.py`, physical RGE implementations, PDF generators or saved historical numerical outputs.

## What the chart now looks like

The **same** header, side panel, quantity dropdown, model checkboxes, scale zoom buttons, 16 model colours and original `T3-A-m2`-style model ordering, canvas plotting geometry, legends, y-axis controls and browser wheel/double-click interactions as the uploaded file. It is not a newly designed graph. No inserted warning banners or large extra panels change the reference layout.

The selectable data are genuine saved numerical final-EFT `C5` matrix norm and six independent absolute components, neutrino mass eigenvalues, two mass-squared splittings, and three mixing angles. They run from **1000 GeV to 100 GeV**. The LLSS selector is retained as a clearly identified *analytic kernel*, NOT a computed C12 trajectory. No intermediate C12 coefficient or hard+direct C5 continuation is fabricated. The original JS zoom controls remain, but the initial zoom is set to the range where numerical data actually exist.

The header reports candidate matching scope and mass-scale outliers rather than presenting all 64 numerical solutions as physical neutrino fits. Existing unstable low-lambda predictions are not rescaled, replaced or hidden.

## Narrow sensitivity gamma

New 11-point grid:

`0.95,0.96,0.97,0.98,0.99,1.00,1.01,1.02,1.03,1.04,1.05`

Each of `T1`, `T2`, `T3` is recomputed at each grid point (33 points per model/benchmark, **2112 evaluations** across 64 slots, when all succeed). The generator and `BuildFittedModelEnvelopes.py` remain unchanged. The HTML adapter **rejects old 0.5–1.5 data** and checks exact gamma grid, per-row identities, source hashes and reproduction of the `gamma=1` central neutrino output. New scan ID: `gamma_095_105_v1`. The old scan `gamma_scans_v1` is neither overwritten nor used.

These sensitivity ovals describe local responses, not confidence intervals. Smaller gamma perturbations lead to smaller envelopes in general, but no improved experimental agreement is implied.

## PowerShell (local checkout)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline

git status --short
Expand-Archive `
  -LiteralPath "$HOME\Downloads\T3_frozen64_reference_RGE_gamma095_105.zip" `
  -DestinationPath . -Force

python -X utf8 -m pytest -q tests\test_t3_frozen64_weinberg_report.py

# Uses existing frozen_general64_v1 numerical source manifest and central results.
.\scripts\RunT3Frozen64NarrowGammaReport.ps1

# Open sensitivity for the same benchmark
Start-Process .\Reports\output\display_general_scalar_frozen_64\frozen_general64_v1\optimised_low\sensitivity_ovals.html
```

Alternatively run the narrow numerical scan explicitly, **without modifying** the original gamma generator:

```powershell
python -u -X utf8 -m scripts.GenerateT3FrozenFullScalarGamma `
  --run-id frozen_general64_v1 --scan-id gamma_095_105_v1 `
  --gamma-grid "0.95,0.96,0.97,0.98,0.99,1.00,1.01,1.02,1.03,1.04,1.05" `
  --resume

python -X utf8 -m Numerical.plotting.BuildT3FullScalarFrozen64HTML `
  --run-id frozen_general64_v1 --scan-id gamma_095_105_v1 `
  --weinberg-output-id weinberg_running_v1
```

The actual frozen-64 data and scientific statuses stay unchanged. The initial C5 endpoints can be numerically large because historical low-λ Yukawa seeds are not full-general scalar refits; the six self-bilinear cases omit a physical matching contribution. The numerical intermediate C12 evolution and complete general-scalar matching are not claimed or included.

## Validation

The local isolated contract tests verify layout preservation, numeric data mapping, gamma-grid rejection of old scans, and outlier classification. The rendered reference-template JavaScript is syntax-checked with Node. The actual project-wide 64 numerical running and sensitivity evaluations must be performed in the user's checkout; this environment has only the uploaded HTML and earlier generated source package, not the full local repository.
