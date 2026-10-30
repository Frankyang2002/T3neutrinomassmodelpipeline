# Zero-extra all-16: original interactive RGE and gamma 0.95–1.05

This is a local-only **additive** adapter for the saved `RunT3AllZeroExtraScalarControl` data. No numerical scalar RGE files, PDFs, or existing frozen-64 reports are modified.

## Source verified from GitHub, 2026-10-30, `main` f68ee94cd8ee

`output/studies/general_scalar_extended/zero_extra_all16/zero_extra_all16_low_v1/`

- `summary.json` has 16 models x 2 modes (32 rows).
- `optimised_low/held_zero/<tag>/continuous_trajectory.json` and matching `neutrino_v1/<tag>_neutrino.json`.
- `optimised_low/boundary_zero/<tag>/continuous_trajectory.json` and matching `neutrino_v1/<tag>_neutrino.json`.
- **Extra scalar couplings and holomorphic bilinear start at zero at 100 TeV**. Original `lambdaT3`, ordinary scalar interactions, diagonal masses and Yukawas remain nonzero.
- `held_zero`: prohibit radiative generation by zeroing the extra beta coordinates (restricted-model control).
- `boundary_zero`: allow generation using the candidate full RGBeta beta coordinates.
- Historical individually fitted, low-coupling, **frozen UV**, no full sequential matching or refit.

Do **not** confuse this with `output/scalar_relevance/zero_boundary_low_energy/`, which covers only B(-1)/C(-1) and uses different UV conditions and matching.

## Output files and project paths (seven new files, no replacements)

- `scripts/T3ZeroExtraReportSources.py`: exact 32-cell source reader and numerical provenance validation.
- `scripts/GenerateT3ZeroExtraWeinbergRunning.py`: computes **genuine final EFT** C5(μ), m1–m3 and five oscillation observables 1 TeV → 100 GeV with the *existing* Weinberg numerical solver; verifies saved 100 GeV endpoint.
- `scripts/GenerateT3ZeroExtraGamma.py`: new 33-point T1/T2/T3 scan for each model/mode; uses the *existing* frozen-64 perturbative C5/low-energy evaluator.
- `Numerical/plotting/BuildT3ZeroExtraReferenceReports.py`: reuses `BuildT3FullScalarFrozen64HTML.operator_html` and its **actual uploaded** `T3OptimisedReferenceRGE.html` canvas, plus `BuildFittedModelEnvelopes.HTML` for ovals; no new plotting UI.
- `scripts/RunT3ZeroExtraReferenceReports.ps1`: run all steps in order, reuse compatible cached results, open report.
- `tests/test_t3_zero_extra_reference_reports.py`: isolated read/plot/grid contract tests.
- `docs/T3_ZERO_EXTRA_REFERENCE_REPORTS.md`: these instructions.

## Install and run (PowerShell)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\T3_zero_extra_original_reference_reports.zip" -DestinationPath . -Force
python -X utf8 -m pytest -q tests\test_t3_zero_extra_reference_reports.py
.\scripts\RunT3ZeroExtraReferenceReports.ps1
```

If the source `zero_extra_all16_low_v1` is absent locally, run the *existing* campaign first:

```powershell
python -u -X utf8 -m scripts.RunT3AllZeroExtraScalarControl `
  --run-id zero_extra_all16_low_v1 --benchmark optimised_low `
  --models all --modes held_zero,boundary_zero
```

Generated source-level data (new folders **inside** the zero-extra campaign, no original data overwritten):

```text
output/studies/general_scalar_extended/zero_extra_all16/zero_extra_all16_low_v1/
  weinberg_running_v1/{held_zero,boundary_zero}/<model>.json
  gamma_095_105_v1/{held_zero,boundary_zero}/<model>.csv
  gamma_095_105_v1/{held_zero,boundary_zero}/<model>.json
```

Generated HTML:

```text
Reports/output/display_zero_extra_all16/zero_extra_all16_low_v1/
  index.html
  held_zero/interactive_rge.html
  held_zero/sensitivity_ovals.html
  held_zero/all_results.html
  boundary_zero/interactive_rge.html
  boundary_zero/sensitivity_ovals.html
  boundary_zero/all_results.html
  coverage.csv
  provenance.json
```

Scan γ = 0.95, 0.96, ..., 1.05 (11 values); 3 directions, 16 models, 2 modes = **1,056 actual neutrino evaluations**. An old 0.5–1.5 scan will never be read. The oval renderer is the same `BuildFittedModelEnvelopes.HTML`, including model selection, NuFIT/JUNO approximate envelopes, zoom, x/y observable controls and γ selection.

**Physics caveats:** The intermediate EFT `C12(mu)` Wilson trajectory is **not saved nor numerically solved** here. The original LLSS η selection remains an explicitly analytic **kernel**, not `C12(mu)`. B(-1)/C(-1) use conditional mass-basis one-scale C5 matching; where self-bilinears are regenerated, their C5 contributions are omitted by the existing restricted kernel; Matchete/RGBeta coupling mapping and full sequential EFT matching are not certified. RGE and ovals display actual conditional numerical results, not refits.

## Scope of validation

Source schema and the two distinct boundary modes were inspected at GitHub commit `f68ee94cd8ee`. Local isolated contract tests do not replace an actual 32-case SciPy Weinberg run and 1,056 gamma evaluations. Both are performed only by the provided local PowerShell command.
