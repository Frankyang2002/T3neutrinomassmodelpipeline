# Frozen-64 RGE display — original interactive report (corrected)

Source verified: `Frankyang2002/T3neutrinomassmodelpipeline`, GitHub `main`, 30 October 2026 (`519a0b3fa9965f161b9b044384bb41243cff8103`). Work locally; **no GitHub modifications**.

## What this ZIP replaces

- **Full replacement:** `Numerical/plotting/BuildT3FullScalarFrozen64HTML.py`
- **New** (unchanged from the previous numerical package): `scripts/GenerateT3Frozen64WeinbergRunning.py`
- **Full replacement:** `tests/test_t3_frozen64_weinberg_report.py`
- **Full replacement:** `docs/T3_FROZEN64_NEUTRINO_RGE_REPLACEMENT.md`

The wrong bespoke HTML canvas and scalar-trajectory renderer are removed **from the frozen-64 report generator**. Its `interactive_rge.html` is now the **original** `Numerical/plotting/InteractiveModelComparison.py` HTML/CSS/JavaScript canvas, extracted directly from the `write_dashboard()` literal `template` using Python AST. The full original interaction logic remains: model checkboxes, show all/none, quantity selector, zoom buttons, scroll/double-click zoom and dynamic/fixed y-axis range. **Do not delete or overwrite `InteractiveModelComparison.py`.** Existing scalar numerical RGE machinery and PDF generators remain untouched. The unrelated historical scalar reporter `BuildGeneralScalarLambdaReports.py` is *not called* or imported by the corrected frozen-64 renderer; other existing workflows may still depend on it.

The only adapted inputs are the real saved numerical `C5` 3×3 symmetric complex magnitude matrix, its Frobenius norm, three neutrino masses and five oscillation observables. The default x-axis is set to the saved numerical interval **1 TeV → 100 GeV** instead of 10^7 GeV → 100 GeV; the original interactive controls can still zoom out. There are 16 selectable models per benchmark and four benchmark pages.

The existing frozen-64 model-source manifest, input hashes and numerical `C5(1 TeV)` / `C5(100 GeV)` endpoint checks are retained. The preceding `GenerateT3Frozen64WeinbergRunning.py` uses the project's existing `Numerical.running.WeinbergTrajectory.run_weinberg_trajectory` and reproduces previously saved 100 GeV outputs; it does not estimate a trajectory by interpolating between the two saved endpoints.

**Physics limitations:** intermediate EFT `C12`/`LLSS` running is **not available** in the frozen-64 results and is not synthesised from an analytic `η` kernel. The original interactive report's old hard+direct *display reconstruction* is likewise disabled for these cases: only the final 1 TeV → 100 GeV SM+Weinberg evolution is saved. A second, genuinely coupled `C12`/`C5` EFT1 numerical calculation would be necessary to plot `C12` between 100 TeV and 1 TeV. Full general-scalar matching remains candidate-level; six self-bilinear controls omit their mixing and B(-1)/C(-1) have conditional one-scale matching. The previously found extreme low-λ neutrino mass predictions remain labelled as outliers, not successful fits.

## Install and run — PowerShell

From the project root (change this path if your checkout is elsewhere):

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline

git status --short

Expand-Archive -LiteralPath "$HOME\Downloads\T3_frozen64_original_interactive_RGE.zip" -DestinationPath . -Force

python -X utf8 -m pytest -q tests\test_t3_frozen64_weinberg_report.py

# Required only if you have not yet calculated the 64 final-EFT trajectories.
python -u -X utf8 -m scripts.GenerateT3Frozen64WeinbergRunning `
    --run-id frozen_general64_v1 --output-id weinberg_running_v1

# Reuse previously calculated trajectories; render in the ORIGINAL interactive RGE canvas.
python -X utf8 -m Numerical.plotting.BuildT3FullScalarFrozen64HTML `
    --run-id frozen_general64_v1 `
    --scan-id gamma_scans_v1 `
    --weinberg-output-id weinberg_running_v1

Start-Process .\Reports\output\display_general_scalar_frozen_64\frozen_general64_v1\optimised_low\interactive_rge.html
```

If the numerical `weinberg_running_v1` output is already present and verified, **skip** regeneration (running the generator without `--resume` intentionally refuses overwriting existing output). Use `--resume` if you need to validate an existing numerical run, provided all source hashes and configuration match.

All four `interactive_rge.html` files are overwritten with the **original interactive report** upon rebuilding. Old scalar RGE report HTML is therefore replaced in place. No historical output data are deleted. This ZIP is safe to extract without overwriting the original interactive renderer.

## Verification scope

- Python syntax and isolated adapter tests can be checked here. The full existing repository is not mounted locally in the execution environment, so actual 64-case RGE integration and browser execution against repository-wide files must be run in the working checkout.
- This is **not** a fully physical C12+C5 general scalar EFT matching solution; no such claim is made.
