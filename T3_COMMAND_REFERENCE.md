# T3 neutrino-mass model pipeline — complete command reference

**Repository:** [Frankyang2002/T3neutrinomassmodelpipeline](https://github.com/Frankyang2002/T3neutrinomassmodelpipeline)  
**Audited revision:** [`7969f23e782549913d3bff640bdbc62179e3db51`](https://github.com/Frankyang2002/T3neutrinomassmodelpipeline/commit/7969f23e782549913d3bff640bdbc62179e3db51), `main`, committed 2026-10-19 22:34 UTC.  
**Intended location in your local project:** `docs/COMMAND_REFERENCE.md`  
**Basis:** Read-only inspection of the GitHub tree, current Python entry points, command-line parser definitions, and repository documentation on 2026-10-20. **Commands are source-verified, not runtime-tested on your Windows machine.** GitHub was not modified.

> **Important:** The current `pipeline.py --full` runs **16 frozen fitted-reference models** (`studies/FullT3Study.py`), not the former four-scenario 16×4 fixed-input sweep. Historical fixed comparison outputs may still exist under `output/full/comparison/`. `--full --numerical ...` is rejected. Some old README/PowerShell snippets point to modules that have been moved; see [Legacy pitfalls](#14-legacy-paths-and-known-pitfalls).

## Contents

1. [Shell setup, help, and conventions](#1-shell-setup-help-and-conventions)
2. [Most useful commands](#2-most-useful-commands)
3. [Master pipeline: symbolic matching and full runs](#3-master-pipeline-symbolic-matching-and-full-runs)
4. [Rebuild reports and dashboards from saved data](#4-rebuild-reports-and-dashboards-from-saved-data)
5. [Analytical, LaTeX and PDF reports](#5-analytical-latex-and-pdf-reports)
6. [Parameter scans, fitting, benchmark diagnostics](#6-parameter-scans-fitting-benchmark-diagnostics)
7. [Fitted RGE, threshold and mass-scale studies](#7-fitted-rge-threshold-and-mass-scale-studies)
8. [Common-input and common-Yukawa studies](#8-common-input-and-common-yukawa-studies)
9. [Common-scalar and beta-function studies](#9-common-scalar-and-beta-function-studies)
10. [Numerical diagnostics, audits, extraction](#10-numerical-diagnostics-audits-extraction)
11. [Group theory, RGE and matching research CLIs](#11-group-theory-rge-and-matching-research-clis)
12. [Tests and regression](#12-tests-and-regression)
13. [Input/output layout and safe workflow order](#13-inputoutput-layout-and-safe-workflow-order)
14. [Legacy paths and known pitfalls](#14-legacy-paths-and-known-pitfalls)
15. [File and command index](#15-file-and-command-index)

## 1. Shell setup, help, and conventions

All examples below are **PowerShell**, run **from the repository root** in the configured Python environment:

```powershell
Set-Location "C:\path\to\T3neutrinomassmodelpipeline"
python --version
python pipeline.py --help
python -m pytest -q
```

- `python -m package.module` executes a module as a CLI **only if** it has a runnable `__main__` block; do not try to execute implementation-only modules.
- Replace placeholder paths such as `PATH_TO_SCAN_CONFIG.json` with real local files. They are **templates, not filenames verified to exist**.
- `--run`, `--run-missing` and `--apply` generally switch dry-run/inventory commands into generating or changing files. `--force` can recompute or overwrite cached output; use deliberately.
- **Presentation-only** means it reads saved numerical physics results and renders images/HTML/TeX; it does *not* compute fresh ODE trajectories, fits, RGBeta or matching.
- **Numerical-only** generally means it uses saved matching/RGBeta JSON and evolves numerical states; it can be time-consuming but does *not* rerun the external symbolic toolchain.
- **Full symbolic** may require Matchete, Mathematica/WolframScript, RGBeta and LaTeX on the machine, besides the Python environment.
- There is **no pinned `requirements.txt`/`pyproject.toml` in the audited root**; this guide does not pretend that Python dependencies are installable with a verified single command.
- Script options are the source of truth. Run `python -m <module> --help` (or `python <script.py> --help`) before unusual invocations.
- Output paths are relative to repository root unless written otherwise. Some tools make absolute paths internally. Commands may fail if required saved artifacts are missing; such failure is **not** proof that the physics model fails.

### Useful help commands

```powershell
python pipeline.py --help
python scripts/run_comparison_study.py --help
python -m Numerical.plotting.BuildDisplayResults --help
python -m tools.RebuildReports --help
python -m Numerical.fitting.ScanCLI --help
python -m Numerical.orchestration.studies.fitted_rge.RunFittedSensitivityComparison --help
```

## 2. Most useful commands

| Goal | Command (from repository root) | Compute new physics? |
|---|---|---|
| Rebuild fixed-grid model display from saved diagnostics | `python -m Numerical.plotting.BuildDisplayResults --output Reports/output/display` | **No** |
| Rebuild fixed-grid dashboard **and** fitted correlation report | `python -m Numerical.plotting.BuildCombinedInteractiveReports` | **No** |
| Rebuild only staged-RG saved figures | `python tools/RebuildSavedStudyPlots.py --study staged-rg` | **No** |
| Rebuild only final-Weinberg saved figures | `python tools/RebuildSavedStudyPlots.py --study weinberg` | **No** |
| Rebuild available figures/reports in one wrapper | `python -m tools.RebuildReports` | **No** |
| Regenerate fitted gamma envelopes | `python -m Numerical.plotting.BuildFittedModelEnvelopes` | **No** |
| Extract CSV/JSON observables from fixed-grid runs | `python -m Numerical.diagnostics.ExtractFullComparisonResults` | **No** |
| Audit saved Weinberg continuation | `python -m Numerical.diagnostics.VerifyFullContinuation` | **No** |
| Inventory numerical-only resume | `python -m Numerical.orchestration.support.ResumePipelineNumericalResults` | **No** (default) |
| Resume missing legacy fixed-grid numerical jobs | `python -m Numerical.orchestration.support.ResumePipelineNumericalResults --run` | **Yes**, numerical |
| Inventory common-input comparison | `python -m Numerical.orchestration.studies.common_inputs.RunCommonInputRGE` | **No** (default) |
| Evaluate common-input comparison | `python -m Numerical.orchestration.studies.common_inputs.RunCommonInputRGE --run` | **Yes**, numerical |
| Inventory fitted sensitivity scans | `python -m Numerical.orchestration.studies.fitted_rge.RunFittedSensitivityComparison` | **No** (default) |
| Evaluate fitted sensitivity | `python -m Numerical.orchestration.studies.fitted_rge.RunFittedSensitivityComparison --run` | **Yes**, numerical |
| Fresh symbolic 5-model smoke run | `python pipeline.py --smoke` | **Yes**, full symbolic |
| Run all 16 frozen fitted-reference models | `python pipeline.py --full` | **Yes**, full symbolic + numerical |
| Test Python code | `python -m pytest -q` | Tests; no routine study run |

**Recommended for figures first:** use the presentation-only commands before considering `pipeline.py --full`. These read what is already saved. Avoid full symbolic runs merely to regenerate plots.

## 3. Master pipeline: symbolic matching and full runs

**Entry point:** [`pipeline.py`](https://github.com/Frankyang2002/T3neutrinomassmodelpipeline/blob/7969f23/pipeline.py), parser implemented in [`common/PipelineCLI.py`](https://github.com/Frankyang2002/T3neutrinomassmodelpipeline/blob/7969f23/common/PipelineCLI.py).

### Core modes

```powershell
# Help / inspect available selectors
python pipeline.py --help

# Five historical T3 regression models
python pipeline.py --smoke

# Explore hypercharge alpha=-4 ... 2 across T3-A...E
python pipeline.py --hypercharge-comparison

# Compare dimensions at alpha=0, with neutrality filtering
python pipeline.py --dimension-comparison

# Full: 16 frozen fitted-reference runs (current implementation)
python pipeline.py --full

# Ordinary T3, one representation: dS1=2, dS2=2, dF=1, alpha=-1
python pipeline.py --dims 2 2 1 --alpha -1 --study B_example

# Shared-scalar mode: two dimensions (dS, dF), alpha fixed to -1
python pipeline.py --dims 2 1 --study scotogenic_example
```

**`--full` specifics (current code):** `studies/FullT3Study.py` reads `Numerical/benchmarks/fitted_reference/benchmark_manifest.json`, checks all 16 frozen input files and SHA-256 values, then calls `pipeline.py` once per model with `--reset-numerical-configs` and fermion-first thresholds by default. It writes `output/full/fitted_reference/full_fitted_reference_summary.json`, plus individual model data/reports. It does **not** optimize new Yukawa textures and it refuses an accompanying `--numerical` argument.

### Numerical and threshold options

```powershell
# Use default numerical config for a single model
python pipeline.py --dims 2 2 1 --alpha -1 --numerical `
  --threshold F --threshold S1 S2

# Supply your own numerical config explicitly
python pipeline.py --dims 2 2 1 --alpha -1 `
  --study B_numeric `
  --numerical "configs\t3_numerical_default_B_alpha_m1.json" `
  --threshold F --threshold S1 S2

# Regenerate model-specific configs from supplied template
# WARNING: replaces generated numerical configs
python pipeline.py --dims 2 2 1 --alpha -1 `
  --numerical "configs\t3_numerical_default_B_alpha_m1.json" `
  --reset-numerical-configs `
  --threshold F --threshold S1 S2

# Explicit symbolic threshold scales, one scale per threshold group
python pipeline.py --dims 2 2 1 --alpha -1 `
  --threshold F --threshold S1 S2 `
  --threshold-scale MF --threshold-scale MS

# Verbose expression reports and Wolfram logs (expensive)
python pipeline.py --dims 2 2 1 --alpha -1 --debug-reports
```

| Argument | Verified purpose |
|---|---|
| `--smoke` | Five regression models (mutually exclusive with the other run selectors). |
| `--hypercharge-comparison` | Scan alpha `-4` to `2` over A–E representations; normal selection requires neutral BSM content. |
| `--dimension-comparison` | Compare A–E at alpha=0 with neutrality filter. |
| `--full` | Call the **frozen fitted-reference 16-model** full-study runner. |
| `--dims DS1 DS2 DF` | Ordinary three-irrep T3 model. |
| `--dims DS DF` | Physical shared-scalar mode; requires `alpha=-1`. |
| `--alpha N` | Hypercharge parameter (default 0 for ordinary three-dimension mode). |
| `--study NAME` | Subfolder for raw and report outputs. |
| `--numerical [JSON]` | Include numerical evolution. No value uses `configs/t3_numerical_default_B_alpha_m1.json`; omitted entirely means symbolic only. |
| `--reset-numerical-configs` | Recreate per-model configs from template; **overwrites/replaces previously generated configs**. |
| `--threshold F --threshold S1 S2` | Ordered decoupling groups: heavy fermion, then both scalars together; the documented fermion-first numerical production hierarchy. |
| `--threshold-scale SCALE` | Repeat to supply one matching scale per threshold group. Supports symbolic/numeric values. |
| `--debug-reports` | Keep logs and produce extended UV/EFT expression reports. |
| `--force` | Bypass neutral-component and SU(2)-dimension support filters; **not** a guarantee that unsupported EFT sequences will work. |

**Threshold caution:** `--threshold` can be repeated; fields in one occurrence decouple simultaneously. Without a threshold declaration the parser groups all fields together. The execution planner rejects unsupported production ordering; do not assume that arbitrary threshold permutations are implemented. The ordinary fermion-first numerical path is `--threshold F --threshold S1 S2`.

**Output:** `output/<study>/` and `Reports/output/<study>/`; named outputs below `full/` when `--full` or a `--study full/...` path is supplied.

### One model / one old fixed-grid scenario

**Script:** [`scripts/run_comparison_study.py`](https://github.com/Frankyang2002/T3neutrinomassmodelpipeline/blob/7969f23/scripts/run_comparison_study.py).

```powershell
python scripts/run_comparison_study.py `
  --dims 3 3 2 --alpha 0 --scenario largeY_largeL

# Select thresholds explicitly and keep additional debugging reports
python scripts/run_comparison_study.py `
  --dims 3 3 2 --alpha 0 --scenario smallY_smallL `
  --threshold F --threshold S1 S2 --debug-reports
```

Required: `--dims DS1 DS2 DF`, `--alpha`, `--scenario` in `smallY_smallL`, `smallY_largeL`, `largeY_smallL`, `largeY_largeL`. Optional: `--numerical JSON`, repeated `--threshold`, repeated `--threshold-scale`, `--debug-reports`, `--force`. **Unlike the full fitted-reference mode**, this script constructs one *fixed comparison* point from the historical `(Y, lambda)` scenario grid, saves its config under `configs/generated_model_sets/comparison_test/<scenario>/`, and launches the pipeline with `--reset-numerical-configs`. Uses fermion-first thresholds by default. This performs new physics calculations.

### Auxiliary full-study wrappers

| Command | Behavior / caution |
|---|---|
| `python scripts/FullStudyOrganised.py` | Runs `pipeline.py --full`, then moves/consolidates symbolic reports into model-centric layout. |
| `python scripts/FullStudyOrganised.py --keep-originals` | Same, but preserve original layout when migrating reports. |
| `python scripts/FullStudyOrganised.py -- --debug-reports` | Pass additional pipeline arguments after `--`. |
| `python scripts/FullStudyShortNames.py` | Runs full pipeline with short study names. **Caution:** it monkey-patches a legacy `_model_key` attribute while the current full-study implementation uses frozen manifest paths; do not rely on it to rename outputs without local verification. |
| `powershell -ExecutionPolicy Bypass -File scripts/FullRunWithWithinModelReports.ps1` | **Broken legacy post-step:** invokes `Numerical.orchestration.BuildComparisonReports`, which no longer exists at that path. Do **not** use as-is. Use the standalone correct module in §4 instead. |

## 4. Rebuild reports and dashboards from saved data

### Fixed-grid display and cross-model dashboards

```powershell
# Rebuild saved fixed-grid comparison plots and interactive dashboards
python -m Numerical.plotting.BuildDisplayResults `
  --data-input "output\full" `
  --output "Reports\output\display"

# Build display + saved fitted-sensitivity correlation dashboard
python -m Numerical.plotting.BuildCombinedInteractiveReports `
  --data-input "output\full" `
  --report-input "Reports\output\full" `
  --correlation-source "Reports\output\full\fitted_sensitivity" `
  --output "Reports\output\display_combined"

# Single saved-model study directory -> HTML
python -m Numerical.plotting.InteractiveModelComparison `
  "output\full\comparison\smallY_smallL\MODEL_FOLDER" `
  "Reports\output\display\single_model.html"

# Refresh all legacy fixed-grid cross-model comparisons from saved diagnostics
python -m Numerical.plotting.InteractiveModelComparison --refresh-full-comparison
```

`BuildDisplayResults` defaults: `--data-input output/full`, `--report-input Reports/output/full` (**accepted for compatibility, but not read**), `--output Reports/output/display`. The source is `output/full/comparison/<scenario>/**/running_diagnostics.json`; it does **not** generate new physics. `BuildCombinedInteractiveReports` additionally reads fitted-sensitivity scans, which must be available or its correlation stage can fail. `InteractiveModelComparison` accepts either two positional paths (`study_dir`, `output`) or `--refresh-full-comparison`, plus `--title`, `--input-root`, `--report-root`, `--display-root`, `--no-display-copy`.

### One wrapper for presentation rebuilds

```powershell
# Rebuild all studies for which saved inputs exist
python -m tools.RebuildReports

# Only selected reports
python -m tools.RebuildReports --study comparison --study correlations
python -m tools.RebuildReports --study staged-rg --study weinberg

# Change the scenario and location of the rebuilt outputs
python -m tools.RebuildReports `
  --scenario smallY_smallL `
  --data-root "output\full" `
  --correlation-source "Reports\output\full\fitted_sensitivity" `
  --output "Reports\output\display_stage5"
```

Choices for repeated `--study`: `comparison`, `staged-rg`, `weinberg`, `correlations`. Default output is `Reports/output/display_stage5/`. Absent optional saved studies are skipped when running without explicit selection; explicitly requested missing studies cause a nonzero exit. Outputs never replace raw input `output/full`.

### Rebuild *only* staged-RG / Weinberg figures

```powershell
python tools/RebuildSavedStudyPlots.py --study staged-rg
python tools/RebuildSavedStudyPlots.py --study weinberg

python tools/RebuildSavedStudyPlots.py `
  --study staged-rg --scenario smallY_smallL `
  --source "Reports\output\full\staged_rg_comparison\smallY_smallL\scan.json" `
  --output-dir "Reports\output\display\rebuilt_studies\staged-rg\smallY_smallL"
```

`--study` is required (`staged-rg` or `weinberg`); `--scenario`, `--source`, `--output-dir` are optional. Reads existing `scan.json` and refuses empty/non-success scans. It **never re-evaluates** a scan.

### Saved gamma scans: envelopes, correlations, comparisons

```powershell
# Fitted model envelope dashboard with browser-side gamma range filtering
python -m Numerical.plotting.BuildFittedModelEnvelopes `
  --scenario smallY_smallL `
  --output "Reports\output\display\fitted_model_envelopes_all.html"

# Explicit scan input; include warning rows if wanted
python -m Numerical.plotting.BuildFittedModelEnvelopes `
  --source "Reports\output\full\fitted_sensitivity\smallY_smallL" `
  --output "Reports\output\display\fitted_envelopes_with_warnings.html" `
  --include-warnings

# Fitted sensitivity HTML correlation explorer
python -m Numerical.plotting.BuildInteractiveCorrelations `
  --source "Reports\output\full\fitted_sensitivity" `
  --target-config "configs\nufit_6_1_ic24_no.json" `
  --output "Reports\output\display\interactive_correlations.html"

# Study A/B comparison dashboard from saved running-vs-fixed CSVs
python -m Numerical.plotting.BuildStudyABInteractive --study both
python -m Numerical.plotting.BuildStudyABInteractive --study A `
  --scenario smallY_smallL --study-a-csv "PATH_TO_STUDY_A.csv" `
  --output-root "Reports\output\display\study_comparisons"

# Saved gamma sensitivity as NuFIT sigma-unit shifts
python -m Numerical.plotting.NormalizedFittedSensitivity --scenario smallY_smallL

# Plot per-model fitted gamma scan summaries
python -m Numerical.plotting.PlotFittedSensitivityFigures

# Plot existing observable-plane sensitivity diagnostics (T1/T2/T3)
python -m Numerical.plotting.ObservableSensitivityDashboard
```

| Module | Optional arguments | Default input/output | Key behavior |
|---|---|---|---|
| `BuildFittedModelEnvelopes` | `--scenario`, `--source`, `--output`, `--include-warnings` | `Reports/output/full/fitted_sensitivity/<scenario>` → `Reports/output/display/fitted_model_envelopes_all.html` | Saved gamma envelopes; no solver. |
| `BuildInteractiveCorrelations` | `--source`, `--output`, `--target-config` | `Reports/output/full/fitted_sensitivity` → `Reports/output/display/interactive_correlations.html` | Interactive correlations with NuFIT target where available. |
| `BuildStudyABInteractive` | `--study A\|B\|both`, `--scenario`, `--study-a-csv`, `--study-b-csv`, `--output-root`, `--reference-json` | `Reports/output/display/study_comparisons/<scenario>/StudyA or StudyB/` | Endpoint/prediction-space Study A/B comparisons; reference ellipses only if reference JSON supplied. |
| `NormalizedFittedSensitivity` | `--scenario`, `--input-dir`, `--output-dir` | Input `Reports/output/full/fitted_sensitivity/<scenario>` → `Reports/output/full/normalized_fitted_sensitivity/<scenario>` | Finite-difference shifts in NuFIT-sigma units. |
| `PlotFittedSensitivityFigures` | `--input`, `--output` | Input `Reports/output/full/fitted_sensitivity/smallY_smallL`; output chosen by script if omitted | PNG plots + `sensitivity_plot_summary.csv`. |
| `ObservableSensitivityDashboard` | `--input-root`, `--output-root` | `output/full/comparison` → `Reports/output/full/observable_sensitivity` | Reads `sensitivity_diagnostics.json`; does not generate missing points. |

### Within-model and three-panel displays

```powershell
# Four-scenario within-model reports from saved comparison diagnostics
python -m Numerical.orchestration.archive_candidates.BuildComparisonReports

# Require all four fixed-grid scenarios for each model
python -m Numerical.orchestration.archive_candidates.BuildComparisonReports --strict

# Alternate saved input/output roots
python -m Numerical.orchestration.archive_candidates.BuildComparisonReports `
  --input "output\full\comparison" `
  --output "Reports\output\full\within_model_comparison"

# Saved common-Yukawa three-panel visualization (reads a pre-existing CSV)
python -m Numerical.orchestration.presentation.BuildThreePanelInteractive
python -m Numerical.orchestration.presentation.BuildThreePanelInteractive `
  --csv "PATH_TO_scan_points.csv" `
  --output "display\my_three_panel.html"

# Plot scalar-only LLSS Wilson fixed-order diagnostic (one-generation caveat)
python -m Numerical.plotting.LLSSWilsonRunningReport --scenario smallY_smallL
```

`BuildComparisonReports` is located under **`archive_candidates`**, but the CLI is still implemented; it is **presentation-only**. Its defaults are `output/full/comparison` → `Reports/output/full/within_model_comparison`. `BuildThreePanelInteractive` defaults to `display/T3_common_yukawa_three_panel_interactive.html` and refuses an already-existing output file; supply a new `--output` name if needed. `LLSSWilsonRunningReport` additionally accepts `--input-root`, `--output-root`, `--model`, `--mu-f` (default `1e5` GeV), `--mu-s` (default `1e3` GeV), `--values-json`, `--force`. Its LLSS plot is a **fixed-order component diagnostic, not the complete matched `C5` evolution**.

### Thesis figure indexing

```powershell
# Index chosen existing figures, without copying
python -m Numerical.orchestration.archive_candidates.CollectThesisFigures

# Index all eligible figures
python -m Numerical.orchestration.archive_candidates.CollectThesisFigures --all

# Copy selected existing figures into a thesis review location
python -m Numerical.orchestration.archive_candidates.CollectThesisFigures `
  --copy --destination "Reports\figure_review"
```

This tool **locates/copies existing figures**, not generates them. `--all` includes broader selections; `--copy` changes the filesystem.

## 5. Analytical, LaTeX and PDF reports

These rebuild or modify **derived report files**, not raw RGE/matching outputs. LaTeX compilation generally needs `pdflatex` accessible on `PATH`. The main pipeline also generates stage-dependent analytical reports through internal `Reports/` functions; those internal report classes have **no standalone CLI**.

```powershell
# FIRST: check that the 16 old fixed-grid analytical records are present
python scripts/BuildCombinedAnalyticalReports.py --check

# Build cross-model analytical reports from saved data
python scripts/BuildCombinedAnalyticalReports.py `
  --scenario smallY_smallL `
  --output "Reports\output\analytical"

# Another scenario or an explicit previously saved folder
python scripts/BuildCombinedAnalyticalReports.py `
  --scenario largeY_smallL `
  --source "output\full\comparison\largeY_smallL" `
  --output "Reports\output\analytical"

# Legibility-oriented C5-only TeX report from existing Matchete report
python scripts/BuildCompactC5Report.py
python scripts/BuildCompactC5Report.py --compile

# Apply LaTeX presentation cleanup to existing report sources
python scripts/PolishAnalyticalReports.py --root "Reports\output\analytical"
python scripts/PolishAnalyticalReports.py --compile

# Insert vetted mathematical definitions and audit unresolved symbols
python scripts/AddMathematicalDefinitions.py
python scripts/AddMathematicalDefinitions.py --compile
```

| Script | Arguments | What it changes |
|---|---|---|
| `scripts/BuildCombinedAnalyticalReports.py` | `--scenario` (one of four historical fixed scenarios), `--source`, `--output` (default `Reports/output/analytical`), `--check` | `--check` **only validates** saved model records; otherwise generates cross-model analytical TeX/PDF reports and **automatically archives** selected temporary/superseded report files after generation. |
| `scripts/BuildCompactC5Report.py` | `--root` (default `Reports/output/analytical/Lagrangian`), `--compile` | Reconstructs a readable C5-only report from saved Matchete source; refuses to conflate unequal kernels. |
| `scripts/PolishAnalyticalReports.py` | `--root` (default `Reports/output/analytical`), `--compile` | Edits **presentation TeX files in place**; optional PDF compilation. |
| `scripts/AddMathematicalDefinitions.py` | `--root` (default `Reports/output/analytical`), `--compile` | Edits generated TeX, writes `mathematical_definitions_audit.json`, and explicitly records unresolved definitions instead of inventing formulae. |

**Automatic report cleanup:** `scripts/AnalyticalReportHousekeeping.py` implements an importable `clean()` function but has **no standalone CLI**; `BuildCombinedAnalyticalReports.py` invokes it after successful regeneration. It limits cleanup to `Reports/output/analytical`, moves build trash/backups into `Reports/report_cleanup_archive/`, and leaves raw `output/` alone. Finished TeX/PDF sources are retained.

### Folder naming / migration (not physics)

```powershell
# Preview rename of existing *presentation* folders (dry-run)
python scripts/RenameT3ModelFolders.py

# Actually rename eligible presentation folders
python scripts/RenameT3ModelFolders.py --apply
```

This changes **folder layout** (`T3-A-m4`, etc.), not model representation data. Use `--apply` only when downstream readers support the new path layout. `FullStudyOrganised.py` (§3) can also move/organize reports but launches a fresh full study first.

## 6. Parameter scans, fitting, benchmark diagnostics

### Numerical scan from a JSON configuration

```powershell
# Validate scan configuration, paths, parameter bindings, threshold and JSON contracts
python -m Numerical.fitting.ScanCLI "configs\t3_scan_example.json" --dry-run

# Execute the scan defined in the JSON (computes numerical trajectories)
python -m Numerical.fitting.ScanCLI "configs\t3_scan_example.json"
```

`ScanCLI` accepts **one required positional** `config` JSON and optional `--dry-run`. The JSON controls model representation, base UV state, parameter bounds/grid or random strategy, matching/RGBeta artifacts, oscillation target, numerical scales, checkpoint and output. CLI **rejects shared-scalar scans** in the current adapter rather than silently using ordinary-model formulas. The existing `configs/t3_scan_example.json` is an example config, **not a guarantee all external artifact paths are available on every machine**. `--dry-run` first is strongly advised.

**Output:** JSON result and checkpoint at paths **defined inside that scan config**. Saved successful points can be visualized without recomputation:

```powershell
python -m Numerical.plotting.PlotScanResults `
  "PATH_TO_scan_result.json" `
  --target "configs\nufit_6_1_ic24_no.json"

# Optional alternate plot location
python -m Numerical.plotting.PlotScanResults `
  "PATH_TO_scan_result.json" `
  --target "configs\nufit_6_1_ic24_no.json" `
  --output-dir "output\plots\my_scan"
```

`PlotScanResults` generates `best_chi2_progress.png`, `chi2_distribution.png`, `best_point_pulls.png`, `best_point_observable_ratios.png`, parameter-vs-chi2 plots and `summary.txt`. Default plot root: `output/plots/<scan-file-stem>/`.

### Local oscillation optimizer

```powershell
# Start from an existing successful scan result
python -m Numerical.fitting.OscillationOptimizer `
  "configs\t3_scan_example.json" `
  --seed-result "PATH_TO_scan_result.json" `
  --active-count 5 --sensitivity-step 0.02 --max-nfev 50 `
  --output "output\my_optimization.json"

# Restart from the best point of a previous optimizer run
python -m Numerical.fitting.OscillationOptimizer `
  "configs\t3_scan_example.json" `
  --seed-optimizer-result "output\my_optimization.json" `
  --output "output\my_optimization_restart.json"

# Plot an existing optimizer result without rerunning the fit
python -m Numerical.plotting.PlotOptimizerResults "output\my_optimization.json"
```

| Optimizer argument | Meaning / default |
|---|---|
| `scan_config` (positional) | Scan config defining valid model, parameter bindings and optimization bounds. |
| `--seed-result JSON` | Start from previously successful scan point. |
| `--seed-optimizer-result JSON` | Warm-start from previous optimizer best fit (alternative seed mode). |
| `--seed-index N` | Select successful seed point from `--seed-result`; default saved best point. Incompatible with optimizer seed. |
| `--active-count N` | Optimize top N locally sensitive parameters; **default 5**. |
| `--sensitivity-step X` | Finite-difference step in normalized [0,1] coordinate box; **default 0.02**. |
| `--max-nfev N` | SciPy bounded least-squares evaluation limit; **default 50**. |
| `--output JSON` | **Required** output path. |
| `--checkpoint JSON` | Optional; defaults to `<output-stem>.checkpoint.json` alongside output. |

`PlotOptimizerResults` accepts `optimizer_json` and optional `--output-dir` (default `output/plots/<file-stem>/`). It creates chi-square trajectories, final pulls, sensitivity rankings, active parameter trajectories and `summary.txt`; no physics rerun.

### Full fitted-point physics and running diagnostics

```powershell
# Re-evaluate an optimized point and make neutrino / threshold diagnostics
python -m Numerical.diagnostics.BestFitDiagnostics `
  "configs\t3_scan_example.json" `
  --optimizer-result "output\my_optimization.json" `
  --output-dir "output\my_best_fit_physics" `
  --n-scale-points 120

# Running trajectory diagnostics in one JSON
python -m Numerical.diagnostics.RunningDiagnostics `
  "configs\t3_scan_example.json" `
  --optimizer-result "output\my_optimization.json" `
  --output "output\my_running_diagnostics.json"

# Create paper/thesis summary figures FROM existing best-fit diagnostic JSON
python -m Numerical.plotting.ThesisResultFigures `
  "PATH_TO_best_fit_physics.json" `
  --output-dir "output\thesis_results\my_fit"

# Create scale-resolved running figures FROM existing running diagnostics
python -m Numerical.plotting.RunningResultFigures `
  "output\my_running_diagnostics.json" `
  --output-dir "output\plots\my_running"
```

`BestFitDiagnostics` and `RunningDiagnostics` **recalculate the numerical trajectory** using stored optimizer parameters and existing symbolic artifacts; `--n-scale-points` defaults to 120. The two plotting modules are **read-only physics consumers**. `RunningResultFigures` requires both a diagnostics JSON and `--output-dir`. `ThesisResultFigures` accepts a best-fit JSON and optional `--output-dir` (default `output/thesis_results/<input-stem>`).

### Benchmark provenance and grouping

```powershell
# Analyze similarity of 16 frozen fitted Yukawa reference matrices
python -m Numerical.orchestration.RunYukawaGrouping

# Tune similarity grouping threshold
python -m Numerical.orchestration.RunYukawaGrouping --threshold 0.1

# OPTIONAL: write unvalidated substitution candidates; not physical fits
python -m Numerical.orchestration.RunYukawaGrouping --write-candidates
```

Options `--manifest` (default frozen fitted-reference manifest), `--output` (default script-defined grouping directory), `--threshold` (default 0.1), `--write-candidates`. Candidate files are deliberately placed under `UNVALIDATED_reference_candidates`; **do not** treat them as independently fitted benchmark inputs.

## 7. Fitted RGE, threshold and mass-scale studies

**Common convention:** Unless stated otherwise, these commands default to an **inventory / dry-run**. They use stored fitted candidates and saved RGBeta/Matchete expressions, not fresh symbolic derivations. The default scenario in this group is `smallY_smallL`.

### Standard fitted-RGE display (scale-resolved trajectories)

```powershell
# Inventory completed/available fitted models
python -m Numerical.orchestration.studies.fitted_rge.RunFittedRGEDisplay

# Numerical-only pilot of first two models
python -m Numerical.orchestration.studies.fitted_rge.RunFittedRGEDisplay --run --limit 2

# Run all missing scale-resolved trajectories; include HTML/figures
python -m Numerical.orchestration.studies.fitted_rge.RunFittedRGEDisplay --run
```

Options: `--scenario`, `--run`, `--limit` (0 = all), `--max-fit-chi2` (default 0.01), `--force` (recompute cached results), `--no-figures`. Reuses `fitted_candidate_config.json` and comparison matching/RGBeta artifacts; never refits. Writes numerical results under `Reports/output/full/fitted_rge/<scenario>/` and interactive dashboards under `Reports/output/display/`.

### Fitted T1/T2/T3 gamma-sensitivity comparison

```powershell
# Inventory
python -m Numerical.orchestration.studies.fitted_rge.RunFittedSensitivityComparison

# Fast pilot, including gamma=1 reference
python -m Numerical.orchestration.studies.fitted_rge.RunFittedSensitivityComparison `
  --run --limit 2 --gamma-min 0.5 --gamma-max 2 --points 7

# Full fitted sensitivity study (default 25 log-spaced gamma values)
python -m Numerical.orchestration.studies.fitted_rge.RunFittedSensitivityComparison --run
```

Options: `--scenario`, `--run`, `--limit` (0 all), `--force`, `--gamma-min` (0.001), `--gamma-max` (1000), `--points` (25), `--spacing log|linear` (log), `--max-baseline-chi2` (0.01). Gamma=1 is included; baseline chi-square is a **verification threshold**, not a new fitting constraint. Failed nonbaseline deformations remain marked failed rather than synthesized. Outputs: `Reports/output/full/fitted_sensitivity/<scenario>/`, including `comparison_summary.json`, `.csv`, per-model results and plots. Feed these to §4 plotters.

### Heavy-fermion mass and threshold scan around fitted points

```powershell
python -m Numerical.orchestration.studies.fermion_scales.RunFittedFermionScaleScan

python -m Numerical.orchestration.studies.fermion_scales.RunFittedFermionScaleScan `
  --run --limit 2 --factors 0.5 1.0 2.0

python -m Numerical.orchestration.studies.fermion_scales.RunFittedFermionScaleScan --run
```

Options: `--scenario`, `--factors` (default `0.5 0.75 1.0 1.5 2.0`), `--run`, `--limit` (0 all). Scales **both** the heavy fermion mass matrix `MF` and its fermion decoupling scale together, while preserving other UV inputs and lower thresholds. This is **not** a pure RGE scale-choice uncertainty study. Produces stored scan JSON/CSV and figures (consult the study's printed `Wrote:` path).

### Staged RGE ablation: UV, scalar-only EFT1, SM Weinberg

```powershell
# Inventory, do not integrate
python -m Numerical.orchestration.studies.fitted_rge.RunFittedStagedRGComparison

# Pilot
python -m Numerical.orchestration.studies.fitted_rge.RunFittedStagedRGComparison `
  --run --limit 2 --factors 1.0

# Full study over default heavy-mass factors
python -m Numerical.orchestration.studies.fitted_rge.RunFittedStagedRGComparison --run
```

Options: `--scenario`, `--factors`, `--limit` (0 all), `--run`. Writes `Reports/output/full/staged_rg_comparison/<scenario>/scan.json` and plots/CSVs. Defines the controlled cumulative states:

| Label | UV RGEs | Intermediate scalar-only RGEs | Final SM+Weinberg RGEs |
|---|---|---|---|
| `000` | frozen | frozen | frozen |
| `100` | running | frozen | frozen |
| `110` | running | running | frozen |
| `111` | running | running | running |

**Physical caveat:** Disabling one RGE region still retains **threshold matching**. Increments depend on this ordering and shared low-energy charged-lepton basis; do not call them separate diagram orders or independent physical predictions.

### Isolate the final Weinberg running contribution

```powershell
python -m Numerical.orchestration.studies.fitted_rge.RunFittedWeinbergRunningComparison

python -m Numerical.orchestration.studies.fitted_rge.RunFittedWeinbergRunningComparison `
  --run --limit 2 --factors 0.5 1.0 2.0

python -m Numerical.orchestration.studies.fitted_rge.RunFittedWeinbergRunningComparison --run
```

Options: `--scenario`, `--factors`, `--run`, `--limit` (0 all). Requires relevant fitted/artifact-aware mass-scale inputs; compares evolved vs frozen final `C5` while keeping the UV/EFT1 running, thresholds and charged-lepton basis consistent. Writes `Reports/output/full/weinberg_running_comparison/<scenario>/scan.json` with comparison figures. It isolates **final Weinberg coefficient running**, not all loops/RGE effects.

## 8. Common-input and common-Yukawa studies

The **common-input** studies apply a shared numerical input across different representations and compare output differences. They are not independent per-model fits. Current common-input RGE uses the saved fitted **A, alpha=-2** point as a reference and preserves the Class-E `lambdaT3` benchmark factor 1/2 by default.

### Common-UV-input RGE across models

```powershell
# Read-only inventory
python -m Numerical.orchestration.studies.common_inputs.RunCommonInputRGE

# Two-model numerical trial
python -m Numerical.orchestration.studies.common_inputs.RunCommonInputRGE --run --limit 2

# Run all eligible models
python -m Numerical.orchestration.studies.common_inputs.RunCommonInputRGE --run

# Recompute existing results (use carefully)
python -m Numerical.orchestration.studies.common_inputs.RunCommonInputRGE --run --force
```

Options: `--scenario` (**only `smallY_smallL` accepted**), `--run`, `--limit` (0 all), `--force`, `--no-figures`, `--raw-class-e-lambda`. Default E coupling is 0.05 when the reference is 0.1; `--raw-class-e-lambda` bypasses that normalization. Output `Reports/output/full/common_input_rge/<scenario>/`, plus `Reports/output/display/common_input_rge_<scenario>.html`.

### Sensitivity envelopes about common-input points

```powershell
# Inventory
python -m Numerical.orchestration.studies.common_inputs.RunCommonInputEnvelopes

# Pilot numerical scans at a few gamma values
python -m Numerical.orchestration.studies.common_inputs.RunCommonInputEnvelopes `
  --run --limit 2 --points 7 --gamma-min 0.5 --gamma-max 2

# Full numerical envelope scan
python -m Numerical.orchestration.studies.common_inputs.RunCommonInputEnvelopes --run

# Rebuild HTML only, using already saved points
python -m Numerical.orchestration.studies.common_inputs.RunCommonInputEnvelopes --display-only
```

Requires `RunCommonInputRGE` results. Options: `--run`, `--display-only`, `--limit` (0 all), `--force`, `--gamma-min` (0.1), `--gamma-max` (10), `--points` (25), `--spacing log|linear` (log). `--display-only` **does not integrate**. Saved scan output: `Reports/output/full/common_input_sensitivity/`; dashboard `Reports/output/display/common_input_model_envelopes_all.html`.

### Common-input Yukawa and lambda coupling scan

```powershell
python -m Numerical.orchestration.studies.common_inputs.RunCommonInputCouplingScan

# Separate one-parameter deformation families around common points
python -m Numerical.orchestration.studies.common_inputs.RunCommonInputCouplingScan `
  --run --limit 2 --modes yukawa lambda --factors 0.5 1.0 2.0

# Include all common-input models, instead of selected subset
python -m Numerical.orchestration.studies.common_inputs.RunCommonInputCouplingScan `
  --run --all-models
```

Options: `--scenario`, `--run`, `--limit` (0 all), `--all-models`, `--factors` (default `0.5 0.75 1 1.5 2`), `--modes` (`yukawa` and/or `lambda`), `--force`. `yukawa` multiplies **both `y1` and `y2`**, keeping lambda fixed; `lambda` scales only `lambdaT3`. Includes a frozen-parameter control that **retains one-loop matching**; not a true standalone tree/fixed-order result. Output `Reports/output/full/common_input_coupling_scan/<scenario>/`.

### Common-input fermion mass scan

```powershell
python -m Numerical.orchestration.studies.common_inputs.RunCommonInputFermionScale

python -m Numerical.orchestration.studies.common_inputs.RunCommonInputFermionScale `
  --run --limit 2 --factors 0.5 1 2

python -m Numerical.orchestration.studies.common_inputs.RunCommonInputFermionScale --run --all-models
```

Options: `--scenario`, `--run`, `--limit` (default 0, selected subset), `--all-models`, `--factors` (default `0.5 0.75 1 1.5 2`; automatically includes 1), `--force`. Scales **MF and the matching scale together**, not the scalar threshold. Writes `Reports/output/full/common_input_fermion_scale/<scenario>/`.

### Cross-representation common input evaluation (16×16)

```powershell
# v1 historical/archived candidate module: inventory by default
python -m Numerical.orchestration.archive_candidates.RunCommonParameterComparison

# v1 pilot / all candidates
python -m Numerical.orchestration.archive_candidates.RunCommonParameterComparison `
  --run --limit-candidates 2

# v2 representation-aware cross-evaluation: inventory
python -m Numerical.orchestration.studies.common_inputs.RunCommonParameterComparisonV2

# v2 pilot / full 256-pair evaluation
python -m Numerical.orchestration.studies.common_inputs.RunCommonParameterComparisonV2 `
  --run --limit-candidates 2
python -m Numerical.orchestration.studies.common_inputs.RunCommonParameterComparisonV2 --run
```

Both expose `--run` and `--limit-candidates` (default 0 = all when running), but use different saved outputs and conventions. **v2** saves under `Reports/output/full/common_parameter_comparison_v2/smallY_smallL/`, including candidate rankings and best common input; v1 is an archived candidate branch, **not** the authoritative fitted-reference full-study runner.

### Optimization of shared Yukawa inputs

```powershell
# Default inventory/no fitting
python -m Numerical.orchestration.studies.common_inputs.OptimizeCommonYukawas

# SPSA search, saved per-iteration checkpoint
python -m Numerical.orchestration.studies.common_inputs.OptimizeCommonYukawas `
  --run --iterations 60 --seed 20261013

# Resume prior common-Yukawa optimization checkpoint
python -m Numerical.orchestration.studies.common_inputs.OptimizeCommonYukawas --run --resume
```

Options: `--run`, `--iterations` (60), `--seed` (20261013), `--resume`. Bounded shared-Yukawa SPSA optimization across available 16 representations. Uses results of `RunCommonParameterComparisonV2`, does not alter the 16 **frozen fitted-reference configurations**. Outputs `Reports/output/full/common_yukawa_optimization/smallY_smallL/`.

### Four benchmark correlations and two-lambda scans

```powershell
# LL / HL / LH / HH fixed benchmark comparisons (no optimization)
python -m Numerical.orchestration.studies.common_inputs.RunFourBenchmarkCorrelations
python -m Numerical.orchestration.studies.common_inputs.RunFourBenchmarkCorrelations `
  --run --limit-models 2
python -m Numerical.orchestration.studies.common_inputs.RunFourBenchmarkCorrelations --display-only

# Two lambdaT3 regimes, gamma sensitivity and cache-aware resumption
python -m Numerical.orchestration.studies.common_inputs.RunScalarCouplingCorrelations
python -m Numerical.orchestration.studies.common_inputs.RunScalarCouplingCorrelations `
  --run-missing --limit-models 2 --points 7 `
  --gamma-min 0.9 --gamma-max 1.1
```

`RunFourBenchmarkCorrelations`: `--run`, `--display-only`, `--retry-failed` (with `--run`), `--limit-models` (0 all). Reads the optimized common input at `Reports/output/full/common_yukawa_least_squares/<scenario>/optimized_common_input.json`; if absent, first produce it by the appropriate upstream workflow. Output `Reports/output/full/four_benchmark_correlations/<scenario>/four_benchmark_predictions.csv` and `interactive_comparison/index.html`. Benchmark factors: `LL=(1,1)`, `HL=(100,1)`, `LH=(1,100)`, `HH=(100,100)`.

`RunScalarCouplingCorrelations`: `--run-missing`, `--points` (7), `--gamma-min` (0.9), `--gamma-max` (1.1), `--limit-models` (0 all). Uses two absolute `lambdaT3` references, **0.0005** and **0.05**, for sensitivity scans using saved optimized common Yukawas; caches per point. Output `Reports/output/full/scalar_coupling_correlations/<scenario>/`. Default invocation inventories, `--run-missing` computes only uncached points.


## 9. Common-scalar and beta-function studies

These are **numerical studies using saved symbolic artifacts**, not a request to regenerate the complete Lagrangian or RGBeta expressions. The default mode of most study scripts is inventory/dry-run; `--run` enables computation.

### Common scalar coupling across representations

```powershell
# Inspect eligible saved fitted candidates
python -m Numerical.orchestration.studies.scalar_rge.RunCommonScalarComparison

# Pilot: two models and reduced Yukawa-target grid
python -m Numerical.orchestration.studies.scalar_rge.RunCommonScalarComparison `
  --run --limit 2 --target-yukawas 0.001 0.01 --max-yukawa 1.0

# All eligible models, default target Yukawa values
python -m Numerical.orchestration.studies.scalar_rge.RunCommonScalarComparison --run
```

**Flags:** `--scenario` (default `smallY_smallL`), `--run`, `--limit` (0 = all), `--reference-model`, `--target-yukawas` (default `1e-4 1e-3 1e-2 1e-1 1.0`), `--max-yukawa` (1.0), `--force` (re-evaluate completed results). The lambda grid is calibrated on **one reference representation** and then used for the others; Yukawas are rescaled to preserve the leading radiative mass combination. Equal scalar input does not imply equal Yukawa singular values between models.

**Outputs:** `Reports/output/full/common_scalar_comparison/<scenario>/`, including model-level `common_scalar_scan.json`, overall `comparison_summary.csv` and `comparison_summary.json`, and plots.

### Common-scalar running-versus-frozen controls

**Prerequisite:** successful `RunCommonScalarComparison` outputs.

```powershell
# Input availability only
python -m Numerical.orchestration.studies.scalar_rge.RunCommonScalarRGControls --dry-run

# Pilot then complete
python -m Numerical.orchestration.studies.scalar_rge.RunCommonScalarRGControls --run --limit 2
python -m Numerical.orchestration.studies.scalar_rge.RunCommonScalarRGControls --run
```

**Flags:** `--scenario`, `--run`, `--dry-run`, `--limit` (0 = all), `--force`. **Output:** `Reports/output/full/common_scalar_rg_controls/<scenario>/running_vs_fixed.json`, `running_vs_fixed.csv` and `plots/`. The frozen control **retains one-loop neutrino-mass matching** while freezing running parameters; it is not a full fixed-order calculation.

### UV beta-sector dominance after common-scalar scans

**Prerequisite:** common-scalar scan results and underlying saved RGBeta expressions.

```powershell
python -m Numerical.orchestration.studies.scalar_rge.RunCommonScalarBetaDominance --dry-run
python -m Numerical.orchestration.studies.scalar_rge.RunCommonScalarBetaDominance `
  --run --limit 2 --points 9
python -m Numerical.orchestration.studies.scalar_rge.RunCommonScalarBetaDominance --run
```

**Flags:** `--scenario`, `--source-root` (default `Reports/output/full/common_scalar_comparison`), `--output-root` (default `Reports/output/full/common_scalar_beta_dominance`), `--points` (9), `--limit` (0 = all), `--run`, `--dry-run`. This computes the gauge/scalar/Yukawa contribution to **UV Yukawa beta functions** by algebraic ablations of the exported beta expressions. It does not calculate independent full-theory physics corrections.

### Heavy-fermion mass scan with common-scalar benchmarks (Study B)

**Prerequisite:** results from `RunCommonScalarComparison`.

```powershell
python -m Numerical.orchestration.studies.fermion_scales.RunHeavyFermionScaleComparison
python -m Numerical.orchestration.studies.fermion_scales.RunHeavyFermionScaleComparison `
  --run --limit 2 --mass-factors 0.5 1.0 2.0 --reference-yukawas 0.01 1.0
python -m Numerical.orchestration.studies.fermion_scales.RunHeavyFermionScaleComparison --run
```

**Flags:** `--scenario` (`smallY_smallL`), `--run`, `--limit` (0 all), `--mass-factors` (default `0.5 1.0 2.0`), `--reference-yukawas` (default `0.01 1.0`), `--force`. **Output:** `Reports/output/full/heavy_fermion_scale_comparison/<scenario>/running_vs_fixed_mass_scan.json`, matching CSV, and plots. The scan scales the physical fermion mass **and the fermion matching scale together**, with the scalar threshold unchanged. Its frozen-parameter control is a diagnostic, not an independent UV completion.

---

## 10. Numerical diagnostics, audits, extraction

### Resume old fixed-scenario numerical results

```powershell
# Read existing status and report missing cases; does not recompute numerics
python -m Numerical.orchestration.support.ResumePipelineNumericalResults

# Resume up to two missing cases from precomputed artifacts
python -m Numerical.orchestration.support.ResumePipelineNumericalResults --run --limit 2

# Recompute already-completed cases only if intentionally desired
python -m Numerical.orchestration.support.ResumePipelineNumericalResults --run --force
```

**Flags:** `--run`, `--force`, `--limit` (at least 1). The default *inventory* still writes a resume summary JSON/CSV. Operates on previously saved **fixed comparison-grid** artifacts (`output/full/comparison`), **not** the new `pipeline.py --full` frozen-fitted run list. In `--run` mode it reuses saved matching/RGBeta and computes numerical continuation without symbolic rematching.

### Extract low-energy observables and verify saved continuation

```powershell
# Saved fixed-grid diagnostic JSON -> flat CSV/JSON (no ODE work)
python -m Numerical.diagnostics.ExtractFullComparisonResults
python -m Numerical.diagnostics.ExtractFullComparisonResults `
  --input-root "output\full\comparison" --output-root "Reports\output\full"

# Inspect saved Weinberg matching/continuation consistency
python -m Numerical.diagnostics.VerifyFullContinuation
python -m Numerical.diagnostics.VerifyFullContinuation `
  --input "output\full\comparison" `
  --output "Reports\output\full\continuation_verification"
```

**Extraction output:** `Reports/output/full/t3_full_comparison_observables.csv` and `.json`; `--input-root`, `--output-root` specify alternative locations. **Continuation audit output:** `Reports/output/full/continuation_verification/continuation_verification.csv` and `.json`; flags `--input`, `--output`. Both use stored results; they do not refit, rematch or rerun the ODE.

### Fast structural configuration validation

```powershell
python -m Numerical.diagnostics.NumericalConfigCheck
```

**No flags.** Checks dimensions/structural suitability of the 16 ordinary-T3 numerical configurations and saved optimal configurations where present. Does **not** run Matchete, RGBeta, matching, ODEs, or plots.

### UV one-loop beta-term and trajectory dominance

```powershell
# Algebraic UV beta-term decomposition at saved benchmark states
python -m Numerical.diagnostics.BetaTermDominance
python -m Numerical.diagnostics.BetaTermDominance `
  --scenario smallY_smallL --strict

# UV-scale trajectory decomposition: REINTEGRATES the UV ODE only
python -m Numerical.diagnostics.BetaTrajectoryDominance `
  --models representative --scenario smallY_smallL --points 33
python -m Numerical.diagnostics.BetaTrajectoryDominance `
  --models all --points 33
```

`BetaTermDominance` options: `--configs` (default `configs/generated_model_sets/comparison`), `--raw` (default `output/full/comparison`), `--output` (default `Reports/output/full/beta_dominance`), repeatable `--scenario` (one of the four old benchmark scenarios), repeatable `--model` (full `T3_dS1_...` key), `--strict`. It evaluates saved RGBeta UV beta terms but does **not** integrate RGEs.

`BetaTrajectoryDominance` accepts the same `--configs`, `--raw`, `--model`, `--scenario` and `--strict`, plus `--models representative|all` (default `representative`), `--points` (33), `--output` (default `Reports/output/full/beta_dominance/uv_trajectory`). **This one reintegrates only the UV numerical ODE using saved RGBeta expressions.** It does not perform new Matchete matching or run the full pipeline.

### SU(2) component, charge and multiplicity audit workflow

These independent scripts operate on exported Matchete CG/generator data. Run only when the prerequisites are present; **charge-allowed channels are not automatically physical loop multiplicities**.

```powershell
# 1. Inspect Matchete-exported component support
python -m Numerical.diagnostics.AuditT3ComponentMultiplicity `
  --root "output\full\hypercharge"

# 2. Build charge spectra and charge-conserving channel list
python -m Numerical.diagnostics.AuditT3ChargeSpectrum
python -m Numerical.diagnostics.AuditT3ChargeChannels

# 3. Explore contractions (no physical multiplicity assignment)
python -m Numerical.diagnostics.ContractT3LoopComponents

# 4. Generator analysis and SU(2) covariance verification
#    (requires previously exported matchete_generator_probe.json)
python -m Numerical.diagnostics.AnalyseMatcheteGenerators
python -m Numerical.diagnostics.VerifyT3CGCovariance

# 5. Resolve charge-basis weights only after all four inputs are available
python -m Numerical.diagnostics.ResolveT3ChargeWeights
```

| Entry point | Principal CLI options | Reads / writes (defaults) | Important scope |
|---|---|---|---|
| `AuditT3ComponentMultiplicity` | `--root`, `--output` | Reads `output/full/hypercharge`; writes `Reports/output/full/multiplicity/component_support_audit.json` | CG support only, no multiplicity correction |
| `AuditT3ChargeSpectrum` | `--input`, `--output` | Component audit → `charge_spectrum_audit.json` | SU(2) weight/electric charge spectra; unresolved triplet-basis warning |
| `AuditT3ChargeChannels` | `--input`, `--output` | Charge spectra → `charge_channels.json` | Tests `Q_F=Q_S1=Q_S2` for neutral external legs |
| `ContractT3LoopComponents` | `--input`, `--output` | Component support → exploratory contraction JSON | Does not certify physical multiplicity |
| `AnalyseMatcheteGenerators` | `--input`, `--output` | `matchete_generator_probe.json` → `matchete_generator_analysis.json` | Generator ordering/Cartan spectrum; not alone CG-basis certification |
| `VerifyT3CGCovariance` | None | Component support and generator probe → `cg_covariance_audit.json` | Checks infinitesimal SU(2) covariance; exits nonzero on failed verification |
| `ResolveT3ChargeWeights` | `--directory` | Four prerequisites → `charge_resolved_loop_weights.json` | Uses verified CG data to derive charge-resolved weights |

Except where explicit paths are supplied, these files are located under `Reports/output/full/multiplicity/`. A successful charge/channel inventory alone **must not** be used to alter physical `lambdaT3` or assign a numerical enhancement factor.

### Diagnostic layers that are **not** CLI programs

`Numerical/diagnostics/IntermediateWilsonDiagnostics.py`, `IntermediateWeinbergDiagnostics.py` and `FinalC5ContributionDiagnostics.py` expose imported calculation routines, **not standalone `python -m ...` entry points**. They are used downstream by the pipeline or other diagnostic/report scripts. The runnable `RunningDiagnostics` and `BestFitDiagnostics` commands are documented in [§6](#6-parameter-scans-fitting-benchmark-diagnostics).

---

## 11. Group theory, RGE and matching research CLIs

The commands in this section are **specialist physics-development/validation entry points**, rather than ordinary plotting/report generation. Many require **previously exported Matchete/RGBeta JSON**, and some write intermediate scientific artifacts. `PATH_TO_*.json` and `PATH_TO_*.txt` below are **placeholders**, not repository files claimed to exist at those names. Check the source/help before applying them to a different EFT or scalar representation.

### Analytic group-factor validators

```powershell
# Direct scalar-EFT C12 -> Weinberg group factor vs analytic SU(2) 6j
python -m RGE.group_factors.validation.ValidateDirectWeinbergGroupFactors
python -m RGE.group_factors.validation.ValidateDirectWeinbergGroupFactors `
  --root "output\full\hypercharge" `
  --output "output\group_factors\direct_weinberg_group_factor_validation.json"

# Heavy-fermion mass and scalar-mass beta coefficients: separate modes
python -m RGE.group_factors.validation.ValidateMassGroupFactors fermion
python -m RGE.group_factors.validation.ValidateMassGroupFactors scalar

# Scalar mixing-quartic beta(lambdaT3) factors
python -m RGE.group_factors.validation.ValidateMixingQuarticGroupFactors

# Yukawa beta group factors for a chosen (dS1,dS2,dF,alpha)
python -m RGE.group_factors.validation.YukawaBetaGroupFactors factors 2 2 1 -1

# Validate beta coefficients against saved RGBeta outputs
python -m RGE.group_factors.validation.YukawaBetaGroupFactors validate
```

| Module suffix | Arguments | Output / meaning |
|---|---|---|
| `ValidateDirectWeinbergGroupFactors` | `--root`, optional `--input` (legacy diagnostic JSON), `--output`, `--diagnostic-output` | Default `output/group_factors/direct_weinberg_group_factor_validation.json`; SU(2) direct-Weinberg comparison |
| `ValidateMassGroupFactors fermion` | `--root`, `--json-output` | Default `output/group_factors/fermion_mass_validation.json` |
| `ValidateMassGroupFactors scalar` | `--root`, `--json-output` | Default `output/group_factors/scalar_mass_validation.json` |
| `ValidateMixingQuarticGroupFactors` | `--root`, `--json-output` | Default `output/group_factors/lambdaT3_analytic_validation.json`; raw RGBeta non-singlet terms are convention-sensitive diagnostics |
| `YukawaBetaGroupFactors factors` | Required `dS1 dS2 dF alpha` | Prints representation-level Yukawa beta factors |
| `YukawaBetaGroupFactors validate` | `--root`, `--input`, `--output`, `--gauge-output`, `--diagnostic-output` | Default `output/group_factors/y2_nongauge_branch_validation.json` and optional companion outputs |

**Physics cautions from current code:** The mass and Yukawa validators explicitly preserve **known, unresolved self-conjugate branch/convention issues**. A validator running successfully does not erase these restrictions; inspect its JSON status and excluded cases. No scalar or Yukawa beta-function coefficient should be automatically "corrected" from a failed diagnostic without an independent derivation.

### RGBeta and scalar recoupling research modes

```powershell
# Five separate RGBeta tensor/convention investigations
python -m RGE.group_factors.validation.ValidateRGBetaRecouplings diagnose-bc-pseudoreal
python -m RGE.group_factors.validation.ValidateRGBetaRecouplings bc-ward
python -m RGE.group_factors.validation.ValidateRGBetaRecouplings bc-tensor
python -m RGE.group_factors.validation.ValidateRGBetaRecouplings e-invariant
python -m RGE.group_factors.validation.ValidateRGBetaRecouplings e-cross

# Analytic portal recoupling derivations (optional --output PATH)
python -m RGE.group_factors.recoupling.PortalRecouplingDerivations e-cross
python -m RGE.group_factors.recoupling.PortalRecouplingDerivations adjoint
python -m RGE.group_factors.recoupling.PortalRecouplingDerivations gauge

# Portal-group factors for a representation and validation against saved output
python -m RGE.group_factors.recoupling.PortalQuarticGroupFactors factors 2 2 1 -1
python -m RGE.group_factors.recoupling.PortalQuarticGroupFactors validate
```

`ValidateRGBetaRecouplings` requires one of the five **positional mode** values shown above; some modes also accept optional `--output PATH` (check each mode's `--help`). `PortalRecouplingDerivations` modes output, respectively, `output/group_factors/e_portal_cross_scalar_recouplings.json`, `portal_adjoint_recouplings.json`, and `portal_gauge_quartics.json` by default. `PortalQuarticGroupFactors validate` accepts `--root` (default `output/full/hypercharge`), `--output` (default `output/group_factors/portal_quartic_consolidation.json`), optional `--backbone-output`; `factors` requires four integer representation parameters.

### Scalar quartic basis and direct SU(2) tensor checks

```powershell
# Quartic-basis map, with optional --rgbeta PATH and --output PATH
python -m RGE.group_factors.core.ScalarQuarticBasis `
  basis-map "PATH_TO_QUARTIC_SEED.json"

# Invariant-completeness check for dS1=2,dS2=2
python -m RGE.group_factors.core.ScalarQuarticBasis `
  invariant-completeness "PATH_TO_QUARTIC_SEED.json" 2 2

# Yukawa scalar-wavefunction SU(2) factor from exported UV CG registry
python -m RGE.group_factors.core.YukawaGroupFactors `
  "PATH_TO_UV_CG_REGISTRY.json" --yukawa y1
python -m RGE.group_factors.core.YukawaGroupFactors `
  "PATH_TO_UV_CG_REGISTRY.json" --yukawa y2

# Independent SU(2) mixing-quartic recoupling
python -m RGE.group_factors.recoupling.RGBetaMixingQuarticRecoupling 2 2 `
  --rgbeta-json "PATH_TO_UV_RGBETA.json" --output "output\group_factors\mixing_quartic_check.json"

# Real-basis mixing-quartic tensor recoupling from saved scalar seed
python -m RGE.group_factors.recoupling.MixingQuarticRecoupling `
  "PATH_TO_QUARTIC_SEED.json" --rgbeta "PATH_TO_EFT1_RGBETA.json"
```

`ScalarQuarticBasis` takes a **required positional mode** `basis-map|invariant-completeness`; the completeness mode additionally needs `dS1 dS2`. `YukawaGroupFactors` requires an exported CG JSON and takes `--yukawa y1|y2` (default y1). `RGBetaMixingQuarticRecoupling` takes positional `dS1 dS2` plus `--rgbeta-json`, `--output`. `MixingQuarticRecoupling` takes a positional scalar quartic seed plus optional `--rgbeta` (normally inferred alongside seed) and `--output`. These are checks of mathematical conventions/tensors, not general-purpose model generation commands.

### Hierarchical final Weinberg coefficient and intermediate EFT research

```powershell
# Assemble final C5 from already exported threshold and flavor transport files
python -m RGE.matching.FinalWeinbergCoefficient `
  "PATH_TO_C5_COEFFICIENT.txt" "PATH_TO_FLAVOR_TRANSPORT.json" `
  --output "output\research\final_weinberg_coefficient.json"

# Direct C12 -> Weinberg contribution between fermion and scalar thresholds
python -m RGE.running.intermediate.DirectWeinbergRunning `
  "PATH_TO_FLAVOR_SEED.json" "PATH_TO_COMPONENT_RGE.json" `
  --mu-high MF --mu-low MS --output "output\research\direct_weinberg.json"

# Full scalar-only Wilson tensor RGE: can be costly
python -m RGE.running.intermediate.ScalarOnlyWilsonTensorRGE `
  "PATH_TO_WILSON_SEED.json" "PATH_TO_QUARTIC_SEED.json" `
  "PATH_TO_RGBETA.json" --output "output\research\wilson_rge.json"

# Fast matched-component plumbing check instead of full tensor closure
python -m RGE.running.intermediate.ScalarOnlyWilsonTensorRGE `
  "PATH_TO_WILSON_SEED.json" "PATH_TO_QUARTIC_SEED.json" `
  "PATH_TO_RGBETA.json" --seed-support-only
```

- `FinalWeinbergCoefficient`: two required positional paths; optional `--output` (default `final_weinberg_coefficient.json`) and `--pole-rge-consistency PATH` (otherwise auto-discovered beside the threshold data). This **combines existing** hard and direct contributions; it is not a standalone Matchete run.
- `DirectWeinbergRunning`: positional flavor seed and component RGE JSON, `--mu-high` (default `MF`), `--mu-low` (`MS`), `--output`. Exports the direct one-loop scalar-EFT Wilson contribution.
- `ScalarOnlyWilsonTensorRGE`: required three positional JSON paths, optional `--output`, `--seed-support-only`. By default studies allowed tensor **closure**, not only nonzero boundary entries.

### Scalar-only Wilson transport and basis adapters

```powershell
# Extract the flavor seed used by scalar-only EFT transport
python -m RGE.running.intermediate.ScalarOnlyWilsonFlow `
  flavor-seed "PATH_TO_WILSON_SEED.json" --output "output\research\flavor_seed.json"

# Fixed-order component Wilson transport between specified thresholds
python -m RGE.running.intermediate.ScalarOnlyWilsonFlow `
  transport "PATH_TO_WILSON_SEED.json" "PATH_TO_WILSON_RGE.json" `
  "PATH_TO_RGBETA.json" --mu-high MF --mu-low MS `
  --output "output\research\component_wilson_transport.json"

# Real-field scalar quartic tensor diagnostics from saved seeds
python -m RGE.running.intermediate.ScalarOnlyTensorAdapters `
  "PATH_TO_QUARTIC_SEED.json" "PATH_TO_RGBETA.json"
```

`ScalarOnlyWilsonFlow` uses a **required first positional mode** `transport|flavor-seed`; `transport` requires all three JSON files and `--mu-high/--mu-low`, while `flavor-seed` requires the Wilson seed and optionally `--output`. `ScalarOnlyTensorAdapters`' active `__main__` executes the **quartic tensor** diagnostics with two positional JSONs; its other tensor-construction routines (e.g. Wilson chirality controls) are imported APIs, not additional simultaneous subcommands.

### Explicit lower scalar threshold rerun — external Wolfram dependency

```powershell
# Expert-only: use existing EFT continuation to rerun lower scalar matching
python -m RGE.running.intermediate.ScalarThresholdMatching `
  "PATH_TO_MODEL_OUTPUT_DIR" "PATH_TO_RUNNING_INSERTION.json" `
  --dims 2 2 1 --alpha -1 --eft-order 5 --loop-order 1
```

Requires **positional** `output_dir running_insertion` and **required** `--dims`. Other flags: `--alpha`, `--script` (default `Lagrangian/RunThresholdStage.wl`), `--eft-order` (5), `--loop-order` (1), `--validation` (expensive consistency rematching), `--tree`, `--loop`, `--transition-tree`, `--transition-full`, `--cg-registry`. This can call Wolfram and **write matching outputs**. It is a specialist threshold calculation, not needed to rebuild saved reports.

### Modules intended for import rather than `python -m` usage

Examples: `Lagrangian/Runner.py`, `Lagrangian/T3ModelMatching.py`, `RGE/running/UVRunning.py`, `RGE/running/IntermediateEFTRunning.py`, `RGE/matching/FinalWeinbergAdapter.py`, `RGE/matching/WeinbergFlavorMatching.py`, `Numerical/core/State.py`, `Numerical/running/T3Trajectory.py`, `Numerical/running/FinalC5TrajectoryAdapter.py`, `Numerical/orchestration/PipelineNumericalResults.py`, `Reports/ReportGeneration.py`, `Reports/GroupFactorReports.py` and `Reports/RGEComparison.py` are **implementation components** without standalone main commands. Use their calling pipeline/study entry points above.

---

## 12. Tests and regression

### Python-only contract/unit tests

```powershell
# All Python tests (may take time, but doesn't launch --full as a command)
python -m pytest -q

# Pipeline CLI, execution-plan and full-study contracts
python -m pytest -q tests/test_pipeline_cli.py `
  tests/test_pipeline_plan.py tests/test_full_study_architecture.py

# Presentation and display regeneration
python -m pytest -q tests/test_build_display_results.py `
  tests/test_interactive_model_comparison.py `
  tests/test_within_model_reports.py

# Fitting and scan interfaces
python -m pytest -q tests/test_oscillation_fit.py `
  tests/test_oscillation_optimizer.py tests/test_scan_cli.py

# Numeric status/resume and intermediate EFT boundaries
python -m pytest -q tests/test_resume_pipeline_numerical_results.py `
  tests/test_intermediate_scalar_boundary.py

# Majorana normalization, Weinberg coefficient and trajectory checks
python -m pytest -q tests/python/test_t3_scotogenic_normalization.py `
  tests/test_final_c5_bridge.py tests/test_full_flavor_weinberg.py
```

The individual tests exist on the audited commit. Python-level test success does not independently establish a fresh Wolfram/Matchete/RGBeta full end-to-end run or global phenomenological validity.

### External smoke and Wolfram C5 regression

```powershell
# Regenerate five smoke models; runs external symbolic dependencies
python pipeline.py --smoke

# Reuse existing output/smoke and run Wolfram regression on C5
python scripts/run_regression.py --no-rematch

# Full regression: re-run smoke before checking C5 with wolframscript
python scripts/run_regression.py
```

`run_regression.py` accepts `--no-rematch` and `--verbose`. **Do not use its `--verbose` flag with a rematch on this commit**: the wrapper passes `--verbose` through to `pipeline.py`, but the current pipeline parser in `common/PipelineCLI.py` **does not define `--verbose`**. The Wolfram regression uses `tests/wolfram/RegressionC5.wl` and writes `output/smoke/t3_regression_report.json` when its expected output exists. `--no-rematch` requires suitable existing smoke artifacts and still runs the Wolfram validation.

### Source-level and runtime-validation distinction

The command syntax and meanings in this reference have been **verified against files on GitHub**. The external tools and your local Python dependencies have **not been run** here; where they are required, the commands are examples for your configured machine. A failed run may indicate missing artifacts or prerequisites rather than a failed physics calculation.

---

## 13. Input/output layout and safe workflow order

### Repository locations that matter

| Local project path | Purpose | Handling |
|---|---|---|
| `pipeline.py` | Master symbolic → matching → RGE → numerical → report orchestrator | Runs full matching; expensive |
| `common/PipelineCLI.py` | Authoritative top-level pipeline flags | Source parser, not separate command |
| `studies/FullT3Study.py` | Current 16-model frozen fitted-reference `--full` loop | Does **not** rebuild old 64 fixed scenarios |
| `Numerical/benchmarks/fitted_reference/benchmark_manifest.json` | Manifest, representation entries, frozen input hashes | **Preserve** for reproducibility |
| `Numerical/benchmarks/fitted_reference/*.json` | 16 frozen independently fitted inputs | Do not overwrite with a proposed common fit |
| `configs/nufit_6_1_ic24_no.json` | NuFIT 6.1 normal-ordering oscillation target | Check dataset convention before thesis claims |
| `configs/generated_model_sets/comparison/<scenario>/input/` | Scenario-specific comparison input configurations | Do not silently replace with resolved/fit configs |
| `configs/generated_model_sets/optimal/` | Older optimal numerical model inputs | Distinguish from frozen fitted reference manifest |
| `output/full/` | Saved underlying symbolic/RGE and numerical raw outputs | **Do not delete when only cleaning plots** |
| `output/full/comparison/<scenario>/.../running_diagnostics.json` | Old four-scenario grid saved trajectories | Input to saved-grid plots, audits and resume tools |
| `output/full/fitted_reference/` | New `--full` results for frozen-reference 16-model study | Current full-study output |
| `output/full/hypercharge/` | Saved symbolic exports used in group-factor checks | Prerequisite for some validators |
| `Reports/output/full/` | Study scan JSON/CSV and some per-model reports | Scientific result records; save provenance |
| `Reports/output/analytical/` | Combined TeX/PDF and equation report outputs | Source for analytical report polish/compact tools |
| `Reports/output/display/` | Rebuildable PNG/HTML displays | Safe **output destination** for presentations |
| `Reports/output/display_stage5/` | Default umbrella `tools.RebuildReports` destination | Safe report output; separate from raw data |
| `Reports/figure_review/` | Figure index and optional thesis copy exports | Selection of existing graphics, not computation |
| `Reports/report_cleanup_archive/` | Archived analytical build trash/backups | Managed by analytical housekeeper |
| `output/group_factors/` | Group-factor validation and tensor-derivation JSON | Treat failing/pending results as diagnostics |
| `tests/`, `tests/python/`, `tests/wolfram/` | Python regression and external symbolic validation | Distinguish testing mechanisms |

**Four historical scenarios:** `smallY_smallL`, `smallY_largeL`, `largeY_smallL`, `largeY_largeL`. They represent nominal fixed-grid Yukawa/scalar input scales. Their names are retained in some later fitted/scan paths, but a later fitted candidate need **not** remain in the nominal small/large regime unless its constraints explicitly enforce it.

### Recommended safe order for common tasks

#### A. I only want the figures/reports from results I already have

```powershell
# 1. Rebuild saved old-grid static + interactive comparisons
python -m Numerical.plotting.BuildDisplayResults --output "Reports\output\display"

# 2. Add combined interactive correlation HTML if saved fitted sensitivity exists
python -m Numerical.plotting.BuildInteractiveCorrelations `
  --source "Reports\output\full\fitted_sensitivity" `
  --output "Reports\output\display\interactive_correlations.html"

# 3. Optionally rebuild saved staged-RG and Weinberg figures
python tools/RebuildSavedStudyPlots.py --study staged-rg
python tools/RebuildSavedStudyPlots.py --study weinberg

# 4. Optionally produce a figure index
python -m Numerical.orchestration.archive_candidates.CollectThesisFigures
```

All four commands **reuse saved outputs**. They require the corresponding input files to exist; do not launch full symbolic calculation merely because a downstream scan is missing. For a combined fixed-grid dashboard/correlation build use the shorter `python -m Numerical.plotting.BuildCombinedInteractiveReports` after confirming fitted sensitivity source exists.

#### B. I want fitted sensitivity / RGE-stage plots but the numerical scan is missing

```powershell
# Inventory/locate saved fitted candidates and artifacts
python -m Numerical.orchestration.studies.fitted_rge.RunFittedSensitivityComparison
python -m Numerical.orchestration.studies.fermion_scales.RunFittedFermionScaleScan

# Pilot fresh numerical-only computations
python -m Numerical.orchestration.studies.fitted_rge.RunFittedSensitivityComparison `
  --run --limit 2 --points 7
python -m Numerical.orchestration.studies.fitted_rge.RunFittedStagedRGComparison `
  --run --limit 2 --factors 1.0

# Replot data afterward, without recomputing anything
python -m Numerical.plotting.BuildFittedModelEnvelopes --scenario smallY_smallL
python tools/RebuildSavedStudyPlots.py --study staged-rg
```

**Note:** Fitted sensitivity scripts may expect pre-existing local `fitted_candidate_config.json` in earlier fit outputs. Frozen reference input JSONs do not automatically provide every missing historical study artifact. When a study's inventory reports missing data, use its error to identify the prerequisite rather than constructing fake cases.

#### C. I want a new physics run for one T3 model

```powershell
# Inspect supported representation and threshold constraints
python pipeline.py --help

# Example ordinary model with fermion-first hierarchy
python pipeline.py --dims 2 2 1 --alpha -1 `
  --study "single_B_alpha_minus1" `
  --numerical "configs\t3_numerical_default_B_alpha_m1.json" `
  --threshold F --threshold S1 S2
```

This is **full symbolic matching plus numerical evaluation**, using the current production fermion-first threshold hierarchy. Check external Wolfram/Matchete/RGBeta and LaTeX environment before running.

#### D. I want the current 16-model fitted-reference study

```powershell
# Reproducibility check before a costly run
python -m Numerical.diagnostics.NumericalConfigCheck
python -m pytest -q tests/test_full_study_architecture.py

# Deliberate full symbolic + numerical 16-model execution
python pipeline.py --full
```

Input: the **manifest-controlled frozen reference** in `Numerical/benchmarks/fitted_reference/`. Output: `output/full/fitted_reference/`, `Reports/output/full/fitted_reference/`. The full runner verifies hashes and resets derived model configs to the reference inputs. It is not the former 64-case comparison sweep and is not a new fit.

#### E. I want to investigate common parameters, not independently fitted parameters

```powershell
# First common scalar benchmark scan (pilot)
python -m Numerical.orchestration.studies.scalar_rge.RunCommonScalarComparison --run --limit 2

# Separate running/frozen numerical controls
python -m Numerical.orchestration.studies.scalar_rge.RunCommonScalarRGControls --run --limit 2

# Or first compare common INPUTS across all representations
python -m Numerical.orchestration.studies.common_inputs.RunCommonInputRGE --run --limit 2

# Follow with common-input sensitivity envelopes if baseline exists
python -m Numerical.orchestration.studies.common_inputs.RunCommonInputEnvelopes `
  --run --limit 2 --points 7 --gamma-min 0.5 --gamma-max 2.0
```

The common-input comparison and common-scalar comparison are **different scans**; neither should be labeled an independent oscillation fit for all 16 models.

### PowerShell command-output logging

For an important run, preserve the exact command and output next to your results. Example:

```powershell
# Output log; `Tee-Object` writes the console transcript to file
python -m Numerical.orchestration.studies.fitted_rge.RunFittedStagedRGComparison `
  --run --limit 2 --factors 1.0 2>&1 | Tee-Object -FilePath "staged_rg_console.log"
```

For reproducibility also record `git rev-parse HEAD`, original input config path and hash, and the output scan manifest/JSON. `--force` and `--reset-numerical-configs` can change the data you intend to compare.

---

## 14. Legacy paths and known pitfalls

The following are **actual source/documentation differences** found in the audited GitHub snapshot. They are recorded to prevent confusing help/ImportError messages or results being generated with the wrong assumptions.

| Obsolete / misleading instruction or assumption | Current verified interpretation or replacement |
|---|---|
| `python -m Numerical.orchestration.RunFittedSensitivityComparison` appears in `README.md` and `docs/NUMERICAL.md` | Module is now **`python -m Numerical.orchestration.studies.fitted_rge.RunFittedSensitivityComparison`** |
| `python -m Numerical.orchestration.BuildComparisonReports` is called in `scripts/FullRunWithWithinModelReports.ps1` | Actual file is **`Numerical/orchestration/archive_candidates/BuildComparisonReports.py`**, CLI **`python -m Numerical.orchestration.archive_candidates.BuildComparisonReports`**; the `.ps1` wrapper should not be assumed to work unchanged |
| Old `RunCommonParameterComparison` command inside the v2 script's docstring | Distinct **v1** `Numerical.orchestration.archive_candidates.RunCommonParameterComparison` and **v2** `Numerical.orchestration.studies.common_inputs.RunCommonParameterComparisonV2` exist; choose intentionally |
| `python pipeline.py --full` described as smoke + hypercharge + dimension sweeps in `common/PipelineCLI.py` help | **Actual `studies/FullT3Study.py` implementation uses 16 frozen fitted reference models**; `--full` does not mean a fresh fixed-input 64-case grid |
| `python pipeline.py --full --numerical config.json` | **Rejected**: the full runner owns manifest-specified model configs; use `--dims ... --numerical ...` for a custom individual run |
| `python scripts/run_regression.py --verbose` with rematch | Regression forwards `--verbose` to `pipeline.py`, which **has no `--verbose` option** in its current parser; avoid it. Plain regression or `--no-rematch` is the meaningful command |
| `scripts/FullStudyShortNames.py` claims to override `_model_key` for short full-run study folders | `studies/FullT3Study.py` now enumerates manifest-defined model names; short-name monkey-patch is legacy; do not rely on its earlier effect |
| `BuildDisplayResults --report-input` implies reading the reports input | Source explicitly says `--report-input` is **retained for compatibility but ignored**. `--data-input` is the primary source |
| `Numerical/plotting/BuildFittedModelEnvelopes` produces new gamma data | It only displays **previously saved** sensitivity diagnostics; use a fitted/common-input scan to calculate missing gamma points |
| `--display-only` and `--run` on study CLI are interchangeable | They are different. Display-only uses prior CSV/JSON; `--run` computes missing numerical points; default often inventories |
| `--force` is required every time | No: most numerical study implementations cache completed outputs; `--force` deliberately overwrites/recomputes and can lose comparability |
| `--reset-numerical-configs` is just a plotting reset | **No.** It overwrites model-specific generated numerical configuration files from the selected template |
| A successful RGE trajectory is a good neutrino-oscillation fit | **No.** Numerical integration and fitted chi-square are separate. The frozen fitted configurations were **calibrated**; a small fitted residual is not independent predictive viability |
| Four scenario labels guarantee four disjoint strict coupling ranges during later fits | **No.** Nominal fixed-grid settings need not be enforced in later feasibility optimizations |
| `gamma` filter in fitted sensitivity means absolute Yukawa magnitude | **No.** `gamma` is a prescribed deformation factor; check each T1/T2/T3 sensitivity direction |
| Direct `LLSS -> O5` scalar-only EFT contribution is the whole matched one-loop neutrino mass | **No.** Distinguish direct running, hard matching at scalar threshold and `C5_final`; staged-RG on/off is not a separate loop-order matching calculation |
| Source-level `pytest` equals external symbolic regression | **No.** Python tests check contracts; use the separate `pipeline.py --smoke` / Wolfram regression for external checks |
| SU(2) charge-allowed channels imply a physical multiplicity enhancement | **No.** CG-basis conventions and contractions require independent validated group-theory treatment |
| LaTeX report scripts can safely consume missing Matchete exports | **No.** Combined analytical PDF generation requires saved model records/source expressions, and PDF compilation requires a TeX engine |

**Historical note:** `docs/RESULTS.md` records legacy fixed-grid coverage (including the 64-case numerical comparison); it is a summary of saved results, **not** the current `--full` run definition. Similarly, archived scenario `input/`, `resolved/`, and `numerical_resume/` configs are not automatically interchangeable.

---

## 15. File and command index

This index identifies **runnable CLIs discovered in source**, grouped by purpose. Detailed flags, examples and input/output folders are in the preceding sections. A file existing in the repository does **not** imply it should be run directly: ordinary library modules are deliberately excluded.

| Category | Entrypoint files / modules | Section |
|---|---|---|
| Main symbolic physics pipeline | `pipeline.py`; `scripts/run_comparison_study.py` | [3](#3-master-pipeline-symbolic-matching-and-full-runs) |
| Full-run wrappers and regression | `scripts/FullStudyOrganised.py`; `scripts/FullStudyShortNames.py`; `scripts/FullRunWithWithinModelReports.ps1` (**stale import**); `scripts/run_regression.py` | [3](#3-master-pipeline-symbolic-matching-and-full-runs), [12](#12-tests-and-regression), [14](#14-legacy-paths-and-known-pitfalls) |
| Saved display/rebuild | `Numerical/plotting/BuildDisplayResults.py`; `BuildCombinedInteractiveReports.py`; `InteractiveModelComparison.py`; `BuildFittedModelEnvelopes.py`; `BuildInteractiveCorrelations.py`; `BuildStudyABInteractive.py` | [4](#4-rebuild-reports-and-dashboards-from-saved-data) |
| Saved scientific figures | `Numerical/plotting/PlotScanResults.py`; `PlotOptimizerResults.py`; `PlotFittedSensitivityFigures.py`; `NormalizedFittedSensitivity.py`; `ObservableSensitivityDashboard.py`; `RunningResultFigures.py`; `ThesisResultFigures.py`; `LLSSWilsonRunningReport.py` | [4](#4-rebuild-reports-and-dashboards-from-saved-data), [6](#6-parameter-scans-fitting-benchmark-diagnostics) |
| Display and collection tools | `tools/RebuildReports.py`; `tools/RebuildSavedStudyPlots.py`; `Numerical/orchestration/archive_candidates/BuildComparisonReports.py`; `CollectThesisFigures.py`; `Numerical/orchestration/presentation/BuildThreePanelInteractive.py` | [4](#4-rebuild-reports-and-dashboards-from-saved-data) |
| Analytical TeX/PDF | `scripts/BuildCombinedAnalyticalReports.py`; `BuildCompactC5Report.py`; `AddMathematicalDefinitions.py`; `PolishAnalyticalReports.py` | [5](#5-analytical-latex-and-pdf-reports) |
| Presentation file changes | `scripts/RenameT3ModelFolders.py`; `scripts/FullStudyOrganised.py` | [5](#5-analytical-latex-and-pdf-reports) |
| Scanning and optimization | `Numerical/fitting/ScanCLI.py`; `OscillationOptimizer.py`; `Numerical/diagnostics/BestFitDiagnostics.py`; `RunningDiagnostics.py` | [6](#6-parameter-scans-fitting-benchmark-diagnostics) |
| Fitted RGE | `Numerical/orchestration/studies/fitted_rge/RunFittedRGEDisplay.py`; `RunFittedSensitivityComparison.py`; `RunFittedStagedRGComparison.py`; `RunFittedWeinbergRunningComparison.py` | [7](#7-fitted-rge-threshold-and-mass-scale-studies) |
| Fitted fermion mass | `Numerical/orchestration/studies/fermion_scales/RunFittedFermionScaleScan.py` | [7](#7-fitted-rge-threshold-and-mass-scale-studies) |
| Common-input studies | `Numerical/orchestration/studies/common_inputs/RunCommonInputRGE.py`; `RunCommonInputEnvelopes.py`; `RunCommonInputFermionScale.py`; `RunCommonInputCouplingScan.py`; `RunCommonParameterComparisonV2.py`; `OptimizeCommonYukawas.py`; `RunFourBenchmarkCorrelations.py`; `RunScalarCouplingCorrelations.py` | [8](#8-common-input-and-common-yukawa-studies) |
| Archived common-input cross-evaluation | `Numerical/orchestration/archive_candidates/RunCommonParameterComparison.py`; `Numerical/orchestration/RunYukawaGrouping.py` | [8](#8-common-input-and-common-yukawa-studies) |
| Common-scalar + related | `Numerical/orchestration/studies/scalar_rge/RunCommonScalarComparison.py`; `RunCommonScalarRGControls.py`; `RunCommonScalarBetaDominance.py`; `Numerical/orchestration/studies/fermion_scales/RunHeavyFermionScaleComparison.py` | [9](#9-common-scalar-and-beta-function-studies) |
| Numerical extraction, resume | `Numerical/orchestration/support/ResumePipelineNumericalResults.py`; `Numerical/diagnostics/ExtractFullComparisonResults.py`; `NumericalConfigCheck.py`; `VerifyFullContinuation.py`; `BetaTermDominance.py`; `BetaTrajectoryDominance.py` | [10](#10-numerical-diagnostics-audits-extraction) |
| Matchete-CG and charge audits | `Numerical/diagnostics/AuditT3ComponentMultiplicity.py`; `AuditT3ChargeSpectrum.py`; `AuditT3ChargeChannels.py`; `ContractT3LoopComponents.py`; `AnalyseMatcheteGenerators.py`; `VerifyT3CGCovariance.py`; `ResolveT3ChargeWeights.py` | [10](#10-numerical-diagnostics-audits-extraction) |
| Group-factor validators | `RGE/group_factors/validation/ValidateDirectWeinbergGroupFactors.py`; `ValidateMassGroupFactors.py`; `ValidateMixingQuarticGroupFactors.py`; `YukawaBetaGroupFactors.py`; `ValidateRGBetaRecouplings.py` | [11](#11-group-theory-rge-and-matching-research-clis) |
| SU(2) tensor / recoupling research | `RGE/group_factors/core/ScalarQuarticBasis.py`; `YukawaGroupFactors.py`; `RGE/group_factors/recoupling/PortalRecouplingDerivations.py`; `PortalQuarticGroupFactors.py`; `RGBetaMixingQuarticRecoupling.py`; `MixingQuarticRecoupling.py` | [11](#11-group-theory-rge-and-matching-research-clis) |
| Wilson / Weinberg matching research | `RGE/matching/FinalWeinbergCoefficient.py`; `RGE/running/intermediate/DirectWeinbergRunning.py`; `ScalarOnlyWilsonFlow.py`; `ScalarOnlyWilsonTensorRGE.py`; `ScalarOnlyTensorAdapters.py`; `ScalarThresholdMatching.py` | [11](#11-group-theory-rge-and-matching-research-clis) |
| Tests | `python -m pytest -q`; `python pipeline.py --smoke`; `python scripts/run_regression.py` | [12](#12-tests-and-regression) |

**Scope statement:** All top-level modes, actual runnable report/plotting/study scripts, notable numerical diagnostics, and specialist RGE/group-factor entry points found during the source audit are documented here. The project also contains a large collection of **internal Python implementation classes, Wolfram scripts invoked indirectly, saved JSON data, and archived/generated reports**; these are not all independent shell commands. This reference should be revisited when `main` changes.

**Provenance / implementation versus verification:** Flags and module locations were checked against the `main` commit named at the top. Numerical predictions, Matchete exports, external tool availability, and user-local configuration validity were **not** tested by this documentation task. The repository was read, **not edited**.
