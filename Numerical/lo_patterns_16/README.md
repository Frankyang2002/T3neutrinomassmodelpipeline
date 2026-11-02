# Local installation and execution (PowerShell)

Target local repository: `C:\thesistemp\T3neutrinomassmodelpipeline`.
Extract the accompanying ZIP into that repository **without replacing any production file**. The package contains only new `Numerical/lo_patterns_16/`, `scripts/RunLOPatternStudy16.py`, and `tests/test_lo_pattern_study_16.py`.

```powershell
Set-Location 'C:\thesistemp\T3neutrinomassmodelpipeline'
git rev-parse HEAD
git status --short
.\.venv\Scripts\python.exe -m pytest -q tests\test_lo_pattern_study_16.py
.\.venv\Scripts\python.exe scripts\RunLOPatternStudy16.py
```

The default reads `output\general_observable\lo_direct_one_loop\provisional_inputs\A_-4.json`, etc., plus all frozen `Numerical\benchmarks\fitted_reference` benchmark inputs, and writes **new** `output\general_observable\lo_pattern_study16\provisional_01`. It fails if that output directory already contains files. Use `--output-dir output\general_observable\lo_pattern_study16\provisional_02` for a second run.

Optional custom input provider:

```powershell
.\.venv\Scripts\python.exe scripts\RunLOPatternStudy16.py `
  --input-dir '.\output\general_observable\lo_direct_one_loop\provisional_inputs' `
  --output-dir '.\output\general_observable\lo_pattern_study16\provisional_02'
```

Run the **16 real frozen-input integration tests** (not synthetic tests) explicitly:

```powershell
$env:RUN_LO_PATTERN_INTEGRATION='1'
.\.venv\Scripts\python.exe -m pytest -q tests\test_lo_pattern_study_16.py
Remove-Item Env:RUN_LO_PATTERN_INTEGRATION
```

This integration test writes into and cleans up a uniquely named temporary subdirectory below the owned `lo_pattern_study16` output root. It never edits the frozen inputs or production source.

**Expected status if current local inputs and source contracts match:** 16 `conditional_uv_evaluated`, zero `blocked`, zero `physical_source_certified`. This is an **expectation**, not a test result from this sandbox. If a source file has drifted, an input is missing, a UV pole does not cancel, or independent calculations disagree, the corresponding model is `blocked` and the command exits with status `2`. Scan failures also return exit status `2`; see `scans.csv` for reasons. `--physical-dictionary` intentionally blocks all unverified documents until the coordinator and both physical dictionaries define an independently checkable certificate contract; no JSON boolean is accepted as proof.

The source-code module does not import EFT or threshold machinery, and it never modifies GitHub.
