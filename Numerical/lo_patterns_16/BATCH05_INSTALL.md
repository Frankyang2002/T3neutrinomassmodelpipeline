# Batch 5 — installation and verification

All six files are **new**, relative to the T3 repository root:

1. `Numerical/lo_patterns_16/texture_synthesis.py`
2. `Numerical/lo_patterns_16/BATCH05_THEORY.md`
3. `Numerical/lo_patterns_16/BATCH05_INSTALL.md`
4. `scripts/SynthesiseLOPatterns16.py`
5. `tests/test_lo_texture_synthesis_16.py`
6. `Numerical/lo_patterns_16/BATCH05_RESULTS_INTERPRETATION.md`

**Prerequisites:** Chat A Batches 1, 2, 3, 4 already installed and run; the existing `provisional_01`, `analysis_batch03_01`, `analysis_batch04_01` results exist and source files match their saved SHA256 hashes. No changes to the shared production direct-UV modules. `.venv` with `numpy`, `pytest`, and optionally `matplotlib`.

```powershell
$repo = 'C:\thesistemp\T3neutrinomassmodelpipeline'
$zip = Join-Path $HOME 'Downloads\Chat_A_Batch05_Yukawa_Causality_Synthesis.zip'
Expand-Archive -LiteralPath $zip -DestinationPath $repo
Set-Location $repo
.\.venv\Scripts\python.exe -m pytest -q tests\test_lo_texture_synthesis_16.py
.\.venv\Scripts\python.exe scripts\SynthesiseLOPatterns16.py
```

Expected unit-test result based on isolated synthetic validation: **26 passed, 1 skipped**. The user’s real frozen file set is not mounted here, so this is a prediction, not a claim about a completed local integration run. The integration test is intentionally opt-in and will write to its own fresh `analysis_batch05_pytest_integration_01` directory.

The analytical runner prints one line for each model and creates only:

`output/general_observable/lo_pattern_study16/analysis_batch05_01/`

If it reports “Output directory already populated,” do not use `-Force` or delete the earlier report. Instead:

```powershell
.\.venv\Scripts\python.exe scripts\SynthesiseLOPatterns16.py --output-dir output\general_observable\lo_pattern_study16\analysis_batch05_02
```

To run without matplotlib:

```powershell
.\.venv\Scripts\python.exe scripts\SynthesiseLOPatterns16.py --no-plots
```

To run the **actual** full-suite integration test deliberately:

```powershell
$env:RUN_LO_TEXTURE_INTEGRATION='1'
.\.venv\Scripts\python.exe -m pytest -q tests\test_lo_texture_synthesis_16.py
Remove-Item Env:\RUN_LO_TEXTURE_INTEGRATION
```

The opt-in test uses `analysis_batch05_pytest_integration_01` and cannot overwrite it. Run it only if you need this extra integration gate.

Primary files to share for interpretation: `BATCH05_ANALYTICAL_SYNTHESIS.md`, `batch05_synthesis.json`, and optionally the three PNGs. Physical certification stays **0/16**, regardless of analytical fit quality.
