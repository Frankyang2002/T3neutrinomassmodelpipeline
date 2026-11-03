# Chat A final analytical synthesis — install and verification

## Purpose

Consolidate Batches 1–8 into one **provisional, source-conditional** technical report. All previous numerical calculations must already exist in the local checkout. This stage performs **no additional UV loop calculations, fitted Yukawa scans, Matchete/RGBeta calls, or EFT steps**.

This ZIP installs four new project-relative files:

1. `Numerical/lo_patterns_16/FINAL_TECHNICAL_SYNTHESIS.md` — the standalone analytical narrative with all 16 CG and frozen-mass entries, equations and limitations. Its fixed figures are explicitly identified as reported from the preceding Windows runs.
2. `Numerical/lo_patterns_16/final_synthesis.py` — fail-closed cross-batch source evidence verifier and local report constructor.
3. `Numerical/lo_patterns_16/FINAL_INSTALL.md` — this file.
4. `scripts/BuildLOFinalSynthesis16.py` — one-command runner.
5. `tests/test_lo_final_synthesis_16.py` — 20 deterministic synthetic tests and one opt-in real-checkout test.

## PowerShell

```powershell
$repo = 'C:\thesistemp\T3neutrinomassmodelpipeline'
$zip = Join-Path $HOME 'Downloads\Chat_A_Batch09_Final_Analytical_Synthesis.zip'
Expand-Archive -LiteralPath $zip -DestinationPath $repo
Set-Location $repo

.\.venv\Scripts\python.exe -m pytest -q tests\test_lo_final_synthesis_16.py
.\.venv\Scripts\python.exe scripts\BuildLOFinalSynthesis16.py
```

The fresh output directory is

`output/general_observable/lo_pattern_study16/analysis_final_synthesis_01/`

It contains:

- `FINAL_VERIFIED_TECHNICAL_REPORT.md` — the narrative **plus** tables extracted and verified from saved results at full useful precision;
- `consolidated_evidence.json` — all 16 result rows, factorization, pair comparisons and scalar-splitting diagnostics;
- `all16_verified.csv`, `scenario_maxima.csv`, `matched_pairs.csv`;
- `SOURCE_MANIFEST.json`, `HANDOFF.md` and `HANDOFF.json`.

If `analysis_final_synthesis_01/` is already nonempty, the runner refuses to overwrite it. To create a new revision, choose a fresh **child directory** of the existing Chat A output tree:

```powershell
.\.venv\Scripts\python.exe scripts\BuildLOFinalSynthesis16.py --output-dir (Join-Path $repo 'output\general_observable\lo_pattern_study16\analysis_final_synthesis_02')
```

The verifier checks all 16 model rows, 288 saved Batch-2 scan rows, seven matched pairs × four common textures, all 80 Batch-6 model–scenario coefficients, 16 Batch-7 corrections, 16 Batch-8 convergence summaries, the saved SHA256 chain, the original candidate CG-weight audit and false physical certification flags.

**Important:** A successful verifier proves that the eight **saved, conditional** result sets are complete and internally consistent. It does not independently certify the physical Matchete/RGBeta spinor, coupling or CG normalisation. Physical status remains **0/16 certified**.

For an additional *opt-in* integration test against the actual checkout rather than synthetic fixtures:

```powershell
$env:RUN_LO_FINAL_SYNTHESIS_INTEGRATION='1'
.\.venv\Scripts\python.exe -m pytest -q tests\test_lo_final_synthesis_16.py
Remove-Item Env:RUN_LO_FINAL_SYNTHESIS_INTEGRATION
```

The default tests do not need any frozen project inputs. If a source-hash or schema error occurs, do **not** edit the old reports or disable the guards: supply the exact error and source file involved for reconciliation. Nothing in this package modifies GitHub or previous local batch outputs.
