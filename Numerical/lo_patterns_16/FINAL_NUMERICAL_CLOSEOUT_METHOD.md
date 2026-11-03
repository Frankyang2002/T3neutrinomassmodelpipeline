# Final T3 numerical closeout — one bounded batch

## Purpose and stopping point

This is the **last batch for rough, source-conditional T3 neutrino-mass numerics**. It consolidates *existing* output and closes the numerical milestone if all validation checks pass. It deliberately **does not** attempt B1/B2 physical normalisation, new RGBeta/Matchete SU(2) work, RGE, EFT or threshold matching.

The physics input for every quoted mass is the saved candidate left-Weyl construction, not a signed and phase-certified physical Matchete/RGBeta action dictionary. The frozen results remain **16/16 provisionally evaluated, 0/16 physically source-certified**.

## GitHub evidence inspected

At preparation, the latest GitHub `main` was `2ac58cc3de12c7af11e2e7e214415dadf0138f4f`. The inputs are the tracked frozen outputs under `output/general_observable/lo_pattern_study16`:

- `provisional_01/all16_summary.json`, `all16_summary.csv`, `scans.csv`, `HANDOFF.json`
- `analysis_batch05_01/batch05_synthesis.json`, `SOURCE_MANIFEST.json`
- `analysis_batch06_01/batch06_factorization.json`, `batch06_coefficient_tests.csv`, `SOURCE_MANIFEST.json`

The existing source results support all 16 model masses, 288 scan points, 64 counterfactual controlled-texture cases (16 × 4), and five Batch 6 coefficient scenarios (16 × 5). Their output files already record SHA256 links from the frozen study to Batch 5 and Batch 6. The new script checks those links against actual local file content and the cross-source model, benchmark and candidate-input hashes.

The original frozen study was captured at an earlier Git HEAD. The closeout **does not** enforce the old HEAD or mistake a current HEAD for the frozen capture. It preserves the original revision as historical metadata. Hash validation tolerates CRLF-to-LF changes *only* when the normalized file is byte-identical to the recorded source. It cannot substitute for a fresh native physics-source audit.

## Output

Generated only on successful local validation, in the new directory `output/general_observable/lo_pattern_study16/final_numerical_closeout_01`:

1. `FINAL_NUMERICAL_CLOSEOUT.md`: source-conditional numerical findings, limitations, and explicit stop rule.
2. `all16_provisional_results.csv`: three provisional Takagi masses in eV, hierarchy ratios, Batch 5 smallest-mode geometry/coherence labels, Batch 6 common-mixing coefficient residual.
3. `controlled_texture_comparisons_64.csv`: fixed common-|lambdaT3|=0.10 Yukawa-texture substitutions; **not refitted physical benchmarks**.
4. `coefficient_scenario_summary_5.csv`: maximum coefficient residual and heavy-group spread for the five Batch 6 scenarios. `charge_split_diagnostic` is explicitly outside the restricted SU(2)-degenerate scalar ansatz.

No JSON tracker, no figures, no amplitude recalculation and no original file modification.

## Acceptance checks

- Freeze counts exactly: 16 evaluated, 0 blocked, 0 physical certificates, 288 scans passed.
- Cross-check original 16-model JSON and CSV masses; compare Batch 5 derived masses and hierarchy ratios. Positive ordered singular masses required.
- Exact scan grid: `y1`, `y2`, `lambdaT3` at `(-1,0,0.5,1,2)`, `heavy_mass` at `(0.5,1,2)` for every model; all 48 zero-coupling cases vanish; all 64 `t=1` points agree with baseline. Positive Yukawa rescalings reproduce linear mass scaling.
- Batch 5: 16 records × four common-mixing reference textures, with separate inputs and model-lineage SHA256 consistency.
- Batch 6: 16 × five coefficient comparisons, cross-checked against the CSV; all 80 saved tensor-reconstruction residuals below `1e-8`; baseline coefficient-fit residual below `1e-8`.
- Original source-conditional flags must remain false; no physical certificate can be issued by this batch.

## PowerShell installation and execution

Run from `C:\thesistemp\T3neutrinomassmodelpipeline`. The ZIP is a single archive of **complete new files**; it does not overwrite physics modules.

```powershell
Set-Location "C:\thesistemp\T3neutrinomassmodelpipeline"
$zip = Join-Path $env:USERPROFILE "Downloads\T3_Final_Provisional_Numerical_Closeout_16.zip"
if (!(Test-Path -LiteralPath $zip)) { throw "ZIP missing: $zip" }
Expand-Archive -LiteralPath $zip -DestinationPath "." -Force

.\.venv\Scripts\python.exe -m pytest tests\test_final_numerical_closeout_16.py -q
if ($LASTEXITCODE -ne 0) { throw "Final closeout unit tests failed" }

$out = "output\general_observable\lo_pattern_study16\final_numerical_closeout_01"
if (Test-Path -LiteralPath $out) { throw "Closeout output already exists: $out" }
.\.venv\Scripts\python.exe scripts\CloseOutProvisionalT3Numerics16.py --output $out
if ($LASTEXITCODE -ne 0) { throw "Numerical closeout source checks failed" }

Get-Content (Join-Path $out "FINAL_NUMERICAL_CLOSEOUT.md") -TotalCount 18
```

If that command succeeds, the **rough numerical milestone is complete**. The source-conditional results can be discussed as such in the thesis. Their masses must **not** be presented as externally certified physical predictions.
