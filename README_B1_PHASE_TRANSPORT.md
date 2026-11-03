# B1 third batch — source-spurion phase transport across all 16 models

Repository root: `C:\thesistemp\T3neutrinomassmodelpipeline`. Extract **B1_Phase_Transport_Batch.zip** into this root. This ZIP contains **new** isolated B1 files only, no replacement source or existing audit runners. No new Wolfram run is required.

## Installation and test commands (PowerShell)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline

$zip = Get-ChildItem "$env:USERPROFILE\Downloads" -Filter 'B1_Phase_Transport_Batch*.zip' -File |
    Sort-Object LastWriteTime -Descending | Select-Object -First 1
if (-not $zip) { throw 'B1 phase ZIP not found in Downloads' }
Expand-Archive -LiteralPath $zip.FullName -DestinationPath (Get-Location).Path -ErrorAction Stop

# Tests include 30 new tests; prior local suite had 44.
.\.venv\Scripts\python.exe -m unittest discover -s tests -p 'test_physical_b1_*.py' -v
if ($LASTEXITCODE -ne 0) { throw 'B1 test failure' }

$base = 'output\general_observable\physical_normalisation'
$native = Join-Path $base 'b1_02\native'
$action = Join-Path $base 'b1_action_contract_20261103_103545\native_action_coefficients.json'
$prior = Join-Path $base 'b1_retry_20261103_102017\all16_component_check.json'
$cc = Join-Path $base 'b1_cc_probe_20261103_103545\A_p0.json'
foreach ($p in @($action,$prior,$cc)) {
    if (-not (Test-Path $p)) { throw "Missing previously executed B1 evidence: $p" }
}

$stamp=Get-Date -Format 'yyyyMMdd_HHmmss'
$out=Join-Path $base "b1_phase_local_$stamp"

.\.venv\Scripts\python.exe scripts\RunPhysicalB1PhaseTransport.py `
    --native-dir $native `
    --action-audit $action `
    --conditional-audit $prior `
    --cc-json $cc `
    --output-dir $out
if ($LASTEXITCODE -ne 0) { throw 'B1 phase transport blocked; report the terminal error' }

$r=Get-Content (Join-Path $out 'phase_transport_all16.json') -Raw | ConvertFrom-Json
$r | Select-Object status,source_flavour_benchmarks_checked_count,maximum_numerator_covariance_residual,physical_source_certified
$r.models | Select-Object model,fermion_type,U3_flavour_covariance_residual,symmetric_nu_pair_covariance_residual | Format-Table -AutoSize
```

Optionally supply the **already executed, separate B2 conditional Ward** report, if present: append `--b2-ward-json 'output\general_observable\physical_normalisation\b2_ward_local_01\ward_witness.json'` to the Python runner invocation. This adds only source-linked metadata; it does NOT turn B2's literal-index diagnostics into a physical source certificate. No B2 run is required for B1.

## Expected status — do not confuse this with certification

`All16SourceSpurionFlavourAndPhaseCovarianceVerified_PhysicalDictionaryPending`;
`source_flavour_benchmarks_checked_count=16`;
relative covariance residuals near double-precision rounding;
`physical_source_certified=false` and `physical_source_certified_models=0`.

This is an **expected** local result, not yet observed on Windows. All output is written to a fresh `b1_phase_local_YYYYMMDD_HHMMSS` folder; the tool refuses to overwrite an existing folder and exits 2 on missing/mismatched source, bench, frozen, or native action evidence. Nothing in `Numerical/direct_one_loop`, frozen benchmarks, Matchete model files, GitHub, or EFT is modified.

## Project file paths created

1. `Numerical/lo_certification/physical_normalisation/b1_phase_transport.py`
2. `Numerical/lo_certification/physical_normalisation/B1_PHASE_TRANSPORT_DERIVATION.md`
3. `scripts/RunPhysicalB1PhaseTransport.py`
4. `tests/test_physical_b1_phase_transport.py`
5. `output/general_observable/physical_normalisation/b1_phase_transport_batch_01/HANDOFF.md`
6. `output/general_observable/physical_normalisation/b1_phase_transport_batch_01/evidence.json`
7. `README_B1_PHASE_TRANSPORT.md`

The runtime additionally creates `output/general_observable/physical_normalisation/b1_phase_local_<timestamp>/phase_transport_all16.json` and `HANDOFF.md`.

## Interpretation and next gate

The chosen source monomials imply `y1*y2*conj(MF)*conj(lambdaT3)`; the complete complex flavour contraction is `y1 @ MF.conj().T @ y2.T * conj(lambdaT3)` in the explicitly declared mass index order. An isolated Yukawa sign flip changes this invariant if all other coefficients are held fixed. Coherent field rephasings preserve it, and external neutrino basis rotations produce congruence. These are **operator/basis checks**; source physical normalisation remains blocked by installed Matchete's concrete `GammaCC` sign and the independent B2 numerical SU(2)/RGBeta action dictionary. **0/16 physically certified.**
