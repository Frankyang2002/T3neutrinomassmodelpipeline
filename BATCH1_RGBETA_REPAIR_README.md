# Batch 1 RGBeta native-command repair (v2)

## Root cause

Batch 1 used raw negative integers such as `-4` as positional arguments to
`wolframscript -file`. On the user's Windows WolframScript installation, this
led the RGBeta runner to see a wrong number of `$ScriptCommandLine` entries,
printing the usage error rather than constructing beta functions.

The repository's `T3RGBetaParseIntegerToken` already accepts `m4` for `-4`,
`m3` for `-3`, `m2` for `-2`, and `m1` for `-1`. This release uses that notation
for all negative-alpha WolframScript invocations. It **does not change** the
physics Lagrangian or native RGBeta scripts.

## Full replacement files / paths

- `scripts/RunT3Batch1GeneralContract.py` — complete runner replacement
- `tests/test_t3_batch1_contract.py` — complete test-module replacement
- `BATCH1_RGBETA_REPAIR_README.md` — new instructions

Everything else stays as it was (including the corrected
`Lagrangian/interactions/ScalarPotential.wl`, the original Batch 1 tensor-map
module, and previously computed native outputs).

## PowerShell, run from the repository root

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline

# Keep local copies of the two overwritten files if you changed them yourself.
Copy-Item scripts\RunT3Batch1GeneralContract.py scripts\RunT3Batch1GeneralContract.py.before_v2 -Force
Copy-Item tests\test_t3_batch1_contract.py tests\test_t3_batch1_contract.py.before_v2 -Force

Expand-Archive -LiteralPath "$env:USERPROFILE\Downloads\T3_Batch1_RGBeta_Repair_v2.zip" -DestinationPath . -Force

.\.venv\Scripts\python.exe -X utf8 -m unittest discover -s tests -p 'test_t3_batch1_contract.py' -v

.\.venv\Scripts\python.exe -X utf8 -m scripts.RunT3Batch1GeneralContract `
  --run-id batch1_native_v1 `
  --repair-rgbeta `
  --a4-source 'output\general_observable\fermion_first_hermitian_portal_v1\A_m4'
```

`--repair-rgbeta` requires the existing `summary.json` and `provenance.json`
and compares SHA-256 of the original physics source files before reusing
existing outputs. It re-runs only missing or invalid UV/EFT1 RGBeta JSONs.
It never re-runs `f_threshold`, never changes the model Lagrangians, and
preserves any pre-repair log, invalid JSON, or report as a
`*.before_batch1_repair_N` file. By default, all 16 model reports are updated,
including those whose betas were reused without modification.

### Quick status check after repair

```powershell
$base = 'output\batch1_general_contract\batch1_native_v1'
Get-ChildItem "$base\models\*.json" | ForEach-Object {
    $m = Get-Content $_.FullName -Raw | ConvertFrom-Json
    [PSCustomObject]@{
        Model=$m.model
        UV=$m.rgbeta.uv.status
        EFT1=$m.rgbeta.eft1.status
        UVReady=$m.rgbeta.uv.candidate_ready
        EFT1Ready=$m.rgbeta.eft1.candidate_ready
    }
} | Sort-Object Model | Format-Table -AutoSize
```

Expected **if** the native jobs now complete: all 16 UV and all 16 EFT1
records say `SuccessCandidateMassBeta`. They will *still* carry outstanding
normalization, one-loop anomalous-dimension, and hard-matching conditions;
these are not native betas failing.

## Verification boundary

11/11 Python tests passed in the artifact environment, including invocation
argument tests for every model and a mocked repair-reuse test. WolframScript
itself cannot be run in this environment. Native outputs must be checked on
Windows. This repair neither creates physical C5 nor certifies any unresolved
basis dictionary. GitHub remains untouched.
