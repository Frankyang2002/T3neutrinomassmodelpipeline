param([switch]$NoPreviousReference)
$ErrorActionPreference = 'Stop'
Set-Location (Split-Path -Parent $PSScriptRoot)
python -m unittest tests.test_t3_bc_f_hard_numerical_batch_v1 -v
if ($LASTEXITCODE -ne 0) { throw "F-hard numeric tests failed (exit $LASTEXITCODE)" }
$argsRun = @('-u','-X','utf8','-m','scripts.RunT3BCFHardNumericalBatch')
if ($NoPreviousReference) { $argsRun += '--no-previous-reference' }
python @argsRun
if ($LASTEXITCODE -ne 0) { throw "F-hard numeric batch failed (exit $LASTEXITCODE)" }
