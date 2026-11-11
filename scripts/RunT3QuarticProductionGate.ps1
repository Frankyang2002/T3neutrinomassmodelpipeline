$ErrorActionPreference = 'Stop'
$repo = Split-Path -Parent (Split-Path -Parent $PSCommandPath)
Set-Location $repo
python -X utf8 -m scripts.TestT3QuarticProductionGate
if ($LASTEXITCODE -ne 0) { throw 'Production gate tests failed' }
python -X utf8 -m scripts.TestT3ExtendedSimultaneousEFT1
if ($LASTEXITCODE -ne 0) { throw 'Joint EFT1 tests failed' }
python -X utf8 -m scripts.RunT3SimultaneousEFT1All16 --self-test
if ($LASTEXITCODE -ne 0) { throw 'Ordinary solver regression failed' }
python -X utf8 -m scripts.RunT3SimultaneousEFT1MixedBC --self-test
if ($LASTEXITCODE -ne 0) { throw 'Mixed solver regression failed' }
Write-Host 'PASS: Batch 5 tests; production remains blocked without independent native normalization and response binding.'
