# Project path: scripts/RunT3QuarticPreflight.ps1
# Read-only preflight. Run from repository root.
$ErrorActionPreference = 'Stop'
python -X utf8 -m scripts.RunGeneralT3ScalarTensorCertification --static-only
if ($LASTEXITCODE -ne 0) { throw 'Stage 22 static audit failed' }
python -X utf8 -m scripts.AuditT3SimultaneousQuarticReadiness
if ($LASTEXITCODE -ne 0) { throw 'Quartic inventory audit failed' }
python -X utf8 -m scripts.RunT3SimultaneousEFT1All16 --self-test
if ($LASTEXITCODE -ne 0) { throw 'Ordinary simultaneous self-test failed' }
python -X utf8 -m scripts.RunT3SimultaneousEFT1MixedBC --self-test
if ($LASTEXITCODE -ne 0) { throw 'Mixed simultaneous self-test failed' }
Write-Host 'Existing solver tests passed; full quartic tensor remains uncertified until native normalization evidence and implementation are available.'
