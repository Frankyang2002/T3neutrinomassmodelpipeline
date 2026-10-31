$ErrorActionPreference = 'Stop'
Set-Location (Resolve-Path (Join-Path $PSScriptRoot '..'))
Write-Host 'T3_C_V17: independent triangle analytic and integration tests'
python scripts\CompareT3BCFiniteEnvelopeToTriangle.experimental.py --self-test
if ($LASTEXITCODE -ne 0) { throw 'v17 triangle master self-test failed' }
Write-Host 'T3_C_V17: B/C source-matched signed triangle-master decomposition'
python scripts\CompareT3BCFiniteEnvelopeToTriangle.experimental.py --output output\audits\t3_bc_triangle_master_v17\summary.json
if ($LASTEXITCODE -ne 0) { throw 'v17 source topology check failed; physical C5 remains uncertified' }
Write-Host 'T3_C_V17: complete; physical C5 remains uncertified'
