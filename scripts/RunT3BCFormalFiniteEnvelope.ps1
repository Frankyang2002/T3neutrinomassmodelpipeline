$ErrorActionPreference = 'Stop'
Set-Location (Resolve-Path (Join-Path $PSScriptRoot '..'))
$py = 'scripts\AssembleT3BCFormalFiniteScaleEnvelope.experimental.py'
Write-Host 'T3_C_V16: self-test'
& python $py --self-test
if ($LASTEXITCODE -ne 0) { throw 'v16 self-test failed' }
Write-Host 'T3_C_V16: B/C formal common-scale envelope'
& python $py --output 'output\audits\t3_bc_formal_finite_envelope_v16\summary.json'
if ($LASTEXITCODE -ne 0) { throw 'v16 B/C envelope failed' }
Write-Host 'T3_C_V16: done; physical C5 NOT CERTIFIED'
