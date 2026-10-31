$ErrorActionPreference = 'Stop'
Set-Location -LiteralPath (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
Write-Host 'B/C mixed quadratic: independent mass-insertion/propagator and real-native-source tests'
& python -X utf8 -m unittest discover -s tests -p 'test_t3_bc_mixed_quadratic_kernel_v2.py' -v
if ($LASTEXITCODE -ne 0) { throw 'mixed quadratic regression failed' }
& python -u -X utf8 -m scripts.RunT3BCMixedQuadraticKernelAudit
if ($LASTEXITCODE -ne 0) { throw 'B/C mixed quadratic audit failed; inspect saved summary.json' }
Write-Host 'Read output\studies\t3_bc_mixed_quadratic_kernel_v2\summary.json; physical C5 remains uncertified'
