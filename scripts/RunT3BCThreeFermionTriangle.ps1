$ErrorActionPreference = 'Stop'
$here = (Resolve-Path '.').Path
$script = 'scripts\EvaluateT3BCThreeFermionTriangle.experimental.py'
$output = 'output\audits\t3_bc_three_fermion_triangle_v18\summary.json'
if (-not (Test-Path $script)) { throw "Missing $script; extract v18 ZIP at project root" }
Write-Host 'T3_C_V18: self-test three-heavy-flavour frozen-reference topology'
python $script --self-test
if ($LASTEXITCODE -ne 0) { throw 'v18 self-test failed' }
Write-Host 'T3_C_V18: evaluate B/C current scalar spectrum and fitted-reference frozen Yukawas'
python $script --root "$here" --output "$output"
if ($LASTEXITCODE -ne 0) { throw 'v18 B/C numerical reference failed' }
Write-Host "T3_C_V18: completed; $output"
Write-Host 'T3_C_V18: physical C5 and measured neutrino masses NOT CERTIFIED'
