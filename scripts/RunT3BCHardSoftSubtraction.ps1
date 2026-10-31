$ErrorActionPreference = 'Stop'
Set-Location (Resolve-Path (Join-Path $PSScriptRoot '..'))
$script = 'scripts\VerifyT3BCHardSoftSubtraction.experimental.py'
$dest = 'output\audits\t3_bc_hard_soft_subtraction_v19\summary.json'
Write-Host 'T3_C_V19: source-independent analytic MSbar hard/soft subtraction self-test'
python $script --self-test
if ($LASTEXITCODE -ne 0) { throw 'v19 self-test failed' }
Write-Host 'T3_C_V19: source-anchored F/Y/X pole assignment and frozen 3-flavour matrices'
python $script --output $dest
if ($LASTEXITCODE -ne 0) { throw 'v19 B/C source comparison failed' }
Write-Host "T3_C_V19: finished; $dest"
Write-Host 'T3_C_V19: restricted triangle subtraction only; FULL physical C5 NOT CERTIFIED'
