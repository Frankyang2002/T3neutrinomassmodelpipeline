$ErrorActionPreference = "Stop"
Set-Location (Resolve-Path (Join-Path $PSScriptRoot ".."))
$script = "scripts\BuildT3BCOneLoopInheritance.experimental.py"
Write-Host "T3_C_V15: independent self-test"
python $script --self-test
if ($LASTEXITCODE -ne 0) { throw "v15 self-test failed" }
Write-Host "T3_C_V15: source-verified B/C O(hbar) Y->X inheritance and F/Y/X residues"
python $script --output "output\audits\t3_bc_one_loop_inheritance_v15\summary.json"
if ($LASTEXITCODE -ne 0) { throw "v15 batch failed" }
Write-Host "T3_C_V15: completed; full physical C5 remains uncertified"
