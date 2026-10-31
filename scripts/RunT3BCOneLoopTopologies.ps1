$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $PSScriptRoot
Set-Location $root
Write-Host 'T3_C_V20: synthetic graph and source-census tests'
python 'scripts\ClassifyT3BCOneLoopTopologies.experimental.py' --self-test
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
Write-Host 'T3_C_V20: audit current local UV source and saved F/Y/X one-loop matched EFTs'
python 'scripts\ClassifyT3BCOneLoopTopologies.experimental.py' `
  --output 'output\audits\t3_bc_one_loop_topology_completeness_v20\summary.json'
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
Write-Host 'T3_C_V20: finished; full physical C5 remains uncertified'
