$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$repo = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
Push-Location $repo
try {
    Write-Host 'T3_C_V21: self-test SU2/U1/Z2 two-trivalent topology closure'
    python scripts/CloseT3BCOneLoopTopologies.experimental.py --self-test
    if ($LASTEXITCODE -ne 0) { throw 'v21 self-test failed' }
    Write-Host 'T3_C_V21: source-anchor archived B/C UV and prior F/Y/X matched EFT evidence'
    python scripts/CloseT3BCOneLoopTopologies.experimental.py --output output/audits/t3_bc_topology_closure_v21/summary.json
    if ($LASTEXITCODE -ne 0) { throw 'v21 source verification failed' }
    Write-Host 'T3_C_V21: topology branch closed in restricted archived B/C UV; full physical C5 NOT CERTIFIED'
} finally {
    Pop-Location
}
