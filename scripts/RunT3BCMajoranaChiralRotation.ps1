$ErrorActionPreference = 'Stop'
$root = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
Push-Location $root
try {
    python -u -X utf8 -m scripts.RunT3BCMajoranaChiralRotation
    if ($LASTEXITCODE -ne 0) { throw "Native probe failed (exit $LASTEXITCODE); inspect output/studies/t3_bc_majorana_chiral_rotation_v1/native_matchete.log" }
} finally { Pop-Location }
