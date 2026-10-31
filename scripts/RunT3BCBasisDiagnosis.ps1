$ErrorActionPreference = 'Stop'
$root = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
Push-Location $root
try {
    python -u -X utf8 -m scripts.RunT3BCBasisDiagnosis @args
    if ($LASTEXITCODE -ne 0) { throw "B/C basis diagnosis exited with $LASTEXITCODE" }
} finally { Pop-Location }
