$ErrorActionPreference = 'Stop'
$repo = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
Push-Location $repo
try {
    $python = Join-Path $repo '.venv\Scripts\python.exe'
    if (-not (Test-Path $python)) { $python = 'python' }
    & $python (Join-Path $repo 'scripts\RunGeneralT3Stage19.py')
    $code = $LASTEXITCODE
    if ($code -ne 0) {
        throw "Stage19 source projection failed with exit $code; inspect output/general_t3_stage19/stage19_summary.json"
    }
} finally { Pop-Location }
