param(
    [int]$Samples = 17
)

$ErrorActionPreference = 'Stop'
$repo = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
Push-Location $repo
try {
    $python = Join-Path $repo '.venv\Scripts\python.exe'
    if (-not (Test-Path $python)) { $python = 'python' }
    & $python (Join-Path $repo 'scripts\RunGeneralT3EFT1WilsonSourceTransport.py') --samples $Samples
    $code = $LASTEXITCODE
    if ($code -ne 0) { throw "Stage18 extra-Wilson transport failed with exit $code; inspect output/general_t3_stage18/stage18_summary.json" }
} finally {
    Pop-Location
}
