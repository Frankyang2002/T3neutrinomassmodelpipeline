param(
    [string]$Root = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
)
$ErrorActionPreference = 'Stop'
Push-Location $Root
try {
    python .\scripts\RunGeneralT3Stage17.py
    if ($LASTEXITCODE -ne 0) { throw "Stage 17 failed with exit code $LASTEXITCODE; inspect output/general_t3_stage17/stage17_summary.json" }
} finally { Pop-Location }
