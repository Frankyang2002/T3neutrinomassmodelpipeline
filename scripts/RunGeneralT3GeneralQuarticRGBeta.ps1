param(
    [switch]$ExportBetas,
    [switch]$CheckGauge,
    [string[]]$Models = @(),
    [ValidateSet('UV','EFT1','both')][string]$Stage = 'both'
)
$ErrorActionPreference = 'Stop'
$root = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
Push-Location $root
try {
    $arguments = @('-m', 'scripts.RunGeneralT3GeneralQuarticRGBeta', '--stage', $Stage)
    if ($ExportBetas) { $arguments += '--export-betas' }
    if ($CheckGauge) { $arguments += '--check-gauge' }
    foreach ($model in $Models) { $arguments += @('--model', $model) }
    & python @arguments
    if ($LASTEXITCODE -ne 0) { throw "Stage 20 failed with exit code $LASTEXITCODE" }
}
finally { Pop-Location }
