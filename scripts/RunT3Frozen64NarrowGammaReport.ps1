# Local-only workflow. Run from any directory; repository is the script's parent.
# Uses existing frozen-general scalar results. Does NOT push or modify GitHub.
[CmdletBinding()]
param(
    [string]$RunId = 'frozen_general64_v1',
    [string]$WeinbergOutputId = 'weinberg_running_v1',
    [string]$ScanId = 'gamma_095_105_v1'
)
$ErrorActionPreference = 'Stop'
$root = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
Push-Location $root
try {
    $manifest = Join-Path $root "output/studies/general_scalar_extended/full_scalar_frozen_64/$RunId/manifest.json"
    if (-not (Test-Path -LiteralPath $manifest -PathType Leaf)) {
        throw "Missing frozen-64 manifest: $manifest. Run the numerical campaign first."
    }
    $grid = '0.95,0.96,0.97,0.98,0.99,1.00,1.01,1.02,1.03,1.04,1.05'

    # Compute actual scale-dependent C5 and neutrino observables (if absent).
    python -u -X utf8 -m scripts.GenerateT3Frozen64WeinbergRunning `
        --run-id $RunId --output-id $WeinbergOutputId --resume
    if ($LASTEXITCODE -ne 0) { throw 'Weinberg-running generation failed.' }

    # Fresh, narrow 0.95–1.05 gamma calculations. Old scan files are NOT reused.
    python -u -X utf8 -m scripts.GenerateT3FrozenFullScalarGamma `
        --run-id $RunId --scan-id $ScanId --gamma-grid $grid --resume
    if ($LASTEXITCODE -ne 0) { throw 'Narrow-gamma sensitivity run failed.' }

    # Uses the user's actual optimised_rge.html drawing code and styling.
    python -X utf8 -m Numerical.plotting.BuildT3FullScalarFrozen64HTML `
        --run-id $RunId --scan-id $ScanId `
        --weinberg-output-id $WeinbergOutputId
    if ($LASTEXITCODE -ne 0) { throw 'RGE/sensitivity HTML generation failed.' }

    $report = Join-Path $root "Reports/output/display_general_scalar_frozen_64/$RunId/optimised_low/interactive_rge.html"
    if (-not (Test-Path -LiteralPath $report -PathType Leaf)) {
        throw "The expected RGE report was not generated: $report"
    }
    Write-Host "Reference-style report: $report"
    Start-Process -FilePath $report
}
finally {
    Pop-Location
}
