# Run the existing full T3 study, then generate within-model comparisons.
# No changes to the physics pipeline or its CLI.
$ErrorActionPreference = 'Stop'
$root = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
Push-Location $root
try {
    & python pipeline.py --full @args
    $fullExit = $LASTEXITCODE
    if ($fullExit -ne 0) {
        Write-Error "Full study failed (exit code $fullExit); within-model reports were not generated."
        exit $fullExit
    }

    & python -m Numerical.orchestration.BuildComparisonReports --strict
    $reportsExit = $LASTEXITCODE
    if ($reportsExit -ne 0) {
        Write-Error "Full study succeeded, but within-model report generation failed (exit code $reportsExit)."
        exit $reportsExit
    }
    Write-Host 'Full study and within-model interactive reports completed.'
    exit 0
}
finally {
    Pop-Location
}
