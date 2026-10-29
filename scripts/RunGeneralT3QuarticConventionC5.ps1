<# Stage 15: source-checked scalar quartic convention and B/C(-1) UV sensitivity. #>
param(
    [string]$InputConfig = 'docs\general_c5_toy_benchmark.json',
    [string]$Output = 'output\general_t3_stage15',
    [switch]$ConventionOnly
)
$ErrorActionPreference = 'Stop'
$repo = Split-Path -Parent $PSScriptRoot
Push-Location $repo
try {
    $python = (Get-Command python -ErrorAction Stop).Source
    $cmd = @((Join-Path $repo 'scripts\RunGeneralT3QuarticConventionC5.py'), '--input', $InputConfig, '--output', $Output)
    if ($ConventionOnly) { $cmd += '--convention-only' }
    & $python @cmd
    $code = $LASTEXITCODE
    if ($code -ne 0) {
        throw ('Stage 15 failed with exit code {0}; inspect {1}\stage15_summary.json' -f $code, $Output)
    }
    Write-Host 'Stage 15 source convention and conditional B/C(-1) estimates completed.'
    Write-Host 'Not a threshold-matched or complete two-loop neutrino-mass calculation.'
}
finally { Pop-Location }
