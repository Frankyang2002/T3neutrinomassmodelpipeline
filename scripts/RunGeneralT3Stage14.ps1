<# Stage 14: conditional UV RGBeta -> extra Matchete-C5 numerical size estimate.
   Only reads prior outputs; does not change production models or RGE code. #>
param(
    [string]$Config = 'docs\general_c5_toy_benchmark.json',
    [string]$Output = 'output\general_t3_stage14'
)
$ErrorActionPreference = 'Stop'
$repo = Split-Path -Parent $PSScriptRoot
Push-Location $repo
try {
    $python = (Get-Command python -ErrorAction Stop).Source
    & $python (Join-Path $repo 'scripts\RunGeneralT3Stage14.py') --input $Config --output $Output
    $code = $LASTEXITCODE
    if ($code -ne 0) { throw ('Stage 14 failed with exit code {0}; inspect {1}\stage14_summary.json' -f $code, $Output) }
    Write-Host 'Stage 14 conditional sensitivity estimates completed; physical coupling mapping remains unverified.'
}
finally { Pop-Location }
