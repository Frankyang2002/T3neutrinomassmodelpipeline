<# Stage 13: complete UV RGBeta numerical running for the already-defined B/C(-1) sector.
   No production files are changed and no Matchete couplings are identified with RGBeta couplings. #>
param(
    [string]$Config = 'docs\general_c5_toy_benchmark.json',
    [string]$Output = 'output\general_t3_stage13',
    [switch]$TensorOnly
)
$ErrorActionPreference = 'Stop'
$repo = Split-Path -Parent $PSScriptRoot
Push-Location $repo
try {
    $python = (Get-Command python -ErrorAction Stop).Source
    $argsList = @((Join-Path $repo 'scripts\RunGeneralT3CoupledUVRGE.py'), '--config', $Config, '--output', $Output)
    if ($TensorOnly) { $argsList += '--tensor-only' }
    & $python @argsList
    $code = $LASTEXITCODE
    if ($code -ne 0) { throw ('Stage 13 failed with exit code {0}; see {1}\stage13_summary.json' -f $code, $Output) }
    Write-Host 'Stage 13 coupled RGBeta UV calculation completed; production benchmark unchanged.'
}
finally { Pop-Location }
