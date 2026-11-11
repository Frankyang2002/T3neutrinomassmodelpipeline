$ErrorActionPreference = 'Stop'
Set-Location (Split-Path -Parent $PSScriptRoot)
python -X utf8 -m scripts.TestT3MatcheteCartesianResponses
if ($LASTEXITCODE -ne 0) { throw 'Cartesian exact-algebra unit tests failed' }
$source = 'output/general_t3_stage22/matchete_tensors.json'
if (-not (Test-Path $source)) {
    Write-Host 'Native Matchete evidence missing. Export it with:'
    Write-Host 'wolframscript -file scripts/ExportGeneralT3ScalarTensorMatchete.wl output/general_t3_stage22/matchete_tensors.json'
    Write-Host 'Synthetic algebra tests passed; full source export skipped.'
    exit 0
}
python -X utf8 -m scripts.BuildT3MatcheteCartesianResponses
if ($LASTEXITCODE -ne 0) { throw 'Native Matchete candidate conversion failed' }
