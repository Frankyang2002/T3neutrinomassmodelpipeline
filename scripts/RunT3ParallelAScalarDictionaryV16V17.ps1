# v16/v17 all-16 native Matchete nonportal export and complete conditional real census.
# No Matchete Match[], no RGBeta beta evolution, no EFT/Wilson running.
param([switch]$VerifyOnly)
$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $PSScriptRoot
Set-Location -LiteralPath $root
$python = Join-Path $root '.venv\Scripts\python.exe'
if (!(Test-Path -LiteralPath $python)) { throw 'Missing project .venv Python: ' + $python }
function RequireExit([string]$label) {
    if ($LASTEXITCODE -ne 0) { throw "$label blocked with exit code $LASTEXITCODE" }
}
$base = 'output\parallel_t3\wilson'
$native = Join-Path $base 'matchete_nonportal_native_v16.json'
$mapped = Join-Path $base 'nonportal_dictionary_conditional_v16.json'
$census = Join-Path $base 'full_scalar_census_conditional_v17.json'
Write-Host 'v16: standalone + integration tests'
& $python -m unittest discover -s tests -p 'test_t3_parallel_a_nonportal_dictionary_v16.py' -v
RequireExit 'v16 tests'
& $python -m unittest discover -s tests -p 'test_t3_parallel_a_full_scalar_census_v17.py' -v
RequireExit 'v17 tests'
if ($VerifyOnly -and (!(Test-Path -LiteralPath $mapped) -or !(Test-Path -LiteralPath $census))) {
    throw 'VerifyOnly requires existing v16 dictionary and v17 census JSON files'
}
if (!(Test-Path -LiteralPath $native)) {
    if ($VerifyOnly) { throw 'v16 native Matchete source does not yet exist; cannot verify' }
    if (Test-Path -LiteralPath ($native + '.partial')) { throw 'v16 native partial exists; investigate before rerunning' }
    Write-Host 'v16 native Matchete: exporting 76 registered self/cross CG tensors across 16 models'
    & wolframscript -file scripts\RunT3ParallelAExportNonportalMatcheteV16.wl $native
    RequireExit 'v16 native CG export'
} else {
    Write-Host 'v16 native JSON already exists; revalidate its source and tensor contents'
}
& $python -m scripts.RunT3ParallelANonportalDictionaryV16 --project-root . --resume
RequireExit 'v16 Matchete-to-RGBeta conditional mapping'
& $python -m scripts.RunT3ParallelAFullScalarCensusV17 --project-root . --resume
RequireExit 'v17 full real-scalar census'
Write-Host 'v16/v17 source-conditioned 150-real-slot census complete; NO physical Wilson evolution.'
