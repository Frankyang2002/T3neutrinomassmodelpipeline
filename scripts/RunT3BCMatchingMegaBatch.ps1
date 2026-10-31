param([switch]$SkipConvergence)
$ErrorActionPreference = 'Stop'
Set-Location (Split-Path -Parent $PSScriptRoot)
& python -m unittest -v tests.test_t3_bc_matching_megabatch_v1
if ($LASTEXITCODE -ne 0) { throw "Matching mega-batch tests failed (exit $LASTEXITCODE)" }
$opts=@()
if ($SkipConvergence) { $opts += '--skip-convergence' }
& python -u -X utf8 -m scripts.RunT3BCMatchingMegaBatch @opts
if ($LASTEXITCODE -ne 0) { throw "Matching mega-batch numerical run failed (exit $LASTEXITCODE)" }
