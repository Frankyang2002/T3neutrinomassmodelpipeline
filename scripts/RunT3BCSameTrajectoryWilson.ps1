param([double]$Rtol = 5e-7)
$ErrorActionPreference = 'Stop'
Set-Location (Split-Path -Parent $PSScriptRoot)
python -m unittest discover -s tests -p 'test_t3_bc_same_trajectory_wilson_v1.py' -v
if ($LASTEXITCODE -ne 0) { throw "Tests failed (exit $LASTEXITCODE)" }
python -u -X utf8 -m scripts.RunT3BCSameTrajectoryWilson --rtol $Rtol
if ($LASTEXITCODE -ne 0) { throw "Same-trajectory Wilson run failed (exit $LASTEXITCODE)" }
