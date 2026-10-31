$ErrorActionPreference = 'Stop'
Set-Location (Split-Path -Parent $PSScriptRoot)
python -m unittest discover -s tests -p 'test_t3_bc_fermion_threshold_batch_v1.py'
if ($LASTEXITCODE -ne 0) { throw "F-threshold tests failed ($LASTEXITCODE)" }
python -u -X utf8 -m scripts.RunT3BCFermionThresholdBatch @args
if ($LASTEXITCODE -ne 0) { throw "F-threshold batch failed ($LASTEXITCODE); see native output prerequisites and logs" }
