$ErrorActionPreference = 'Stop'
Set-Location (Split-Path -Parent $PSScriptRoot)
python -m unittest discover -s tests -p 'test_t3_bc_convention_matching_batch_v1.py'
if ($LASTEXITCODE -ne 0) { throw "Unit tests failed (exit $LASTEXITCODE)" }
python -u -X utf8 -m scripts.RunT3BCConventionMatchingBatch
if ($LASTEXITCODE -ne 0) { throw "Convention matching batch failed (exit $LASTEXITCODE)" }
