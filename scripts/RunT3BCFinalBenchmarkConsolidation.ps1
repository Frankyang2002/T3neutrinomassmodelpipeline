$ErrorActionPreference = 'Stop'
Set-Location (Split-Path -Parent $PSScriptRoot)
python -m unittest discover -s tests -p 'test_t3_bc_final_benchmark_consolidation_v1.py'
if ($LASTEXITCODE -ne 0) { throw "Unit tests failed ($LASTEXITCODE)" }
python -u -X utf8 -m scripts.RunT3BCFinalBenchmarkConsolidation @args
if ($LASTEXITCODE -ne 0) { throw "B/C consolidation failed ($LASTEXITCODE)" }
