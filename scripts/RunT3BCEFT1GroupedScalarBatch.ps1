$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $PSScriptRoot
Push-Location $root
try {
    python -m unittest discover -s tests -p 'test_t3_bc_eft1_grouped_scalar_batch_v1.py'
    if ($LASTEXITCODE -ne 0) { throw "Tests failed (exit $LASTEXITCODE)" }
    python -u -X utf8 -m scripts.RunT3BCEFT1GroupedScalarBatch @args
    if ($LASTEXITCODE -ne 0) { throw "EFT1 grouped scalar batch failed (exit $LASTEXITCODE)" }
} finally { Pop-Location }
