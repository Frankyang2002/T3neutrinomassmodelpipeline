$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $PSScriptRoot
Push-Location $root
try {
    python -m unittest discover -s tests -p 'test_t3_bc_majorana_action_batch_v1.py'
    if ($LASTEXITCODE -ne 0) { throw "Unit tests failed: $LASTEXITCODE" }
    python -m scripts.RunT3BCMajoranaActionBatch @args
    if ($LASTEXITCODE -ne 0) { throw "Action batch failed: $LASTEXITCODE" }
} finally { Pop-Location }
