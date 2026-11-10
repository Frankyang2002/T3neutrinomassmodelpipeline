param(
    [switch]$Fresh,
    [switch]$WithTests,
    [switch]$BoundaryOnly,
    [string]$Model = ''
)
$ErrorActionPreference = 'Stop'
Set-Location (Split-Path $PSScriptRoot -Parent)
$py = '.\.venv\Scripts\python.exe'
if (-not (Test-Path $py)) { throw 'Missing project .venv Python' }
if ($WithTests) {
    & $py -m pytest -q tests\test_t3_fermion_hard_triangle_batch1.py
    if ($LASTEXITCODE -ne 0) { throw 'Batch 1 tests failed' }
}
$argsBatch = @('scripts\RunT3FHardTriangleBatch1.py')
if ($Fresh) { $argsBatch += '--overwrite' } else { $argsBatch += '--resume' }
if ($BoundaryOnly) { $argsBatch += '--boundary-only' }
if ($Model) { $argsBatch += @('--model', $Model) }
& $py @argsBatch
if ($LASTEXITCODE -ne 0) { throw "F-hard triangle Batch 1 failed with exit code $LASTEXITCODE" }
& $py scripts\RunT3FHardTriangleBatch1.py --show
if ($LASTEXITCODE -ne 0) { throw 'Batch 1 summary read failed' }
