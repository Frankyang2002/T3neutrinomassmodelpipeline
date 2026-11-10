# Install this complete file as scripts\RunT3NativeFiniteBatch4.ps1
$ErrorActionPreference = 'Stop'
Set-Location (Split-Path -Parent $PSScriptRoot)
$py = '.\.venv\Scripts\python.exe'
& $py -m pytest -q tests\test_t3_native_finite_batch4.py
if ($LASTEXITCODE -ne 0) { throw 'Batch 4 unit tests failed.' }
& $py scripts\RunT3NativeFiniteBatch4.py --resume
if ($LASTEXITCODE -ne 0) { throw 'Batch 4 all-16 run has blocked models.' }
& $py scripts\RunT3NativeFiniteBatch4.py --show
if ($LASTEXITCODE -ne 0) { throw 'Batch 4 summary unavailable.' }
