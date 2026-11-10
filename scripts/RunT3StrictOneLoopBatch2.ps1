# Run from the local repository root. The script never modifies GitHub.
$ErrorActionPreference = 'Stop'
if (-not (Test-Path '.\.venv\Scripts\python.exe')) {
    throw 'Activate/install the project .venv and invoke from the repository root.'
}
& .\.venv\Scripts\python.exe -m pytest -q tests\test_t3_strict_one_loop_batch2.py
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
& .\.venv\Scripts\python.exe scripts\RunT3StrictOneLoopBatch2.py --resume
exit $LASTEXITCODE
