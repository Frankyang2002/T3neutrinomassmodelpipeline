# Run the entire repository test suite with deterministic Unicode text I/O on Windows.
# Does not alter physical model inputs, source files, or saved numerical outputs.
param(
    [switch]$VerboseTrace
)
$ErrorActionPreference = 'Stop'
$root = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
Push-Location $root
try {
    $env:PYTHONUTF8 = '1'
    python -c "import sys; assert sys.flags.utf8_mode == 1, 'Python UTF-8 mode is not enabled'; print('Python UTF-8 mode: enabled')"
    if ($LASTEXITCODE -ne 0) { throw 'Python UTF-8 preflight failed' }
    if ($VerboseTrace) {
        python -m pytest -q --tb=short
    } else {
        python -m pytest -q --tb=line
    }
    if ($LASTEXITCODE -ne 0) {
        throw "Pytest failed with exit code $LASTEXITCODE. Preserve the full output for triage."
    }
}
finally {
    Pop-Location
}
