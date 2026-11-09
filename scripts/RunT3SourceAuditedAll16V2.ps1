# Run from project root; keep the 10 successful v1 results by default.
param(
    [switch]$Fresh,
    [switch]$BoundaryOnly,
    [switch]$WithTests
)
$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
Set-Location $project
$python = if (Test-Path '.\.venv\Scripts\python.exe') { '.\.venv\Scripts\python.exe' } else { 'python' }
if ($WithTests) {
    & $python -m pytest -q tests\test_t3_source_audited_all16.py tests\test_t3_source_audited_all16_v2.py
    if ($LASTEXITCODE -ne 0) { throw 'Source-map tests failed. No numerical batch executed.' }
}
$summary = Join-Path $project 'output\parallel_t3\source_audited_all16_v1\summary.json'
$arguments = @('scripts\RunT3SourceAuditedAll16.py')
if ($Fresh) { $arguments += '--overwrite' }
elseif (Test-Path $summary) { $arguments += '--resume' }
if ($BoundaryOnly) { $arguments += '--boundary-only' }
& $python @arguments
$exitCode=$LASTEXITCODE
& $python scripts\RunT3SourceAuditedAll16.py --show
if ($exitCode -ne 0) { throw "All-16 numerical batch reports blocked models; exit=$exitCode. See summary.json." }
