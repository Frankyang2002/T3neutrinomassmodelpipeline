# One-command T3 benchmark reproducibility wrapper. No GitHub operations.
[CmdletBinding()]
param(
  [switch]$CollectOnly,
  [switch]$VerifyOnly,
  [switch]$SkipTests
)
$ErrorActionPreference = 'Stop'
$root = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$python = Join-Path $root '.venv\Scripts\python.exe'
if (-not (Test-Path -LiteralPath $python)) { throw "Missing Python venv: $python" }
$argsList = @((Join-Path $root 'scripts\RunT3FinalRepro.py'), '--root', $root)
if ($CollectOnly) { $argsList += '--collect-only' }
if ($VerifyOnly) { $argsList += '--verify-only' }
if ($SkipTests) { $argsList += '--skip-tests' }
& $python @argsList
if ($LASTEXITCODE -ne 0) { throw "T3 final reproducibility run failed (exit $LASTEXITCODE)" }
