# Save under scripts/RunT3HardSoftSubtractionBatch5.ps1, execute at project root.
param(
  [switch] $Resume,
  [switch] $Overwrite,
  [switch] $AuditOnly,
  [switch] $Show,
  [string] $Model
)
$ErrorActionPreference = 'Stop'
$project = (Resolve-Path '.').Path
$python = Join-Path $project '.venv\Scripts\python.exe'
if (-not (Test-Path -LiteralPath $python)) { throw 'Missing .venv Python' }
$script = Join-Path $project 'scripts\RunT3HardSoftSubtractionBatch5.py'
$argv = @($script)
if ($Resume) { $argv += '--resume' }
if ($Overwrite) { $argv += '--overwrite' }
if ($AuditOnly) { $argv += '--audit-only' }
if ($Show) { $argv += '--show' }
if ($Model) { $argv += @('--model', $Model) }
& $python @argv
exit $LASTEXITCODE
