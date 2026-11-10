param(
  [switch]$Resume,
  [switch]$Overwrite,
  [switch]$SkipLowFeedback,
  [double]$LambdaHTransferFactor = 1.0
)
$ErrorActionPreference = 'Stop'
$project = Split-Path -Parent $PSScriptRoot
Set-Location $project
$py = Join-Path $project '.venv\Scripts\python.exe'
if (!(Test-Path $py)) { throw 'Missing .venv Python executable' }
$argv = @('scripts\RunT3ScalarCWBatch3.py', '--lambdaH-transfer-factor', "$LambdaHTransferFactor")
if ($Resume) { $argv += '--resume' }
if ($Overwrite) { $argv += '--overwrite' }
if ($SkipLowFeedback) { $argv += '--skip-low-feedback' }
& $py @argv
if ($LASTEXITCODE -ne 0) { throw "Batch3 failed with exit code $LASTEXITCODE" }
& $py 'scripts\RunT3ScalarCWBatch3.py' '--show'
if ($LASTEXITCODE -ne 0) { throw "Batch3 show failed with exit code $LASTEXITCODE" }
