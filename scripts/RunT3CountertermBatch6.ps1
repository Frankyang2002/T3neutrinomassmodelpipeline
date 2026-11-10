param(
    [switch]$Resume,
    [switch]$Overwrite,
    [switch]$Show,
    [switch]$SkipLowEnergyRunning
)
$ErrorActionPreference = 'Stop'
Set-Location (Split-Path -Parent $PSScriptRoot)
$exe = Join-Path (Get-Location).Path '.venv\Scripts\python.exe'
if (-not (Test-Path $exe)) { throw 'Activate/install the project .venv first.' }
$arguments = @('scripts\RunT3CountertermBatch6.py')
if ($Resume) { $arguments += '--resume' }
if ($Overwrite) { $arguments += '--overwrite' }
if ($Show) { $arguments += '--show' }
if ($SkipLowEnergyRunning) { $arguments += '--skip-low-energy-running' }
& $exe @arguments
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
