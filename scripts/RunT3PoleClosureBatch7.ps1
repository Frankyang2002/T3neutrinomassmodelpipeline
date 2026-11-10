param(
    [ValidateSet('All','B_-1','A_-4','A_-2','A_+0','B_-3','B_+1','C_-3','C_-1','C_+1','D_-2','D_+0','D_+2','E_-4','E_-2','E_+0','E_+2')]
    [string]$Model = 'All',
    [switch]$Resume,
    [switch]$Overwrite,
    [switch]$SkipLowEnergyRunning,
    [switch]$Show
)
$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
Set-Location $project
$python = Join-Path $project '.venv\Scripts\python.exe'
if (-not (Test-Path -LiteralPath $python)) { throw "Missing project venv python: $python" }
$runner = Join-Path $project 'scripts\RunT3PoleClosureBatch7.py'
$args7 = @($runner)
if ($Model -ne 'All') { $args7 += @('--model', $Model) }
if ($Resume) { $args7 += '--resume' }
if ($Overwrite) { $args7 += '--overwrite' }
if ($SkipLowEnergyRunning) { $args7 += '--skip-low-energy-running' }
if ($Show) { $args7 += '--show' }
& $python @args7
if ($LASTEXITCODE -ne 0) { throw "Batch7 returned exit code $LASTEXITCODE" }
