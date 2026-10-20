<#
Retire obsolete/generated display folders after the new figures have been built.
Project path: scripts/RetireOldStudyDisplays.ps1

PowerShell from the project root:
    .\scripts\RetireOldStudyDisplays.ps1        # Preview only
    .\scripts\RetireOldStudyDisplays.ps1 -Apply # Actually delete

Does not delete or alter saved numerical scan.json or scan.csv files.
#>
[CmdletBinding(SupportsShouldProcess=$true, ConfirmImpact='Medium')]
param(
    [switch]$Apply,
    [string]$Scenario = 'smallY_smallL'
)

$ErrorActionPreference = 'Stop'
$ProjectRoot = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
if (-not (Test-Path -LiteralPath (Join-Path $ProjectRoot 'pipeline.py') -PathType Leaf)) {
    throw 'Not a T3 repository root: pipeline.py is missing.'
}
if ($Scenario -notmatch '^[A-Za-z0-9_-]+$') {
    throw 'Invalid scenario directory identifier.'
}

$ReportRoot = Join-Path $ProjectRoot 'Reports\output'
$NewStage = Join-Path $ReportRoot ("display\fermion_threshold_variation\{0}" -f $Scenario)
$NewWeinberg = Join-Path $ReportRoot ("display\weinberg_running\{0}" -f $Scenario)

function Get-FigureCount([string]$Path, [string]$Pattern) {
    if (-not (Test-Path -LiteralPath $Path -PathType Container)) { return 0 }
    return @(Get-ChildItem -LiteralPath $Path -Filter $Pattern -File).Count
}
$StageCount = Get-FigureCount $NewStage 'fermion_threshold_running_effect_*.png'
$WeinbergCount = Get-FigureCount $NewWeinberg 'weinberg_running_*.png'
Write-Host ("New figures found: fermion-threshold={0}/5; Weinberg={1}/5" -f $StageCount, $WeinbergCount)
if ($Apply -and ($StageCount -lt 5 -or $WeinbergCount -lt 5)) {
    throw 'Refusing to remove old plots: first regenerate the 5 new figures in each study.'
}

$OldDirectories = @(
    (Join-Path $ReportRoot 'display\study_comparisons'),
    (Join-Path $ReportRoot 'display\rebuilt_studies'),
    (Join-Path $ReportRoot 'display_stage5\rebuilt_studies')
)
$OldPlots = @(
    @{ Folder = (Join-Path $ReportRoot ("full\staged_rg_comparison\{0}" -f $Scenario)); Pattern = 'staged_rg_*.png' },
    @{ Folder = (Join-Path $ReportRoot ("full\weinberg_running_comparison\{0}" -f $Scenario)); Pattern = 'weinberg_running_*.png' }
)
$OldHtml = Join-Path $ReportRoot 'display\interactive_correlations.html'

foreach ($dir in $OldDirectories) {
    if (-not (Test-Path -LiteralPath $dir -PathType Container)) { continue }
    Write-Host "[RETIRE FOLDER] $dir"
    if ($Apply -and $PSCmdlet.ShouldProcess($dir, 'Delete retired generated display directory')) {
        Remove-Item -LiteralPath $dir -Force -Recurse
    }
}
foreach ($entry in $OldPlots) {
    if (-not (Test-Path -LiteralPath $entry.Folder -PathType Container)) { continue }
    foreach ($item in @(Get-ChildItem -LiteralPath $entry.Folder -File -Filter $entry.Pattern)) {
        Write-Host "[RETIRE OLD PLOT] $($item.FullName)"
        if ($Apply -and $PSCmdlet.ShouldProcess($item.FullName, 'Delete old rendered PNG only')) {
            Remove-Item -LiteralPath $item.FullName -Force
        }
    }
}
if (Test-Path -LiteralPath $OldHtml -PathType Leaf) {
    Write-Host "[RETIRE HTML] $OldHtml"
    if ($Apply -and $PSCmdlet.ShouldProcess($OldHtml, 'Delete retired correlation HTML')) {
        Remove-Item -LiteralPath $OldHtml -Force
    }
}
if (-not $Apply) {
    Write-Host 'Preview only. Re-run with -Apply after rebuilding the plots.'
} else {
    Write-Host 'Retired generated displays removed. Original numerical scan JSON/CSV retained.'
}
