<#
Preview or remove ONLY known superseded generated report displays.
Project path: scripts/RetireOldStudyDisplays.ps1

From project root:
  .\scripts\RetireOldStudyDisplays.ps1         # Preview, never deletes
  .\scripts\RetireOldStudyDisplays.ps1 -Apply  # Delete after both new figure sets exist
  .\scripts\RetireOldStudyDisplays.ps1 -Apply -WhatIf

Saved scan JSON/CSV, numerical benchmarks, and symbolic artifacts are protected.
#>
[CmdletBinding(SupportsShouldProcess = $true, ConfirmImpact = 'Medium')]
param(
    [switch]$Apply,
    [ValidatePattern('^[A-Za-z0-9_-]+$')]
    [string]$Scenario = 'smallY_smallL'
)

$ErrorActionPreference = 'Stop'
$ProjectRoot = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
if (-not (Test-Path -LiteralPath (Join-Path $ProjectRoot 'pipeline.py') -PathType Leaf)) {
    throw 'Refusing cleanup: missing pipeline.py (invalid T3 project root).'
}
$ReportRoot = Join-Path $ProjectRoot 'Reports\output'
$NewStage = Join-Path $ReportRoot ("display\fermion_threshold_variation\{0}" -f $Scenario)
$NewWeinberg = Join-Path $ReportRoot ("display\weinberg_running\{0}" -f $Scenario)
$Observables = @('delta_m21_sq_ev2', 'delta_m3l_sq_ev2',
                 'sin2_theta12', 'sin2_theta13', 'sin2_theta23')

function Assert-NewFigures([string]$Folder, [string]$Prefix) {
    foreach ($Observable in $Observables) {
        $Image = Join-Path $Folder ("{0}{1}.png" -f $Prefix, $Observable)
        if (-not (Test-Path -LiteralPath $Image -PathType Leaf)) {
            throw "Refusing cleanup: regenerated figure is missing: $Image"
        }
        if ((Get-Item -LiteralPath $Image).Length -lt 1024) {
            throw "Refusing cleanup: regenerated figure is too small: $Image"
        }
    }
}

Write-Host 'Retired display cleanup - preview by default.'
if ($Apply) {
    Assert-NewFigures $NewStage 'fermion_threshold_running_effect_'
    Assert-NewFigures $NewWeinberg 'weinberg_running_'
    Write-Host 'Verified 5 new fermion-threshold and 5 new Weinberg PNGs.'
} else {
    Write-Host 'Use -Apply only after verifying the newly rendered plots.'
}

# Delete ONLY these generated directories. Never traverse numerical data roots.
$OldDirectories = @(
    (Join-Path $ReportRoot 'display\study_comparisons'),
    (Join-Path $ReportRoot 'display\rebuilt_studies'),
    (Join-Path $ReportRoot 'display_stage5\rebuilt_studies')
)
$OldPlots = @(
    @{ Folder = (Join-Path $ReportRoot ("full\staged_rg_comparison\{0}" -f $Scenario)); Pattern = 'staged_rg_*.png' },
    @{ Folder = (Join-Path $ReportRoot ("full\weinberg_running_comparison\{0}" -f $Scenario)); Pattern = 'weinberg_running_*.png' }
)
$OldHtml = @(
    (Join-Path $ReportRoot 'display\interactive_correlations.html'),
    (Join-Path $ReportRoot 'display_stage5\interactive_correlations.html')
)

$Candidates = New-Object System.Collections.Generic.List[System.IO.FileSystemInfo]
foreach ($Dir in $OldDirectories) {
    if (-not (Test-Path -LiteralPath $Dir -PathType Container)) { continue }
    $Protected = @(Get-ChildItem -LiteralPath $Dir -Recurse -File | Where-Object {
        $_.Name -match '^(scan|comparison_summary)\.(json|csv)$'
    })
    if ($Protected.Count -ne 0) {
        throw "Refusing to retire $Dir because saved numerical scan/summary data was found inside it."
    }
    $Candidates.Add((Get-Item -LiteralPath $Dir))
}
foreach ($Spec in $OldPlots) {
    if (-not (Test-Path -LiteralPath $Spec.Folder -PathType Container)) { continue }
    foreach ($File in @(Get-ChildItem -LiteralPath $Spec.Folder -File -Filter $Spec.Pattern)) {
        $Candidates.Add($File)
    }
}
foreach ($File in $OldHtml) {
    if (Test-Path -LiteralPath $File -PathType Leaf) {
        $Candidates.Add((Get-Item -LiteralPath $File))
    }
}

if ($Candidates.Count -eq 0) {
    Write-Host 'Nothing retired remains. No files changed.'
    return
}
foreach ($Item in $Candidates) {
    Write-Host ("[RETIRE {0}] {1}" -f $(if ($Item.PSIsContainer) {'DIRECTORY'} else {'FILE'}), $Item.FullName)
    if ($Apply -and $PSCmdlet.ShouldProcess($Item.FullName, 'Remove known obsolete generated display')) {
        Remove-Item -LiteralPath $Item.FullName -Recurse -Force
    }
}
if (-not $Apply) {
    Write-Host ("Preview: {0} candidates. Nothing was deleted. Re-run with -Apply." -f $Candidates.Count)
} else {
    Write-Host 'Cleanup complete (subject to -WhatIf/-Confirm). Source scan JSON/CSV was not touched.'
}
