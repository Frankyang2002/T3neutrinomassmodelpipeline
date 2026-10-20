<#
.SYNOPSIS
  Audit Reports/output/full and optionally retire verified old studies, zero-byte
  files, and empty directories. Local only; preview is the default.
.DESCRIPTION
  Retired numerical directories are moved OUT of Reports/output/full into a
  dated .cleanup_quarantine/reports_full folder, not destroyed. This permits
  recovery of historical saved scans. Zero-byte files and superseded rendered
  PNGs are deleted only with -Apply. Current numerical studies are untouched.
  Scan JSON/CSV in retained studies are never targeted by content or extension.
  If Windows locks an otherwise empty directory, warn and continue; do not
  force removal or interrupt cleanup of other eligible empty directories.

  Existing scripts/CleanReportOutput.ps1 is substantially broader: do not use
  it for selective Reports/output/full cleanup.
.EXAMPLE
  .\scripts\AuditAndPruneReportsFull.ps1
.EXAMPLE
  .\scripts\AuditAndPruneReportsFull.ps1 -Apply
.EXAMPLE
  .\scripts\AuditAndPruneReportsFull.ps1 -Apply -WhatIf
#>
[CmdletBinding(SupportsShouldProcess = $true, ConfirmImpact = 'Medium')]
param([switch] $Apply)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'
$root = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
$full = Join-Path $root 'Reports\output\full'
$display = Join-Path $root 'Reports\output\display'
if (-not (Test-Path -LiteralPath (Join-Path $root 'pipeline.py') -PathType Leaf)) {
    throw "Not a T3 repository: missing $root\pipeline.py"
}
if (-not (Test-Path -LiteralPath $full -PathType Container)) {
    throw "Missing Reports/output/full at $full"
}
if (((Get-Item -LiteralPath $full).Attributes -band [IO.FileAttributes]::ReparsePoint) -ne 0) {
    throw 'Refusing to process a symbolic-link Reports/output/full.'
}

# ONLY the old four-benchmark, obsolete E-p0 shared-parameter work,
# and superseded four-scenario comparison report trees.
# No active study read by tools/RebuildCurrentReports.py is listed here.
$retiredNames = @(
    'four_benchmark_correlations',
    'common_parameter_comparison',
    'common_parameter_comparison_v2',
    'common_yukawa_optimization',
    'common_yukawa_least_squares',
    'common_yukawa_final_summary',
    'common_yukawa_correlation_scans',
    'common_yukawa_three_panel_scans',
    'scalar_coupling_correlations',
    'benchmark_summary',
    'yukawa_grouping',
    'comparison',
    'within_model_comparison'
)

$retired = @($retiredNames | ForEach-Object {
    $p = Join-Path $full $_
    if (Test-Path -LiteralPath $p -PathType Container) {
        $item = Get-Item -LiteralPath $p
        if (($item.Attributes -band [IO.FileAttributes]::ReparsePoint) -ne 0) {
            throw "Refusing symbolic-link legacy directory: $p"
        }
        $item
    }
})

# Retained numerical studies may include older, misleading PNG renderings.
# Only retire these after validating all 5 CURRENT display PNGs for that study.
$observables = @('delta_m21_sq_ev2', 'delta_m3l_sq_ev2',
                 'sin2_theta12', 'sin2_theta13', 'sin2_theta23')
$studies = @(
    @{ Source='staged_rg_comparison'; Destination='fermion_threshold_variation'; Old='staged_rg_*.png'; New='fermion_threshold_running_effect_' },
    @{ Source='weinberg_running_comparison'; Destination='weinberg_running'; Old='weinberg_running_*.png'; New='weinberg_running_' }
)
$oldPngs = New-Object System.Collections.Generic.List[System.IO.FileInfo]
$deferred = New-Object System.Collections.Generic.List[string]
foreach ($study in $studies) {
    $source = Join-Path $full $study.Source
    if (-not (Test-Path -LiteralPath $source -PathType Container)) { continue }
    foreach ($scenarioDir in @(Get-ChildItem -LiteralPath $source -Directory -Force)) {
        if (($scenarioDir.Attributes -band [IO.FileAttributes]::ReparsePoint) -ne 0) { continue }
        $matches = @(Get-ChildItem -LiteralPath $scenarioDir.FullName -File -Filter $study.Old)
        if ($matches.Count -eq 0) { continue }
        $newDir = Join-Path (Join-Path $display $study.Destination) $scenarioDir.Name
        $ready = $true
        foreach ($observable in $observables) {
            $png = Join-Path $newDir ($study.New + $observable + '.png')
            if (-not (Test-Path -LiteralPath $png -PathType Leaf)) { $ready = $false; break }
            if ((Get-Item -LiteralPath $png).Length -lt 1024) { $ready = $false; break }
        }
        if ($ready) {
            foreach ($png in $matches) { $oldPngs.Add($png) }
        } else {
            $deferred.Add("$($study.Source)/$($scenarioDir.Name): new five-figure set missing; old PNGs kept")
        }
    }
}

# Zero-byte files are never usable saved numerical data. Do not follow links.
$emptyFiles = @(Get-ChildItem -LiteralPath $full -Recurse -File -Force |
    Where-Object { $_.Length -eq 0 -and
        (($_.Attributes -band [IO.FileAttributes]::ReparsePoint) -eq 0) })
# Pre-existing empty directories; more may become empty after -Apply.
$emptyDirs = @(Get-ChildItem -LiteralPath $full -Directory -Recurse -Force |
    Where-Object {
        (($_.Attributes -band [IO.FileAttributes]::ReparsePoint) -eq 0) -and
        ($null -eq (Get-ChildItem -LiteralPath $_.FullName -Force | Select-Object -First 1))
    })

Write-Host "Reports/full audit at: $full"
Write-Host $(if ($Apply) { 'Mode: APPLY (explicit -Apply requested)' } else { 'Mode: PREVIEW ONLY' })
Write-Host "Retired folder candidates: $($retired.Count)"
foreach ($dir in $retired) { Write-Host "  [QUARANTINE] $($dir.FullName)" }
Write-Host "Superseded PNG candidates: $($oldPngs.Count)"
foreach ($png in $oldPngs) { Write-Host "  [DELETE OLD PNG] $($png.FullName)" }
Write-Host "Zero-byte file candidates: $($emptyFiles.Count)"
foreach ($file in $emptyFiles) { Write-Host "  [DELETE EMPTY FILE] $($file.FullName)" }
Write-Host "Already-empty folder candidates: $($emptyDirs.Count)"
foreach ($dir in $emptyDirs) { Write-Host "  [PRUNE EMPTY DIR] $($dir.FullName)" }
foreach ($warning in $deferred) { Write-Warning "DEFERRED $warning" }
if (-not $Apply) {
    Write-Host 'No changes made. Use -Apply to quarantine retired studies and prune empties.'
    return
}
if ($WhatIfPreference) {
    Write-Host 'WhatIf: no changes or quarantine created.'
    return
}

$archiveRoot = Join-Path $root '.cleanup_quarantine\reports_full'
$stamp = (Get-Date).ToString('yyyyMMdd_HHmmss')
$archive = Join-Path $archiveRoot $stamp
$moved = 0
foreach ($dir in $retired) {
    if (-not (Test-Path -LiteralPath $dir.FullName -PathType Container)) { continue }
    if ($PSCmdlet.ShouldProcess($dir.FullName, 'Move retired historical report folder into local quarantine')) {
        if (-not (Test-Path -LiteralPath $archive -PathType Container)) {
            New-Item -ItemType Directory -Path $archive -Force | Out-Null
        }
        $destination = Join-Path $archive $dir.Name
        if (Test-Path -LiteralPath $destination) { throw "Quarantine collision: $destination" }
        Move-Item -LiteralPath $dir.FullName -Destination $destination -ErrorAction Stop
        $moved++
    }
}
$deleted = 0
foreach ($f in @($oldPngs) + @($emptyFiles)) {
    # Files under retired directories have already moved; skip them.
    if (-not (Test-Path -LiteralPath $f.FullName -PathType Leaf)) { continue }
    if ($PSCmdlet.ShouldProcess($f.FullName, 'Remove superseded plot or zero-byte file')) {
        Remove-Item -LiteralPath $f.FullName -Force -ErrorAction Stop
        $deleted++
    }
}
$pruned = 0
$skippedEmptyDirs = New-Object System.Collections.Generic.List[string]
# Delete ONLY directories which are genuinely empty AFTER removing files.
# Deepest first. Never use -Recurse on Remove-Item for this step.
$allDirs = @(Get-ChildItem -LiteralPath $full -Directory -Recurse -Force |
    Sort-Object { $_.FullName.Length } -Descending)
foreach ($dir in $allDirs) {
    $path = $dir.FullName
    if (-not (Test-Path -LiteralPath $path -PathType Container)) { continue }
    try {
        # Refresh metadata: another process may have changed the directory.
        $current = Get-Item -LiteralPath $path -Force -ErrorAction Stop
        if (($current.Attributes -band [IO.FileAttributes]::ReparsePoint) -ne 0) { continue }
        $first = Get-ChildItem -LiteralPath $path -Force -ErrorAction Stop |
            Select-Object -First 1
    }
    catch {
        $skippedEmptyDirs.Add($path)
        Write-Warning ("Could not inspect directory; kept unchanged: {0} ({1})" -f
            $path, $_.Exception.Message)
        continue
    }
    if ($null -ne $first) { continue }
    if (-not $PSCmdlet.ShouldProcess($path, 'Remove now-empty folder')) { continue }

    # Windows may temporarily lock an empty directory (Explorer, VS Code,
    # an open terminal, antivirus, or another process). A short retry handles
    # transient locks. A persistent lock is NOT overridden or treated as a
    # fatal error: report it for a later run instead.
    for ($attempt = 1; $attempt -le 2; $attempt++) {
        try {
            Remove-Item -LiteralPath $path -Force -ErrorAction Stop
            $pruned++
            break
        }
        catch {
            if ($attempt -lt 2) {
                Start-Sleep -Milliseconds 200
            }
            else {
                $skippedEmptyDirs.Add($path)
                Write-Warning ("Empty directory remains locked/inaccessible; skipped: {0} ({1})" -f
                    $path, $_.Exception.Message)
            }
        }
    }
}
Write-Host "Moved $moved retired folders to quarantine: $archive"
Write-Host "Deleted $deleted obsolete PNGs/zero-byte files; pruned $pruned empty folders."
if ($skippedEmptyDirs.Count -gt 0) {
    Write-Warning (("{0} directory/directories could not be removed or inspected. " +
        "Close programs using those paths and rerun -Apply later.") -f $skippedEmptyDirs.Count)
    foreach ($path in $skippedEmptyDirs) { Write-Host "  [LEFT UNCHANGED] $path" }
} else {
    Write-Host 'All accessible empty folders were pruned.'
}
Write-Host 'Retained numerical JSON/CSV and fitted configuration folders were not touched.'
Write-Host 'No GitHub commands were used.'
