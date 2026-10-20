<#
.SYNOPSIS
    Preview or archive superseded per-model fitted candidates without changing
    which candidate best_candidate() selects.

.DESCRIPTION
    Reads the saved feasibility_summary.json best.chi2 and requires the matching
    fitted_candidate_config.json. Uses the exact stage priority implemented in
    Numerical/orchestration/studies/fitted_rge/RunFittedSensitivityComparison.py:
        fit_feasibility -> fit_feasibility_retry -> fit_feasibility_focused
    The lowest finite chi2 wins; ties favour the earlier stage.

    Without -Apply, does not write anything. With -Apply, moves ONLY superseded
    model directories (two known JSON files) into a dated archive below the
    project root's .cleanup_quarantine directory. All winning candidate files
    remain at their original paths, and saved batch/retry summaries are kept.

    No fits, RGEs, or symbolic work are run. No GitHub commands are run.

.EXAMPLE
    .\scripts\ArchiveSupersededFits.ps1
.EXAMPLE
    .\scripts\ArchiveSupersededFits.ps1 -Apply -WhatIf
.EXAMPLE
    .\scripts\ArchiveSupersededFits.ps1 -Apply
#>
[CmdletBinding(SupportsShouldProcess = $true, ConfirmImpact = 'Medium')]
param(
    [switch] $Apply,
    [ValidatePattern('^[A-Za-z0-9_-]+$')]
    [string] $Scenario = 'smallY_smallL'
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'
$ProjectRoot = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
$ReportFull = Join-Path $ProjectRoot 'Reports\output\full'
if (-not (Test-Path -LiteralPath (Join-Path $ProjectRoot 'pipeline.py') -PathType Leaf)) {
    throw 'Not a valid T3 project root: pipeline.py is missing.'
}
if (-not (Test-Path -LiteralPath $ReportFull -PathType Container)) {
    throw "Saved numerical results missing: $ReportFull"
}

$Stages = @('fit_feasibility', 'fit_feasibility_retry', 'fit_feasibility_focused')
$Candidates = [System.Collections.Generic.List[object]]::new()

function Read-FitChi2([string] $SummaryPath) {
    $payload = Get-Content -LiteralPath $SummaryPath -Raw -Encoding UTF8 | ConvertFrom-Json
    if ($null -eq $payload -or $null -eq $payload.best -or $null -eq $payload.best.chi2) {
        throw "Missing best.chi2: $SummaryPath"
    }
    $chi2 = 0.0
    $ok = [double]::TryParse(
        [string] $payload.best.chi2,
        [System.Globalization.NumberStyles]::Float,
        [System.Globalization.CultureInfo]::InvariantCulture,
        [ref] $chi2
    )
    if (-not $ok -or [double]::IsNaN($chi2) -or [double]::IsInfinity($chi2)) {
        throw "Nonfinite or invalid best.chi2: $SummaryPath"
    }
    return $chi2
}

foreach ($stageIndex in 0..($Stages.Count - 1)) {
    $stage = $Stages[$stageIndex]
    $scenarioDir = Join-Path (Join-Path $ReportFull $stage) $Scenario
    if (-not (Test-Path -LiteralPath $scenarioDir -PathType Container)) { continue }
    foreach ($directory in @(Get-ChildItem -LiteralPath $scenarioDir -Directory -Force)) {
        if ($directory.Name -notmatch '^T3_dS1_\d+_dS2_\d+_dF_\d+_alpha_[mp]\d+$') { continue }
        if (($directory.Attributes -band [System.IO.FileAttributes]::ReparsePoint) -ne 0) {
            throw "Refusing to traverse symlink/junction: $($directory.FullName)"
        }
        $summary = Join-Path $directory.FullName 'feasibility_summary.json'
        $config = Join-Path $directory.FullName 'fitted_candidate_config.json'
        if (-not (Test-Path -LiteralPath $summary -PathType Leaf) -or
            -not (Test-Path -LiteralPath $config -PathType Leaf)) {
            Write-Warning "Incomplete fit candidate (left unchanged): $($directory.FullName)"
            continue
        }
        try {
            $chi2 = Read-FitChi2 $summary
        } catch {
            Write-Warning "Unusable fit candidate (left unchanged): $($directory.FullName): $_"
            continue
        }
        $Candidates.Add([pscustomobject]@{
            Model = $directory.Name
            Stage = $stage
            StageIndex = $stageIndex
            Chi2 = [double] $chi2
            Path = $directory.FullName
        })
    }
}

if ($Candidates.Count -eq 0) {
    throw 'No valid fitted candidates found; nothing was archived.'
}

$Groups = @($Candidates | Group-Object Model)
if ($Groups.Count -ne 16) {
    throw "Safety check: expected 16 models with at least one valid fit, found $($Groups.Count). No changes made."
}

$Winners = [System.Collections.Generic.List[object]]::new()
$Losers = [System.Collections.Generic.List[object]]::new()
foreach ($group in $Groups) {
    $ordered = @($group.Group | Sort-Object Chi2, StageIndex)
    $winner = $ordered[0]
    $Winners.Add($winner)
    if ($ordered.Count -lt 2) { continue }

    # If an existing final running diagnostic explicitly cites a different
    # fitting stage, do not silently archive its provenance source.
    $runningDiagnostic = Join-Path (Join-Path (Join-Path $ReportFull 'fitted_rge') $Scenario) (Join-Path $winner.Model 'running_diagnostics.json')
    if (Test-Path -LiteralPath $runningDiagnostic -PathType Leaf) {
        try {
            $diagnostic = Get-Content -LiteralPath $runningDiagnostic -Raw -Encoding UTF8 | ConvertFrom-Json
            $recordedStage = $null
            if ($null -ne $diagnostic.fitted_rge) {
                $recordedStage = $diagnostic.fitted_rge.fit_stage
            }
            if ($null -ne $recordedStage -and [string] $recordedStage -ne $winner.Stage) {
                Write-Warning "Skipping superseded candidates for $($winner.Model): saved fitted_rge diagnostics record stage '$recordedStage' instead of current best '$($winner.Stage)'."
                continue
            }
        } catch {
            Write-Warning "Cannot verify saved RGE provenance for $($winner.Model); leaving its superseded candidates unchanged: $_"
            continue
        }
    }

    foreach ($candidate in $ordered[1..($ordered.Count - 1)]) {
        # Move only well-known two-file fitted-candidate directories.
        $contents = @(Get-ChildItem -LiteralPath $candidate.Path -Force)
        $expectedNames = @('feasibility_summary.json', 'fitted_candidate_config.json')
        if ($contents.Count -ne 2 -or @($contents | Where-Object {
            $_.PSIsContainer -or $_.Name -notin $expectedNames -or
            (($_.Attributes -band [System.IO.FileAttributes]::ReparsePoint) -ne 0)
        }).Count -ne 0) {
            Write-Warning "Extra or unexpected files in $($candidate.Path); leaving this candidate unchanged."
            continue
        }
        $Losers.Add([pscustomobject]@{
            Model = $candidate.Model
            Stage = $candidate.Stage
            Chi2 = $candidate.Chi2
            Path = $candidate.Path
            WinnerStage = $winner.Stage
            WinnerChi2 = $winner.Chi2
        })
    }
}

Write-Host "Fit candidate audit for $Scenario"
Write-Host "Valid model candidates: $($Candidates.Count); models: $($Winners.Count)"
Write-Host "Mode: $(if ($Apply) { 'APPLY (move to archive)' } else { 'PREVIEW ONLY' })"
foreach ($stage in $Stages) {
    $count = @($Winners | Where-Object Stage -eq $stage).Count
    Write-Host "  Selected from $stage : $count"
}
Write-Host "Superseded model directories eligible for archive: $($Losers.Count)"
foreach ($loser in $Losers) {
    Write-Host ('  [ARCHIVE] {0} / {1}: chi2={2:G7} -> selected {3}, chi2={4:G7}' -f `
        $loser.Model, $loser.Stage, $loser.Chi2, $loser.WinnerStage, $loser.WinnerChi2)
}
if (-not $Apply) {
    Write-Host 'Preview complete; nothing moved. Run with -Apply after checking the selected fits.'
    return
}
if ($Losers.Count -eq 0) {
    Write-Host 'Nothing superseded remains. No changes made.'
    return
}

$stamp = (Get-Date).ToUniversalTime().ToString('yyyyMMddTHHmmssZ')
$ArchiveRoot = Join-Path $ProjectRoot ('.cleanup_quarantine\superseded_fit_candidates\' + $stamp)
$Moved = [System.Collections.Generic.List[object]]::new()
foreach ($loser in $Losers) {
    $archivePath = Join-Path $ArchiveRoot ("Reports\output\full\{0}\{1}\{2}" -f $loser.Stage, $Scenario, $loser.Model)
    if (Test-Path -LiteralPath $archivePath) {
        throw "Archive path already exists; refusing overwrite: $archivePath"
    }
    if (-not $PSCmdlet.ShouldProcess($loser.Path, "Move superseded fit to $archivePath")) { continue }
    $destinationParent = Split-Path -Path $archivePath -Parent
    New-Item -ItemType Directory -Path $destinationParent -Force | Out-Null
    Move-Item -LiteralPath $loser.Path -Destination $archivePath -ErrorAction Stop
    $Moved.Add([pscustomobject]@{
        Model = $loser.Model
        PreviousStage = $loser.Stage
        PreviousChi2 = $loser.Chi2
        SelectedStage = $loser.WinnerStage
        SelectedChi2 = $loser.WinnerChi2
        OriginalPath = $loser.Path
        ArchivePath = $archivePath
    })
    Write-Host "  [MOVED] $($loser.Stage) / $($loser.Model)"
}
if ($Moved.Count -gt 0) {
    $manifest = Join-Path $ArchiveRoot 'archive_manifest.csv'
    $Moved | Export-Csv -LiteralPath $manifest -NoTypeInformation -Encoding UTF8
    Write-Host "Moved $($Moved.Count) superseded candidates. Manifest: $manifest"
    Write-Host 'The best candidates and batch/retry JSON/CSV summaries remain in Reports/output/full.'
    Write-Host 'Historical previous_summary paths inside some fitting JSON may now refer to archived folders; use the archive manifest for provenance.'
} else {
    Write-Host 'No changes made.'
}
