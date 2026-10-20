<#
.SYNOPSIS
    Preview or remove regenerable and superseded Reports/output artifacts.

.DESCRIPTION
    Honours T3 project cleanup, audited against GitHub main at
    b75a1408307fe738b0dca25412011788b836eec0 (2026-10-20).

    By default PREVIEWS only. -Apply deletes matching local files.
    No GitHub commands are run; no source code, benchmark configurations,
    output/full raw artifacts, or numerical JSON/CSV in retained studies
    are deleted.

    Scope:
      - Superseded four-scenario/four-coupling benchmark reports;
      - Superseded E-p0-based shared-Yukawa optimisation report chain;
      - Historical symbolic smoke/dimension/hypercharge report outputs;
      - All generated PNG, HTML, and ZIP files under Reports/output;
      - Figure-folder README.txt files and stale generated report manifests;
      - Any remaining per-model _scenario_variants folders.

    By default, canonical Reports/output/models PDFs and TeX, combined
    analytical PDFs and TeX, and retained numerical JSON/CSV stay intact.
    -IncludeCompiledPdfs ALSO deletes all PDFs under Reports/output.

.EXAMPLE
    .\scripts\CleanReportOutput.ps1
.EXAMPLE
    .\scripts\CleanReportOutput.ps1 -ListFiles
.EXAMPLE
    .\scripts\CleanReportOutput.ps1 -Apply
.EXAMPLE
    .\scripts\CleanReportOutput.ps1 -Apply -IncludeCompiledPdfs
#>
[CmdletBinding()]
param(
    [switch] $Apply,
    [switch] $IncludeCompiledPdfs,
    [switch] $ListFiles
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

$projectRoot = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
$reportRoot = [System.IO.Path]::GetFullPath((Join-Path $projectRoot 'Reports/output'))
$reportPrefix = $reportRoot + [System.IO.Path]::DirectorySeparatorChar

if (-not (Test-Path -LiteralPath (Join-Path $projectRoot 'pipeline.py') -PathType Leaf)) {
    throw "Expected project root containing pipeline.py: $projectRoot. Place this file in scripts/."
}
if (-not (Test-Path -LiteralPath $reportRoot -PathType Container)) {
    throw "Reports/output not found: $reportRoot"
}
if (((Get-Item -LiteralPath $reportRoot).Attributes -band [System.IO.FileAttributes]::ReparsePoint) -ne 0) {
    throw 'Refusing to operate on a symbolic-link Reports/output directory.'
}

$targetFiles = [System.Collections.Generic.Dictionary[string, object]]::new(
    [System.StringComparer]::OrdinalIgnoreCase
)

function Add-TargetFile {
    param(
        [Parameter(Mandatory)] [System.IO.FileInfo] $File,
        [Parameter(Mandatory)] [string] $Reason
    )
    $path = [System.IO.Path]::GetFullPath($File.FullName)
    if (-not $path.StartsWith($script:reportPrefix, [System.StringComparison]::OrdinalIgnoreCase)) {
        throw "Outside Reports/output: $path"
    }
    if (($File.Attributes -band [System.IO.FileAttributes]::ReparsePoint) -ne 0) {
        Write-Warning "Skipping symbolic-link file: $path"
        return
    }
    if (-not $script:targetFiles.ContainsKey($path)) {
        $relative = $path.Substring($script:projectRoot.Length + 1).Replace('\', '/')
        $script:targetFiles.Add($path, [pscustomobject]@{
            Path         = $path
            RelativePath = $relative
            Bytes        = [long] $File.Length
            Reason       = $Reason
        })
    }
}

function Add-TargetDirectory {
    param(
        [Parameter(Mandatory)] [string] $RelativePath,
        [Parameter(Mandatory)] [string] $Reason
    )
    $path = [System.IO.Path]::GetFullPath((Join-Path $script:reportRoot $RelativePath))
    if (-not $path.StartsWith($script:reportPrefix, [System.StringComparison]::OrdinalIgnoreCase)) {
        throw "Directory outside Reports/output: $path"
    }
    if (-not (Test-Path -LiteralPath $path -PathType Container)) { return }
    if (((Get-Item -LiteralPath $path).Attributes -band [System.IO.FileAttributes]::ReparsePoint) -ne 0) {
        throw "Refusing to traverse symbolic-link directory: $path"
    }
    foreach ($file in Get-ChildItem -LiteralPath $path -Recurse -File -Force) {
        Add-TargetFile -File $file -Reason $Reason
    }
}

# Full retired studies; none are input sources for the current A-m2 common
# RGE output or per-model fitted-RGE report set. Existing legacy Python
# entrypoints may still reference these files and should not be run unchanged.
$retired = @(
    'full/beta_dominance',
    'full/observable_sensitivity',
    'full/four_benchmark_correlations',
    'full/common_parameter_comparison',
    'full/common_parameter_comparison_v2',
    'full/common_scalar_pilot',
    'full/common_yukawa_optimization',
    'full/common_yukawa_least_squares',
    'full/common_yukawa_final_summary',
    'full/common_yukawa_correlation_scans',
    'full/common_yukawa_three_panel_scans',
    'full/scalar_coupling_correlations',
    'full/benchmark_summary',
    'full/yukawa_grouping',
    'full/normalized_fitted_sensitivity',
    'full/comparison',
    'full/within_model_comparison',
    'full/hypercharge',
    'full/dimensions',
    'full/smoke',
    'display'
)
foreach ($relative in $retired) {
    Add-TargetDirectory -RelativePath $relative -Reason 'Retired study or generated display'
}

# Older per-model symbolic variants, if still present in the local checkout.
$modelsRoot = Join-Path $reportRoot 'models'
if (Test-Path -LiteralPath $modelsRoot -PathType Container) {
    foreach ($model in Get-ChildItem -LiteralPath $modelsRoot -Directory -Force) {
        Add-TargetDirectory -RelativePath "models/$($model.Name)/_scenario_variants" `
            -Reason 'Obsolete four-scenario PDF variants'
    }
}

# Metadata and outputs whose old paths/lineage no longer describe the current
# canonical report organisation.
$staleFiles = @(
    'full/numerical_resume_summary.csv',
    'full/numerical_resume_summary.json',
    'models/report_migration_manifest.json',
    'analytical/combined_analytical_manifest.json'
)
foreach ($relative in $staleFiles) {
    $path = Join-Path $reportRoot $relative
    if (Test-Path -LiteralPath $path -PathType Leaf) {
        Add-TargetFile -File (Get-Item -LiteralPath $path) -Reason 'Stale derived summary or manifest'
    }
}

# Aggressive presentation cleanup, but without discarding saved numerical
# JSON/CSV or the canonical symbolic source TeX files.
foreach ($file in Get-ChildItem -LiteralPath $reportRoot -File -Recurse -Force) {
    $ext = $file.Extension.ToLowerInvariant()
    if ($ext -in @('.png', '.html', '.zip')) {
        Add-TargetFile -File $file -Reason 'Regenerable rendered output'
    }
    elseif ($file.Name -eq 'README.txt' -and $file.Directory.Name -eq 'figures') {
        Add-TargetFile -File $file -Reason 'Generated plot README'
    }
    elseif ($IncludeCompiledPdfs -and $ext -eq '.pdf') {
        Add-TargetFile -File $file -Reason 'Regenerable compiled PDF (optional)'
    }
}

$items = @($targetFiles.Values | Sort-Object RelativePath)
$totalBytes = [long] 0
foreach ($item in $items) { $totalBytes += $item.Bytes }
$remainingBytes = [long] 0
$remainingFiles = 0
foreach ($file in Get-ChildItem -LiteralPath $reportRoot -File -Recurse -Force) {
    if (-not $targetFiles.ContainsKey([System.IO.Path]::GetFullPath($file.FullName))) {
        $remainingBytes += [long] $file.Length
        $remainingFiles++
    }
}

Write-Host "Project: $projectRoot"
Write-Host "Mode: $(if ($Apply) { 'APPLY - DELETE' } else { 'PREVIEW ONLY' })"
Write-Host "Candidate files: $($items.Count)"
Write-Host ('Candidate size: {0:N2} MB' -f ($totalBytes / 1MB))
Write-Host "Remaining files: $remainingFiles"
Write-Host ('Remaining size: {0:N2} MB' -f ($remainingBytes / 1MB))
Write-Host ''
if ($items.Count -gt 0) {
    $items | Group-Object Reason | ForEach-Object {
        $bytes = [long] 0
        foreach ($row in $_.Group) { $bytes += $row.Bytes }
        [pscustomobject]@{ Category = $_.Name; Files = $_.Count; MB = [math]::Round($bytes / 1MB, 2) }
    } | Format-Table -AutoSize | Out-String | Write-Host
}
if ($ListFiles) {
    $items | Select-Object RelativePath, Bytes, Reason | Format-Table -AutoSize -Wrap |
        Out-String -Width 220 | Write-Host
}
if (-not $Apply) {
    Write-Host 'Dry run complete. No files changed. To delete: .\scripts\CleanReportOutput.ps1 -Apply'
    Write-Host 'Optional: add -IncludeCompiledPdfs to remove even canonical compiled PDFs (TeX stays).'
    return
}

# Record the exact pre-deletion file selection outside Reports/output.
$manifestDir = Join-Path $projectRoot '.cleanup_quarantine'
New-Item -ItemType Directory -Path $manifestDir -Force | Out-Null
$stamp = (Get-Date).ToString('yyyyMMdd_HHmmss')
$manifestPath = Join-Path $manifestDir "report_cleanup_$stamp.csv"
$items | Select-Object RelativePath, Bytes, Reason |
    Export-Csv -LiteralPath $manifestPath -NoTypeInformation -Encoding UTF8
Write-Host "Pre-deletion manifest: $manifestPath"

$deleted = 0
foreach ($item in $items) {
    if (Test-Path -LiteralPath $item.Path -PathType Leaf) {
        Remove-Item -LiteralPath $item.Path -Force -ErrorAction Stop
        $deleted++
    }
}

# Prune empty directories only. A directory containing anything not selected
# is retained; no recursive directory deletion is performed.
$allDirs = @(Get-ChildItem -LiteralPath $reportRoot -Directory -Recurse -Force |
    Sort-Object { $_.FullName.Length } -Descending)
foreach ($directory in $allDirs) {
    if (-not (Test-Path -LiteralPath $directory.FullName -PathType Container)) { continue }
    if (($directory.Attributes -band [System.IO.FileAttributes]::ReparsePoint) -ne 0) { continue }
    $first = Get-ChildItem -LiteralPath $directory.FullName -Force | Select-Object -First 1
    if ($null -eq $first) {
        Remove-Item -LiteralPath $directory.FullName -Force -ErrorAction Stop
    }
}
Write-Host "Deleted $deleted report files. No files outside Reports/output were deleted."
Write-Host 'Only the cleanup manifest was created in .cleanup_quarantine.'
Write-Host 'Review tracked deletions with: git status --short'
