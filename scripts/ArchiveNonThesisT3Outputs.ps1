# Preview by default. This script NEVER modifies physics source code.
# Run from project root: .\scripts\ArchiveNonThesisT3Outputs.ps1
# To archive ONLY the listed newer output directories:
#   .\scripts\ArchiveNonThesisT3Outputs.ps1 -Apply
[CmdletBinding()]
param([switch]$Apply)
$ErrorActionPreference = 'Stop'
$root = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$chosen = Join-Path $root 'Reports/output/thesis_qualitative/representation_running'
if (-not (Test-Path -LiteralPath (Join-Path $chosen 'csv/frozen_stage_observable_changes.csv') -PathType Leaf)) {
    throw 'Original thesis dataset missing; refusing to archive anything.'
}
$oldDatasetPath = (Resolve-Path $chosen).Path
$archiveBase = Join-Path $root 'output/archive/non_authoritative_thesis_2026-11-11'
# Deliberately narrow: no background/, rough_inputs/, RGE/, Lagrangian/,
# Numerical/, scripts/, original Reports/output/ or any shared runtime files.
$candidates = @(
    'output/parallel_t3/thesis_A4_Dm2_pilot_v1',
    'output/parallel_t3/thesis_generalised_running_v1',
    'output/parallel_t3/thesis_original_style_figures_v2',
    'output/parallel_t3/thesis_original_style_figures_v3',
    'output/parallel_t3/thesis_running_BC_repaired_v1'
)
$items = @()
foreach ($rel in $candidates) {
    $source = Join-Path $root $rel
    if (-not (Test-Path -LiteralPath $source -PathType Container)) { continue }
    $destination = Join-Path $archiveBase $rel
    if (Test-Path -LiteralPath $destination) {
        throw "Refusing to overwrite existing archive path: $destination"
    }
    $items += [pscustomobject]@{ Source = $source; Destination = $destination; Relative = $rel }
}
if ($items.Count -eq 0) {
    Write-Host 'No listed newer output directories found. Nothing to archive.'
    return
}
Write-Host 'SELECTED original dataset (will not be moved):' $oldDatasetPath
$items | Format-Table Relative,Source,Destination -AutoSize
if (-not $Apply) {
    Write-Host 'PREVIEW ONLY. No files moved. Use -Apply to archive listed outputs.'
    return
}
foreach ($item in $items) {
    $parent = Split-Path -Parent $item.Destination
    New-Item -Path $parent -ItemType Directory -Force | Out-Null
    Move-Item -LiteralPath $item.Source -Destination $item.Destination -ErrorAction Stop
    Write-Host 'ARCHIVED:' $item.Relative
}
$record = [ordered]@{
    schema = 't3-non-authoritative-output-archive-v1'
    created_at_utc = (Get-Date).ToUniversalTime().ToString('o')
    authoritative_dataset = 'Reports/output/thesis_qualitative/representation_running'
    moved = @($items | ForEach-Object { $_.Relative })
    source_code_moved = $false
}
$record | ConvertTo-Json -Depth 4 |
    Set-Content -LiteralPath (Join-Path $archiveBase 'ARCHIVE_MANIFEST.json') -Encoding utf8
Write-Host 'Archive completed. No source code or original thesis outputs were moved.'
