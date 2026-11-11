# Archive historical numerical branches without changing any project source/output.
# Intended project path: scripts/ArchiveLegacyT3Pipeline.ps1
[CmdletBinding()]
param(
    [string]$Root = (Split-Path -Parent $PSScriptRoot),
    [switch]$IncludeLegacyOutputs
)
$ErrorActionPreference = 'Stop'
$Root = (Resolve-Path -LiteralPath $Root).Path
$archiveDir = Join-Path $Root 'archives'
New-Item -ItemType Directory -Path $archiveDir -Force | Out-Null
$stamp = Get-Date -Format 'yyyyMMdd_HHmmss'
$staging = Join-Path $env:TEMP ("t3_legacy_archive_$stamp")
$archivePath = Join-Path $archiveDir ("t3_legacy_numerical_$stamp.zip")
if (Test-Path -LiteralPath $staging) { throw "Staging directory already exists: $staging" }
if (Test-Path -LiteralPath $archivePath) { throw "Archive already exists: $archivePath" }

# Explicit allowlist: old/alternative numerical branches, NOT the selected coupled-IVP scripts.
$paths = @(
    'Numerical/parallel_t3',
    'scripts/RunT3ThesisTwoModelPilot.py',
    'scripts/RunT3StrictOneLoopBatch2.py',
    'scripts/RunT3FHardTriangleBatch1.py',
    'scripts/RunT3FinalRepro.py',
    'scripts/RunT3FinalRepro.ps1',
    'scripts/RunT3PoleClosureBatch7.py',
    'scripts/RunT3ReducedScalarSobolComparison.py',
    'README_T3_THESIS_TWO_MODEL_PILOT.md',
    'README_T3_FINAL_REPRO_BATCH8.md',
    'README_SOURCE_AUDITED_ALL16.md',
    'docs/IMPLEMENTATION_STATUS.md',
    'Numerical/benchmarks/fitted_reference'
)
if ($IncludeLegacyOutputs) {
    $paths += @(
        'output/parallel_t3',
        'output/batch13', 'output/batch14', 'output/batch15',
        'output/batch16', 'output/batch17'
    )
}
$copied = New-Object System.Collections.Generic.List[string]
$missing = New-Object System.Collections.Generic.List[string]
try {
    New-Item -ItemType Directory -Path $staging -Force | Out-Null
    foreach ($relative in $paths) {
        $source = Join-Path $Root $relative
        if (-not (Test-Path -LiteralPath $source)) { $missing.Add($relative); continue }
        $dest = Join-Path $staging $relative
        New-Item -ItemType Directory -Path (Split-Path -Parent $dest) -Force | Out-Null
        Copy-Item -LiteralPath $source -Destination $dest -Recurse -Force
        $copied.Add($relative)
    }
    if ($copied.Count -eq 0) { throw 'No legacy files found; refusing to create empty archive.' }
    $entries = @(Get-ChildItem -LiteralPath $staging -File -Recurse | ForEach-Object {
        [pscustomobject]@{
            path = $_.FullName.Substring($staging.Length + 1).Replace('\','/')
            bytes = $_.Length
            sha256 = (Get-FileHash -LiteralPath $_.FullName -Algorithm SHA256).Hash.ToLowerInvariant()
        }
    })
    $manifest = [ordered]@{
        archive_kind = 'historical_numerical_branches_read_only_snapshot'
        generated_at = (Get-Date).ToString('o')
        source_root = $Root
        source_head = (& git -C $Root rev-parse HEAD 2>$null)
        include_legacy_outputs = [bool]$IncludeLegacyOutputs
        selected_coupled_ivp_paths_intentionally_excluded = @(
            'scripts/RunT3SimultaneousEFT1A4.py',
            'scripts/RunT3SimultaneousEFT1All16.py',
            'scripts/RunT3SimultaneousEFT1MixedBC.py',
            'scripts/RunT3SimultaneousEFT1StageFreezeHeatmaps.py',
            'output/studies/general_scalar_extended/experimental_simultaneous_eft1_*'
        )
        included_roots = @($copied)
        missing_optional_roots = @($missing)
        file_count = $entries.Count
        files = $entries
        note = 'This archive preserves selected legacy code/results for reference; it is not guaranteed to be independently executable without shared repo modules and environments. Nothing was moved or deleted.'
    }
    ($manifest | ConvertTo-Json -Depth 7) | Set-Content -LiteralPath (Join-Path $staging 'ARCHIVE_MANIFEST.json') -Encoding UTF8
    Add-Type -AssemblyName System.IO.Compression.FileSystem
    [System.IO.Compression.ZipFile]::CreateFromDirectory($staging, $archivePath, [System.IO.Compression.CompressionLevel]::Optimal, $false)
    $hash = (Get-FileHash -LiteralPath $archivePath -Algorithm SHA256).Hash
    Write-Host "CREATED: $archivePath"
    Write-Host "FILES: $($entries.Count) | SHA256: $hash"
    Write-Host 'Historical files were copied. Original paths and coupled-IVP outputs were not modified.'
    if ($missing.Count -gt 0) { Write-Host "Missing optional paths: $($missing -join ', ')" }
}
finally {
    if (Test-Path -LiteralPath $staging) { Remove-Item -LiteralPath $staging -Recurse -Force }
}
