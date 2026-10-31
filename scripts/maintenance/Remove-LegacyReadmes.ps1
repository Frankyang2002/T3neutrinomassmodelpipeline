# Save as scripts/maintenance/Remove-LegacyReadmes.ps1
# Preview by default. -Apply deletes matching README files after a timestamped backup.
[CmdletBinding()]
param(
    [string]$RepoRoot = (Get-Location).Path,
    [switch]$Apply,
    [switch]$IncludeRootReadme
)

$ErrorActionPreference = 'Stop'
$root = (Resolve-Path -LiteralPath $RepoRoot).Path.TrimEnd([IO.Path]::DirectorySeparatorChar)
if (-not (Test-Path -LiteralPath (Join-Path $root '.git'))) {
    throw "Not a Git working tree root: $root"
}
if (-not (Get-Command git -ErrorAction SilentlyContinue)) {
    throw 'Git is required for tracked-file verification.'
}
$excludeParts = @('.git', '.venv', 'venv', 'node_modules', '__pycache__', '.pytest_cache', 'output', 'Reports')
$readmes = @(Get-ChildItem -LiteralPath $root -File -Recurse -Force | Where-Object {
    $candidate = $_
    $parts = $candidate.FullName.Substring($root.Length).TrimStart('\','/') -split '[\\/]'
    $candidate.Name -match '(?i)^readme(?:[._-].*)?\.(?:md|markdown|txt|rst)$' -and
    -not @($parts | Where-Object { $excludeParts -contains $_ }).Count -and
    ($IncludeRootReadme -or $candidate.DirectoryName.TrimEnd('\','/') -ne $root)
} | Sort-Object FullName)
if ($readmes.Count -eq 0) { Write-Host 'No matching README files found.'; exit 0 }
Write-Host "Repository: $root"
Write-Host "README files selected: $($readmes.Count)"
$readmes | ForEach-Object { Write-Host ('  ' + $_.FullName.Substring($root.Length + 1)) }
if (-not $Apply) {
    Write-Host "PREVIEW ONLY. Add -Apply to back up and delete selected files."
    exit 0
}
$dirty = @(git -C $root status --porcelain)
if ($LASTEXITCODE -ne 0 -or $dirty.Count -gt 0) {
    throw 'Working tree must be clean before deletion. Commit/stash changes first.'
}
$backupDir = Join-Path $root '_readme_cleanup_backups'
New-Item -ItemType Directory -Path $backupDir -Force | Out-Null
$stamp = Get-Date -Format 'yyyyMMdd-HHmmss'
$zipPath = Join-Path $backupDir "readmes-$stamp.zip"
$stage = Join-Path ([IO.Path]::GetTempPath()) "t3-readmes-$stamp-$PID"
try {
    New-Item -ItemType Directory -Path $stage -Force | Out-Null
    foreach ($f in $readmes) {
        $relative = $f.FullName.Substring($root.Length + 1)
        $destination = Join-Path $stage $relative
        New-Item -ItemType Directory -Force -Path (Split-Path -Parent $destination) | Out-Null
        Copy-Item -LiteralPath $f.FullName -Destination $destination
    }
    Compress-Archive -Path (Join-Path $stage '*') -DestinationPath $zipPath -Force
    if (-not (Test-Path -LiteralPath $zipPath)) { throw 'Backup archive failed.' }
    foreach ($f in $readmes) { Remove-Item -LiteralPath $f.FullName -Force }
    Write-Host "Deleted $($readmes.Count) README files. Backup: $zipPath"
    Write-Host 'Review deletions: git status --short'
    Write-Host 'Restore tracked files: git restore -- "**/README*"'
} finally {
    if (Test-Path -LiteralPath $stage) { Remove-Item -LiteralPath $stage -Recurse -Force }
}
