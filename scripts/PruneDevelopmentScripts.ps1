# One-time, local-only cleanup of historical T3 development entry points.
# Default: preview. -Apply: make verified archive outside repository, then delete.
# No GitHub writes, branch inspection, HEAD commands, fetch, pull, or commits.
[CmdletBinding()]
param([switch] $Apply)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
$scriptsDir = Join-Path $repo 'scripts'

# Exact filename categories reviewed from GitHub repository contents.
# The expected inventory checksum below prevents expansion to new/unreviewed files.
$patterns = @(
    '(?:Worker|Workstream)',
    '\.experimental\.',
    '^Audit',
    '^(?:RunBatch(?:13|14|15|16|17|34)|CheckBatch34|RunT3.*Batch\d+)',
    '^(?:Archive|Retire|CleanReportOutput|RenameT3ModelFolders|Consolidate)',
    '\.before_v'
)

$candidates = @{}
foreach ($file in @(Get-ChildItem -LiteralPath $scriptsDir -File)) {
    foreach ($pattern in $patterns) {
        if ([regex]::IsMatch($file.Name, $pattern)) {
            $rel = 'scripts/' + $file.Name
            $candidates[$rel] = $file.FullName
            break
        }
    }
}

$names = [string[]]@($candidates.Keys)
[Array]::Sort($names, [StringComparer]::Ordinal)

# FNV-1a checksum of the newline-separated, ordinal-sorted path inventory.
# This is an inventory guard, not a security hash or content-certification claim.
$joined = [string]::Join("`n", $names)
[long]$checksum = 2166136261
foreach ($byte in [Text.Encoding]::UTF8.GetBytes($joined)) {
    $checksum = ((($checksum -bxor [long]$byte) * 16777619) % 4294967296)
}
$inventory = ([uint32]$checksum).ToString('x8')
if ($names.Count -ne 151 -or $inventory -cne '5e917598') {
    throw "The candidate inventory differs from the inspected repository (count=$($names.Count), checksum=$inventory; expected 151 / 5e917598). No files removed."
}

# Scan currently RETAINED executable source, configuration files and the permanent
# scientific documentation for explicit references. Deliberately do not treat old
# reports, historical notes, caches or output ledgers as executable dependencies.
$sourceDirs = @(
    'common', 'model', 'Lagrangian', 'RGE', 'Numerical', 'physics',
    'studies', 'validation', 'tests', 'configs', 'input', 'scripts',
    'tools', 'display', 'Reports'
)
$codeExtensions = @('.py', '.ps1', '.psm1', '.wl', '.wls', '.m', '.sh', '.json', '.toml', '.yaml', '.yml')
$sourceFiles = New-Object 'System.Collections.Generic.List[System.IO.FileInfo]'
foreach ($dir in $sourceDirs) {
    $path = Join-Path $repo $dir
    if (!(Test-Path -LiteralPath $path -PathType Container)) { continue }
    foreach ($file in @(Get-ChildItem -LiteralPath $path -Recurse -File)) {
        if ($codeExtensions -notcontains $file.Extension.ToLowerInvariant()) { continue }
        $rel = $file.FullName.Substring($repo.Length + 1).Replace('\', '/')
        if ($rel -match '(?:^|/)(?:__pycache__|\.venv|\.git|output|\.cleanup_quarantine)/') { continue }
        # Skip generated/archived results: retaining their provenance is a separate review.
        if ($rel -match '^(?:Reports|display)/output/') { continue }
        # A removable script is not an independent reason to keep another removable script.
        if ($candidates.ContainsKey($rel)) { continue }
        $sourceFiles.Add($file)
    }
}
foreach ($rel in @('pipeline.py', 'full16_symbolic.py', 'README.md',
                    'docs/README.md', 'docs/RESULTS.md', 'docs/NUMERICAL.md',
                    'docs/REPRODUCIBILITY.md', 'docs/MODEL_AND_PIPELINE.md',
                    'docs/INSTALLATION.md', 'docs/IMPLEMENTATION_STATUS.md',
                    'docs/THESIS_SCOPE.md')) {
    $full = Join-Path $repo ($rel.Replace('/', [IO.Path]::DirectorySeparatorChar))
    if (Test-Path -LiteralPath $full -PathType Leaf) {
        $sourceFiles.Add((Get-Item -LiteralPath $full))
    }
}

# Build reference tokens. For a backed-up .before_v* file, match the FULL file
# name only; the stem could otherwise collide with the legitimate live module.
$tokens = @{}
foreach ($rel in $names) {
    $filename = [IO.Path]::GetFileName($rel)
    if ($filename -match '\.before_v') {
        $tokens[$rel] = @($filename)
    }
    else {
        $stem = [IO.Path]::GetFileNameWithoutExtension($filename)
        $tokens[$rel] = @($filename, $stem)
    }
}

$protected = @{}
$protectedBy = @{}
$contents = @{}
foreach ($file in $sourceFiles) {
    $rel = $file.FullName.Substring($repo.Length + 1).Replace('\', '/')
    try {
        $contents[$rel] = [IO.File]::ReadAllText($file.FullName)
    }
    catch {
        throw "Cannot safely inspect retained source $rel : $_"
    }
}

function Test-ExplicitReference([string] $haystack, [string[]] $needles) {
    foreach ($needle in $needles) {
        if ($haystack.IndexOf($needle, [StringComparison]::Ordinal) -ge 0) {
            return $true
        }
    }
    return $false
}

# Protect candidate scripts directly referenced by a noncandidate, retained file.
foreach ($rel in $names) {
    foreach ($sourceRel in $contents.Keys) {
        if (Test-ExplicitReference $contents[$sourceRel] $tokens[$rel]) {
            $protected[$rel] = $true
            $protectedBy[$rel] = $sourceRel
            break
        }
    }
}

# Transitive closure: when a protected candidate script calls another candidate,
# preserve that dependency too. A mutually referencing obsolete group is removable.
$changed = $true
while ($changed) {
    $changed = $false
    foreach ($parent in @($protected.Keys)) {
        if (!$candidates.ContainsKey($parent)) { continue }
        $parentText = [IO.File]::ReadAllText($candidates[$parent])
        foreach ($child in $names) {
            if ($protected.ContainsKey($child) -or $child -ceq $parent) { continue }
            if (Test-ExplicitReference $parentText $tokens[$child]) {
                $protected[$child] = $true
                $protectedBy[$child] = $parent
                $changed = $true
            }
        }
    }
}

$removals = @($names | Where-Object { !$protected.ContainsKey($_) })
Write-Host "Repository: $repo"
Write-Host "Reviewed script candidates: $($names.Count)"
Write-Host "Protected by retained-source reference: $($protected.Count)"
Write-Host "Selected for archiving/removal: $($removals.Count)"
Write-Host 'Protected candidates (reason):'
foreach ($rel in @($protected.Keys | Sort-Object)) {
    Write-Host "  KEEP $rel <- $($protectedBy[$rel])"
}
Write-Host 'Planned removals:'
$removals | ForEach-Object { Write-Host "  REMOVE $_" }

if (!$Apply) {
    Write-Host 'PREVIEW ONLY; no files were changed.'
    Write-Host 'Run .\scripts\PruneDevelopmentScripts.ps1 -Apply to archive and delete only the unreferenced candidates.'
    return
}

if ($removals.Count -eq 0) {
    Write-Host 'Nothing eligible for removal.'
    return
}

# Never remove modified/untracked candidates. This inspects local worktree status
# only. It does not query GitHub, a remote branch, or a commit/HEAD reference.
$gitOutput = @(& git -C $repo status --porcelain --untracked-files=all -- @removals)
if ($LASTEXITCODE -ne 0) {
    throw 'Cannot read local candidate status. No files removed.'
}
if ($gitOutput.Count -gt 0) {
    Write-Host 'Local candidate changes:'
    $gitOutput | ForEach-Object { Write-Host $_ }
    throw 'Candidate scripts have uncommitted changes. No files removed.'
}

# Archive outside the repository; record SHA-256 to allow restoration.
$parent = Split-Path -Parent $repo
$stamp = Get-Date -Format 'yyyyMMdd_HHmmss'
$suffix = [Guid]::NewGuid().ToString('N').Substring(0, 8)
$archive = Join-Path $parent "T3_removed_development_scripts_${stamp}_$suffix"
New-Item -ItemType Directory -Path $archive -ErrorAction Stop | Out-Null
$records = New-Object 'System.Collections.Generic.List[object]'
foreach ($rel in $removals) {
    $source = Join-Path $repo ($rel.Replace('/', [IO.Path]::DirectorySeparatorChar))
    $destination = Join-Path $archive ($rel.Replace('/', [IO.Path]::DirectorySeparatorChar))
    $destDir = Split-Path -Parent $destination
    New-Item -ItemType Directory -Path $destDir -Force | Out-Null
    Copy-Item -LiteralPath $source -Destination $destination -ErrorAction Stop
    $hash1 = (Get-FileHash -LiteralPath $source -Algorithm SHA256).Hash
    $hash2 = (Get-FileHash -LiteralPath $destination -Algorithm SHA256).Hash
    if ($hash1 -cne $hash2) {
        throw "Archive hash mismatch for $rel. No files removed."
    }
    $records.Add([pscustomobject]@{ Path = $rel; SHA256 = $hash1 })
}
$records | Export-Csv -LiteralPath (Join-Path $archive 'removed_files.csv') -NoTypeInformation -Encoding UTF8
$protected.Keys | Sort-Object | ForEach-Object {
    [pscustomobject]@{ Path = $_; ReferencedBy = $protectedBy[$_] }
} | Export-Csv -LiteralPath (Join-Path $archive 'protected_files.csv') -NoTypeInformation -Encoding UTF8
([pscustomobject]@{
    CandidateInventory = $inventory
    CandidateCount = $names.Count
    ProtectedCount = $protected.Count
    RemovedCount = $removals.Count
    CreatedAt = (Get-Date).ToString('o')
    Repository = $repo
}) | ConvertTo-Json | Set-Content -LiteralPath (Join-Path $archive 'summary.json') -Encoding UTF8

# Only remove files backed up and SHA-256 verified above.
foreach ($rel in $removals) {
    Remove-Item -LiteralPath (Join-Path $repo ($rel.Replace('/', [IO.Path]::DirectorySeparatorChar))) -ErrorAction Stop
}
Write-Host "Removed $($removals.Count) historical script files. Backup: $archive"
Write-Host 'Verify: .\scripts\RunCoreTests.ps1; python -X utf8 -m pytest -q --tb=short; python -X utf8 pipeline.py --help'
Write-Host 'After verification, remove this one-time cleanup script and its README.'
