# A(-4) Matchete backend update and physical scalar-threshold rerun.
# Full replacement for the v1 launcher; never changes repository physics or GitHub.
# All Wolfram expressions are executed from .wl files to preserve quoting on Windows.

[CmdletBinding()]
param(
    [string] $RunId = 'thesis_a4_m061_scalar_v2',
    [switch] $InstallMatchete,
    [string] $Python = 'python',
    [string] $WolframScript = 'wolframscript'
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
Set-Location -LiteralPath $project

$runner = Join-Path $project 'scripts\RunT3ThesisFinalOneBatchSafe.py'
$scalar = Join-Path $project 'Lagrangian\RunT3AMinus4FullGeneralPhysicalScalar.wl'
if (!(Test-Path -LiteralPath $runner -PathType Leaf) -or
    !(Test-Path -LiteralPath $scalar -PathType Leaf)) {
    throw 'The current local A(-4) Python runner or Wolfram source is missing.'
}
if ($RunId -notmatch '^[a-zA-Z0-9_-]+$') {
    throw 'RunId may contain only letters, digits, underscore and hyphen.'
}
$run = Join-Path $project ("output\parallel_t3\full_general_final_one_batch\" + $RunId)
if (Test-Path -LiteralPath $run) {
    throw "Refusing to overwrite an existing run: $run"
}

# This is deliberately a literal here-string. Do not pass Wolfram code through
# -code on Windows PowerShell: quotation marks can disappear before evaluation.
$versionProbe = @'
f = FindFile["Matchete`"];
If[!StringQ[f],
  Print["T3_MAT_PATH_MISSING"]; Exit[3]
];
paths = NestList[DirectoryName, f, 10];
root = SelectFirst[paths,
  FileExistsQ[FileNameJoin[{#, "version"}]] &, Missing["NotFound"]];
If[MissingQ[root],
  Print["T3_MAT_VERSION_FILE_MISSING;FOUND=" <> ToString[f]]; Exit[4]
];
Print["T3_MAT_VERSION=" <>
  StringTrim[Import[FileNameJoin[{root, "version"}], "Text"]]];
Print["T3_MAT_ROOT=" <> root];
'@

function Invoke-WolframFileCode {
    param([Parameter(Mandatory=$true)][string]$Source)
    $temporary = Join-Path ([System.IO.Path]::GetTempPath()) (
        't3_wolfram_' + [guid]::NewGuid().ToString('N') + '.wl')
    try {
        [System.IO.File]::WriteAllText(
            $temporary, $Source, (New-Object System.Text.UTF8Encoding($false)))
        $lines = @(& $WolframScript -file $temporary 2>&1)
        $exitCode = $LASTEXITCODE
        return [pscustomobject]@{
            ExitCode = $exitCode
            Lines = $lines
        }
    }
    finally {
        Remove-Item -LiteralPath $temporary -Force -ErrorAction SilentlyContinue
    }
}

function Get-InstalledMatchete {
    $probeResult = Invoke-WolframFileCode -Source $versionProbe
    $lines = @($probeResult.Lines)
    if ($probeResult.ExitCode -ne 0) {
        throw ("Matchete location probe failed (Wolfram exit " +
            $probeResult.ExitCode + '): ' + ($lines -join ' | '))
    }
    $ver = @($lines | Where-Object { "$_" -match '^T3_MAT_VERSION=' })
    $dir = @($lines | Where-Object { "$_" -match '^T3_MAT_ROOT=' })
    if ($ver.Count -lt 1 -or $dir.Count -lt 1) {
        throw ("Matchete version file could not be located: " + ($lines -join ' | '))
    }
    $versionString = ([string]$ver[-1]).Substring('T3_MAT_VERSION='.Length)
    $match = [regex]::Match($versionString, '(\d+)\.(\d+)\.(\d+)')
    if (!$match.Success) {
        throw "Unrecognised Matchete version: $versionString"
    }
    return [pscustomobject]@{
        Version = $versionString
        Parsed = [version]::new([int]$match.Groups[1].Value,
                                 [int]$match.Groups[2].Value,
                                 [int]$match.Groups[3].Value)
        Root = ([string]$dir[-1]).Substring('T3_MAT_ROOT='.Length)
    }
}

$minimum = [version]::new(0,6,1)
$installed = Get-InstalledMatchete
Write-Host ("Detected Matchete " + $installed.Version + " in " + $installed.Root)

if ($installed.Parsed -lt $minimum) {
    if (!$InstallMatchete) {
        throw ("Installed Matchete " + $installed.Version +
            " is older than 0.6.1. Use -InstallMatchete to install the " +
            "upstream distribution with an automatic local backup.")
    }

    # Download and validate official upstream distribution BEFORE touching the
    # installed package. This avoids the interactive ChoiceDialog in install.m.
    $scratch = Join-Path ([System.IO.Path]::GetTempPath()) (
        't3_matchete_upgrade_' + [guid]::NewGuid().ToString('N'))
    New-Item -ItemType Directory -Path $scratch -Force | Out-Null
    try {
        $zip = Join-Path $scratch 'upstream.zip'
        $unpacked = Join-Path $scratch 'unpacked'
        $url = 'https://gitlab.com/matchete/matchete/-/archive/master/matchete-master.zip'
        Write-Host "Downloading official Matchete upstream archive..."
        Invoke-WebRequest -Uri $url -OutFile $zip -UseBasicParsing
        Expand-Archive -LiteralPath $zip -DestinationPath $unpacked -Force
        $source = Get-ChildItem -LiteralPath $unpacked -Directory |
            Where-Object {
                (Test-Path -LiteralPath (Join-Path $_.FullName 'version')) -and
                (Test-Path -LiteralPath (Join-Path $_.FullName 'Kernel\init.m'))
            } | Select-Object -First 1
        if ($null -eq $source) {
            throw 'Downloaded Matchete archive lacks version or Kernel/init.m.'
        }
        $newVersionString = (Get-Content -LiteralPath (
            Join-Path $source.FullName 'version') -Raw).Trim()
        $m = [regex]::Match($newVersionString, '(\d+)\.(\d+)\.(\d+)')
        if (!$m.Success) {
            throw "Unrecognised downloaded Matchete version: $newVersionString"
        }
        $newVersion = [version]::new(
            [int]$m.Groups[1].Value, [int]$m.Groups[2].Value,
            [int]$m.Groups[3].Value)
        if ($newVersion -lt $minimum) {
            throw "Upstream archive is too old ($newVersionString); original install unchanged."
        }

        # Preserve a complete copy of the previously installed package.
        $backups = Join-Path $project 'output\local_matchete_backups'
        New-Item -ItemType Directory -Force -Path $backups | Out-Null
        $stamp = (Get-Date).ToUniversalTime().ToString('yyyyMMddTHHmmssZ')
        $backup = Join-Path $backups (
            'Matchete_' + $installed.Version + '_' + $stamp)
        Copy-Item -LiteralPath $installed.Root -Destination $backup -Recurse -Force
        Write-Host "Previous Matchete backed up at $backup"

        $replacementAttempted = $false
        try {
            $replacementAttempted = $true
            Remove-Item -LiteralPath $installed.Root -Recurse -Force
            Copy-Item -LiteralPath $source.FullName `
                -Destination $installed.Root -Recurse -Force
            $installed = Get-InstalledMatchete
            if ($installed.Parsed -lt $minimum) {
                throw "Wolfram still sees Matchete $($installed.Version) after replacement."
            }
        }
        catch {
            if ($replacementAttempted) {
                Write-Warning 'Upgrade failed. Attempting restoration of package backup.'
                if (Test-Path -LiteralPath $installed.Root) {
                    Remove-Item -LiteralPath $installed.Root -Recurse -Force
                }
                Copy-Item -LiteralPath $backup `
                    -Destination $installed.Root -Recurse -Force
            }
            throw
        }
        Write-Host ("Matchete upgraded successfully to " + $installed.Version)
    }
    finally {
        Remove-Item -LiteralPath $scratch -Recurse -Force -ErrorAction SilentlyContinue
    }
}
else {
    Write-Host "Matchete already meets the minimum version; no update required."
}

# Do not substitute an exploratory source for exact saved first-threshold EFT1.
$names = @(
    'T3_A4_EXPLORATORY_SOURCE_MATCH',
    'T3_A4_TRACE_DEFINECG',
    'T3_A4_PROBE_TWO_HEAVY_TRACE'
)
$previous = @{}
foreach ($name in $names) {
    $previous[$name] = [Environment]::GetEnvironmentVariable($name, 'Process')
}
try {
    foreach ($name in $names) {
        [Environment]::SetEnvironmentVariable($name, $null, 'Process')
    }
    & $Python -X utf8 -m scripts.RunT3ThesisFinalOneBatchSafe `
        --model A_-4 --run-id $RunId --debug
    $runExit = $LASTEXITCODE
}
finally {
    foreach ($name in $names) {
        [Environment]::SetEnvironmentVariable(
            $name, $previous[$name], 'Process')
    }
}

$statusPath = Join-Path $run 'A_-4\native_two_threshold_status.json'
if (Test-Path -LiteralPath $statusPath -PathType Leaf) {
    $status = Get-Content -LiteralPath $statusPath -Raw | ConvertFrom-Json
    Write-Host ("Native status: " + $status.status)
    Write-Host ("Native stage: " + $status.stage)
    Write-Host ("Native reason: " + $status.reason)
    if ($status.status -like '*PhysicalScalarSoftOnly*') {
        Write-Host ('Scalar-soft succeeded. F-hard inheritance, EFT1 Wilson RG ' +
            'and finite subtraction are still outstanding.')
    }
}
else {
    Write-Warning "No native scalar-stage status generated. Check $run"
}
Write-Host ("Matchete version used: " + $installed.Version)
Write-Host ("Output directory: " + $run)
# Nonzero is expected for an incomplete physical C5 calculation.
exit $runExit
