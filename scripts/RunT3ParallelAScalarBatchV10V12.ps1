# Workstream A v10-v12: resumable, read-only source-guarded diagnostic batch.
# -NativeOnly: skip already passed v11/v12 math tests; verify existing evidence.
# -SkipNative: just run/check Python diagnostics; do not call Wolfram.
# Compatible with powershell.exe 5.1. No production physics/physical Wilson run.
param([switch]$SkipNative, [switch]$NativeOnly)
$ErrorActionPreference = 'Stop'
if ($SkipNative -and $NativeOnly) { throw '-SkipNative and -NativeOnly cannot be combined.' }
$root = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
Set-Location $root
$python = Join-Path $root '.venv\Scripts\python.exe'
if (-not (Test-Path -LiteralPath $python)) { throw "Missing local venv Python: $python" }

function Test-ExitCode([string]$step) {
    if ($LASTEXITCODE -ne 0) { throw "$step failed with exit code $LASTEXITCODE" }
}
function Verify-Existing([string]$which) {
    & $python -m scripts.RunT3ParallelAScalarBatchCheckV10V12 --project-root . --which $which
    Test-ExitCode "saved $which source/physics evidence verification"
}
function Run-If-Missing {
    param([string]$Path, [string]$Which, [string]$Executable, [string[]]$ArgumentList)
    if (Test-Path -LiteralPath $Path) {
        Verify-Existing $Which
        return
    }
    & $Executable @ArgumentList
    Test-ExitCode "$Which fresh diagnostic"
    Verify-Existing $Which
}

& $python -m unittest discover -s tests -p 'test_t3_parallel_a_native_quartic_projectors_v10.py' -v
Test-ExitCode 'v10 unit tests'
& $python -m unittest discover -s tests -p 'test_t3_parallel_a_scalar_batch_resume_v10_v12.py' -v
Test-ExitCode 'v10-v12 resume unit tests'
if (-not $NativeOnly) {
    & $python -m unittest discover -s tests -p 'test_t3_parallel_a_scalar_nonportal_v11.py' -v
    Test-ExitCode 'v11 unit tests'
    & $python -m unittest discover -s tests -p 'test_t3_parallel_a_scalar_topology_v12.py' -v
    Test-ExitCode 'v12 unit tests'
}
$base = 'output\parallel_t3\wilson'
$v11 = Join-Path $base 'nonportal_scalar_invariants_conditional_v11.json'
$v12 = Join-Path $base 'topology_scalar_action_conditional_v12.json'
if ($NativeOnly) {
    # Requires existing v11/v12 outputs: validates them without re-generating.
    Verify-Existing 'v11'
    Verify-Existing 'v12'
} else {
    Run-If-Missing -Path $v11 -Which 'v11' -Executable $python -ArgumentList @('-m','scripts.RunT3ParallelAScalarNonportalV11','--project-root','.')
    Run-If-Missing -Path $v12 -Which 'v12' -Executable $python -ArgumentList @('-m','scripts.RunT3ParallelATopologyScalarV12','--project-root','.')
}
if ($SkipNative) {
    Write-Host 'v11/v12 independently verified. Native RGBeta Lam not evaluated.'
    exit 0
}

# The native v10 exporter intentionally rejects these flags if enabled.
# This diagnostic runner isolates a baseline EFT1 Wolfram child from the
# caller's optional-extension settings. No flags are changed persistently:
# values are restored in finally. Do not remove the guard in the .wl file.
$flags = @(
    'T3_EXPERIMENTAL_GENERAL_FERMIONS',
    'T3_EXPERIMENTAL_GENERAL_QUADRATICS',
    'T3_EXPERIMENTAL_GENERAL_QUARTICS',
    'T3_EXPERIMENTAL_FLAVOR_MASSES'
)
$old = @{}
foreach ($flag in $flags) {
    $old[$flag] = [Environment]::GetEnvironmentVariable($flag, 'Process')
}
$enabled = @($flags | Where-Object { $old[$_] -eq '1' })
if ($enabled.Count -gt 0) {
    Write-Host ('v10 baseline EFT1 diagnostic: isolating enabled optional flags: ' + ($enabled -join ', '))
}
try {
    foreach ($flag in $flags) {
        [Environment]::SetEnvironmentVariable($flag, $null, 'Process')
    }
    foreach ($flag in $flags) {
        if ([Environment]::GetEnvironmentVariable($flag, 'Process') -eq '1') {
            throw "v10 baseline EFT1 flags could not be cleared in diagnostic process: $flag"
        }
    }
    if (-not (Get-Command 'wolframscript' -ErrorAction SilentlyContinue)) {
        throw 'wolframscript unavailable. v10 native execution not possible.'
    }
    $wl = 'scripts\RunT3ParallelANativeQuarticProjectorsV10.wl'
    $smoke = Join-Path $base 'rgbeta_quartic_projectors_native_v10_smoke.json'
    $full = Join-Path $base 'rgbeta_quartic_projectors_native_v10.json'
    $audit = Join-Path $base 'native_quartic_projectors_conditional_v10.json'
    # Previously failed v10 smoke is immutable diagnostic evidence, not a
    # passed witness. Preserve it and extract the actual per-coupling failures.
    # A fresh blocked smoke also writes JSON before exiting 7; inspect it.
    if (-not (Test-Path -LiteralPath $smoke)) {
        & wolframscript -file $wl $smoke 'SMOKE'
        $nativeExit = $LASTEXITCODE
        if ($nativeExit -ne 0 -and $nativeExit -ne 7) {
            throw "v10-smoke fresh native execution failed before auditable output (exit $nativeExit)"
        }
    }
    $smokeDoc = Get-Content -LiteralPath $smoke -Raw | ConvertFrom-Json
    if ($smokeDoc.status -eq 'Blocked_NativeRGBetaQuarticProjectorResponses') {
        Write-Host 'v10 smoke has recorded real native contractions but did NOT pass.'
        Write-Host 'Running source-guarded v13 per-coupling failure classification.'
        & $python -m scripts.RunT3ParallelAQuarticForensicsV13 --project-root . --smoke $smoke
        Test-ExitCode 'v13 native failure forensic classification'
        $digest = (Get-FileHash -LiteralPath $smoke -Algorithm SHA256).Hash.ToLowerInvariant().Substring(0,12)
        $forensic = Join-Path $base ("rgbeta_quartic_forensics_v13_{0}.json" -f $digest)
        $reprobe = Join-Path $base ("rgbeta_quartic_reprobe_v13_{0}.json" -f $digest)
        if (-not (Test-Path -LiteralPath $reprobe)) {
            & wolframscript -file 'scripts\RunT3ParallelAProbeQuarticResponseV13.wl' $forensic $reprobe
            Test-ExitCode 'v13 focused installed-RGBeta response re-probe'
        } else {
            Write-Host "v13 existing re-probe retained (not overwritten): $reprobe"
        }
        & $python -m scripts.RunT3ParallelACompareQuarticReprobeV13 --forensic $forensic --reprobe $reprobe
        Test-ExitCode 'v13 differential diagnostic classification'
        Write-Host "v13 forensic: $forensic"
        Write-Host "v13 differential native response: $reprobe"
        Write-Host 'v10 remains BLOCKED. Full v10 audit/Wilson physics NOT run.'
        exit 7
    }
    if ($smokeDoc.status -ne 'NativeRGBetaQuarticProjectorSmoke3ModelsPassed') {
        throw "v10-smoke unknown or inconsistent native status: $($smokeDoc.status)"
    }
    Verify-Existing 'v10-smoke'
    Run-If-Missing -Path $full -Which 'v10-native' -Executable 'wolframscript' -ArgumentList @('-file',$wl,$full)
    Run-If-Missing -Path $audit -Which 'v10-final' -Executable $python -ArgumentList @('-m','scripts.RunT3ParallelANativeQuarticProjectorsAuditV10','--project-root','.')
    & $python -m unittest discover -s tests -p 'test_t3_parallel_a_native_quartic_projectors_v10.py' -v
    Test-ExitCode 'v10 completed local evidence integration'
    Write-Host 'v10-v12 conditional scalar batch complete; no physical Wilson evolution authorized.'
}
finally {
    foreach ($flag in $flags) {
        [Environment]::SetEnvironmentVariable($flag, $old[$flag], 'Process')
    }
    Write-Host 'v10 diagnostic environment flags restored to their original process values.'
}
