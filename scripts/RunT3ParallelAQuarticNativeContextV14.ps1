# Source-conditioned native RGBeta PackageScope Lam + quartic projector audit.
# Run on Windows with installed Wolfram RGBeta v1.2.0; preserves v10/v13 evidence.
# A completed but BLOCKED native witness exits 7 after printing failures.
param([switch]$SmokeOnly)
$ErrorActionPreference='Stop'
$root=Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
Set-Location $root
$python=Join-Path $root '.venv\Scripts\python.exe'
if (-not (Test-Path -LiteralPath $python)) { throw "Missing .venv Python $python" }
function Assert-Exit([string]$step) { if($LASTEXITCODE -ne 0) { throw "$step failed with exit code $LASTEXITCODE" } }
& $python -m unittest discover -s tests -p 'test_t3_parallel_a_native_context_v14.py' -v
Assert-Exit 'v14 Python contract tests'
$base='output\parallel_t3\wilson'
$wl='scripts\RunT3ParallelANativeQuarticProjectorsV14.wl'
$smoke=Join-Path $base 'rgbeta_quartic_projectors_native_v14_smoke.json'
$full=Join-Path $base 'rgbeta_quartic_projectors_native_v14.json'
$flags=@('T3_EXPERIMENTAL_GENERAL_FERMIONS','T3_EXPERIMENTAL_GENERAL_QUADRATICS',
         'T3_EXPERIMENTAL_GENERAL_QUARTICS','T3_EXPERIMENTAL_FLAVOR_MASSES')
$previous=@{}
foreach($f in $flags) { $previous[$f]=[Environment]::GetEnvironmentVariable($f,'Process') }
try {
    foreach($f in $flags) { [Environment]::SetEnvironmentVariable($f,$null,'Process') }
    if(-not (Get-Command wolframscript -ErrorAction SilentlyContinue)) { throw 'Missing installed wolframscript' }
    if(-not (Test-Path -LiteralPath $smoke)) {
        & wolframscript -file $wl $smoke 'SMOKE'
        $code=$LASTEXITCODE
        if($code -ne 0 -and $code -ne 7) { throw "native v14 smoke could not produce evidence (exit $code)" }
    } else { Write-Host "Preserving existing native v14 smoke: $smoke" }
    & $python -m scripts.RunT3ParallelAQuarticNativeContextAuditV14 --project-root . --native $smoke --inspect
    Assert-Exit 'v14 source-guarded smoke inspection'
    & $python -m scripts.RunT3ParallelANativeComplexChannelsV15 --project-root . --native $smoke --output (Join-Path $base 'rgbeta_complex_conjugate_native_v15_smoke.json')
    Assert-Exit 'v15 smoke native complex-channel classification'
    $smokeStatus=(Get-Content -LiteralPath $smoke -Raw | ConvertFrom-Json).status
    if($smokeStatus -ne 'NativeRGBetaQuarticProjectorSmoke3ModelsPassed') {
        Write-Host 'v14 smoke is blocked. Old v10/v13 evidence unchanged. No full production audit.'
        exit 7
    }
    if($SmokeOnly) { Write-Host 'v14 smoke passed. Full run skipped by -SmokeOnly.';exit 0 }
    if(-not (Test-Path -LiteralPath $full)) {
        & wolframscript -file $wl $full
        $code=$LASTEXITCODE
        if($code -ne 0 -and $code -ne 7) { throw "native v14 full exporter failed before auditable evidence (exit $code)" }
    } else { Write-Host "Preserving existing native v14 full evidence: $full" }
    & $python -m scripts.RunT3ParallelAQuarticNativeContextAuditV14 --project-root . --native $full --inspect
    Assert-Exit 'v14 source-guarded full inspection'
    & $python -m scripts.RunT3ParallelANativeComplexChannelsV15 --project-root . --native $full --output (Join-Path $base 'rgbeta_complex_conjugate_native_v15.json')
    Assert-Exit 'v15 all-model native complex-channel classification'
    $fullStatus=(Get-Content -LiteralPath $full -Raw | ConvertFrom-Json).status
    if($fullStatus -ne 'NativeRGBetaRealQuarticLamProjectorResponsesObserved_Conditional') {
        Write-Host 'v14 full native responses remain blocked; see above per-coupling evidence.'
        exit 7
    }
    & $python -m scripts.RunT3ParallelAQuarticNativeContextAuditV14 --project-root . --native $full --audit
    Assert-Exit 'v14 full 16-model audit and v9 source recheck'
    Write-Host 'v14 real native responses verified conditionally. No physical Wilson run.'
}
finally {
    foreach($f in $flags) { [Environment]::SetEnvironmentVariable($f,$previous[$f],'Process') }
    Write-Host 'Native diagnostic environment variables restored.'
}
