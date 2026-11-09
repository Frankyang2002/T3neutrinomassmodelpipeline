param([switch]$VerifyOnly)
$ErrorActionPreference='Stop'
$root=Split-Path -Parent $PSScriptRoot
Set-Location -LiteralPath $root
$py=Join-Path $root '.venv\Scripts\python.exe'
if(-not (Test-Path -LiteralPath $py)){throw "Python virtual environment missing: $py"}
& $py -m unittest discover -s tests -p 'test_t3_parallel_a_bc_native_matchete_v25_v27.py' -v
if($LASTEXITCODE -ne 0){throw "v25-v27 regression tests failed with exit code $LASTEXITCODE"}
$flags=@('--project-root','.')
if($VerifyOnly){$flags+=('--verify-only')}else{$flags+=('--resume')}
& $py -m scripts.RunT3ParallelABCMatcheteThresholdBridgeV25V27 @flags
if($LASTEXITCODE -ne 0){throw "v25-v27 native-source evidence audit BLOCKED with exit code $LASTEXITCODE"}
Write-Host 'v25-v27 source audit complete. No Match[], no physical Wilson evolution.'
