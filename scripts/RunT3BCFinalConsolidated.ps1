$ErrorActionPreference = 'Stop'
$root = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
Set-Location -LiteralPath $root
$script = 'scripts\FinalizeT3BCLeadingBenchmark.experimental.py'
$dest = 'output\audits\t3_bc_consolidated_leading_reference_v23\summary.json'
if (-not (Test-Path -LiteralPath $script)) { throw "Missing $script; extract v23 ZIP at repository root" }
Write-Host 'T3_C_V23: standalone integral and evidence-chain self-tests'
& python -u -X utf8 $script --self-test
if ($LASTEXITCODE -ne 0) { throw 'v23 self-test failed' }
Write-Host 'T3_C_V23: consolidated B/C UV topology, direct/sequential, low-energy reference'
& python -u -X utf8 $script --root $root --output $dest --low-scale 100
if ($LASTEXITCODE -ne 0) { throw 'v23 consolidated B/C computation failed' }
Write-Host "T3_C_V23: done; $dest"
Write-Host 'T3_C_V23: frozen benchmark + SM-only low-energy RGE only; physical C5 NOT CERTIFIED'
