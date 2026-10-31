$ErrorActionPreference = 'Stop'
$root = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
Set-Location -LiteralPath $root
$py = 'python'
$log = 'output\studies\t3_full_stack_overnight_v1'
New-Item -ItemType Directory -Force -Path $log | Out-Null
Write-Host 'T3 overnight: Python tests and full model preflight'
& $py -u -X utf8 -m unittest discover -s tests -p 'test_t3_full_stack_overnight_v1.py' -v
if ($LASTEXITCODE -ne 0) { throw 'Full-stack tests failed' }
& $py -u -X utf8 -m scripts.RunT3FullStackOvernight --preflight
if ($LASTEXITCODE -ne 0) { throw 'Full-stack model preflight failed' }
Write-Host 'T3 overnight: all 16 matrix-valued UV RGEs, B/C five-mass trajectories, mass eigenbasis and frozen reference matching'
Write-Host 'Legacy Y/X trajectory runs are separated and are NOT the new B/C five-mass trajectory.'
& $py -u -X utf8 -m scripts.RunT3FullStackOvernight --all --resume `
  --output-id t3_full_stack_overnight_v1 --convergence --include-legacy-sequential
$code = $LASTEXITCODE
Write-Host "T3 overnight: output\studies\t3_full_stack_overnight_v1\summary.json ; exit=$code"
if ($code -ne 0) { throw "One or more T3 overnight stages blocked (exit $code); inspect summary.json and progress.json" }
