param(
    [switch]$Resume,
    [switch]$NoScan
)
$ErrorActionPreference='Stop'
Set-Location -LiteralPath (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
Write-Host 'T3 B/C running-to-source matching: regression tests (includes live Matchete v22 frozen parity when source is present)'
& python -u -X utf8 -m unittest discover -s tests -p 'test_t3_bc_running_one_scale_triangle_v1.py' -v
if ($LASTEXITCODE -ne 0) { throw 'B/C matching adapter tests failed' }
$flags = @('--all')
if (-not $NoScan) { $flags += '--scan' }
if ($Resume) { $flags += '--resume' }
Write-Host 'T3 B/C: common-scale running UV endpoint -> archived restricted Matchete direct 480-term A5'
& python -u -X utf8 -m scripts.RunT3BCRunningOneScaleTriangle @flags
if ($LASTEXITCODE -ne 0) { throw 'B/C running-to-matching campaign blocked; inspect output/studies/t3_bc_running_one_scale_triangle_v1/progress.json' }
Write-Host 'T3 B/C running one-scale reference completed; physical C5 remains NOT CERTIFIED.'
