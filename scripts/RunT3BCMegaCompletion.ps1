param(
    [switch]$Resume,
    [switch]$SkipComplexUV,
    [switch]$SkipSMProxy,
    [string]$OutputId = 't3_bc_mega_completion_v1'
)
$ErrorActionPreference = 'Stop'
$root = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
Set-Location -LiteralPath $root
Write-Host 'T3 B/C mega: complete standalone and source-backed contract tests'
& python -u -X utf8 -m unittest discover -s tests -p 'test_t3_bc_mega_completion_v1.py' -v
if ($LASTEXITCODE -ne 0) { throw 'Mega-batch regression tests failed' }
Write-Host 'T3 B/C mega: complex-heavy RGBeta UV covariance, same-scale source matching, native quartics, SM-only reference'
$argsPython = @('-u','-X','utf8','-m','scripts.RunT3BCMegaCompletion','--all','--output-id',$OutputId)
if ($Resume) {$argsPython += '--resume'}
if ($SkipComplexUV) {$argsPython += '--skip-complex-uv'}
if ($SkipSMProxy) {$argsPython += '--skip-sm-proxy'}
& python @argsPython
if ($LASTEXITCODE -ne 0) { throw 'Mega-batch numerical campaign has blocking results; inspect summary.json' }
Write-Host "T3 B/C mega COMPLETE: output/studies/$OutputId/summary.json (physical matching NOT certified)"
