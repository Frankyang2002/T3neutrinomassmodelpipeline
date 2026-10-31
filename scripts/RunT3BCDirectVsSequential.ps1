$ErrorActionPreference = 'Stop'
$root = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
Set-Location -LiteralPath $root
$script = 'scripts\DirectVsSequentialOneLoop.experimental.py'
$output = 'output\audits\t3_bc_direct_vs_sequential_full_mass_v22\summary.json'
if (-not (Test-Path -LiteralPath $script)) { throw "Missing $script; extract the v22 ZIP at the repository root." }
Write-Host 'T3_C_V22: safe direct-Matchete InputForm parser and normalization self-tests'
& python $script --self-test
if ($LASTEXITCODE -ne 0) { throw 'v22 self-test failed' }
Write-Host 'T3_C_V22: direct simultaneous all-heavy vs sequential F/Y/X full-mass ordered A'
& python $script --root "$root" --output "$output"
if ($LASTEXITCODE -ne 0) {
    throw "v22 finite direct-vs-sequential comparison failed; inspect $output if present. No physical C5 is certified."
}
Write-Host "T3_C_V22: completed; $output"
Write-Host 'T3_C_V22: physical C5 remains uncertified; no full EFT running performed'
