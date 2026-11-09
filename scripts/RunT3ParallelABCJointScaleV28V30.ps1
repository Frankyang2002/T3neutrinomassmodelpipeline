param(
 [ValidateSet('auto','five_mass','quadratic_only')][string]$Backend = 'auto',
 [ValidateSet('geometric','upper')][string]$Scheme = 'geometric',
 [switch]$SkipTests
)
$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $PSScriptRoot
Set-Location $root
$python = Join-Path $root '.venv\Scripts\python.exe'
if (-not (Test-Path -LiteralPath $python)) { $python='python' }
if (-not $SkipTests) {
 & $python -m unittest discover -s tests -p 'test_t3_parallel_a_bc_joint_scale_v28_v30.py' -v
 if ($LASTEXITCODE -ne 0) { throw "v28-v30 source guards failed ($LASTEXITCODE)" }
}
& $python -m scripts.RunT3ParallelABCJointScaleV28V30 --backend $Backend --scheme $Scheme
if ($LASTEXITCODE -ne 0) { throw "v28-v30 candidate EFT1 joint scale BLOCKED ($LASTEXITCODE). Inspect missing source or inconsistent saved beta export." }
Write-Host 'v28-v30 candidate renormalisable scale replay complete; physical Wilson/SMEFT matching NOT performed.'
