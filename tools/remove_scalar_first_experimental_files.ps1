$ErrorActionPreference = "Stop"

$required = @(
    ".\common\PipelineCLI.py",
    ".\common\PipelinePlan.py",
    ".\Lagrangian\Runner.py",
    ".\Lagrangian\T3ModelMatching.py",
    ".\RGE\running\IntermediateEFTRunning.py"
)

$missing = $required | Where-Object { -not (Test-Path $_ -PathType Leaf) }
if ($missing.Count -gt 0) {
    throw "Refusing cleanup because replacement files are missing: $($missing -join ', ')"
}

$retired = @(
    ".\RGE\general\Psi2Phi3RGE.py",
    ".\RGE\running\intermediate\ScalarFirstDimensionSixSeed.py",
    ".\tests\test_psi2phi3_rge.py",
    ".\tests\test_scalar_first_dimension_six_foundation.py"
)

foreach ($path in $retired) {
    if (Test-Path $path) {
        Remove-Item $path -Force
        Write-Host "Removed $path"
    }
}

Get-ChildItem -Recurse -Directory -Filter __pycache__ -ErrorAction SilentlyContinue |
    Remove-Item -Recurse -Force
Get-ChildItem -Recurse -Directory -Filter .pytest_cache -ErrorAction SilentlyContinue |
    Remove-Item -Recurse -Force

Write-Host "Scalar-first experimental source cleanup complete."
