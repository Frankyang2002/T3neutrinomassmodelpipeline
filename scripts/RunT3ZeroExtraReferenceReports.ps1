# Local-only. Uses saved all-16 zero-extra control, never modifies GitHub.
[CmdletBinding()]
param(
  [string]$RunId = 'zero_extra_all16_low_v1',
  [string]$WeinbergOutputId = 'weinberg_running_v1',
  [string]$ScanId = 'gamma_095_105_v1',
  [ValidateSet('boundary_zero','held_zero','all')]
  [string]$Mode = 'boundary_zero'
)
$ErrorActionPreference = 'Stop'
$root = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
Push-Location $root
try {
  $source = Join-Path $root "output/studies/general_scalar_extended/zero_extra_all16/$RunId/summary.json"
  if (-not (Test-Path -LiteralPath $source -PathType Leaf)) {
    throw "Missing saved zero-extra all16 summary: $source"
  }
  python -u -X utf8 -m scripts.GenerateT3ZeroExtraWeinbergRunning `
    --run-id $RunId --output-id $WeinbergOutputId --resume
  if ($LASTEXITCODE -ne 0) { throw 'Failed generating verified Weinberg running.' }

  python -u -X utf8 -m scripts.GenerateT3ZeroExtraGamma `
    --run-id $RunId --scan-id $ScanId --mode $Mode --resume
  if ($LASTEXITCODE -ne 0) { throw 'Failed generating narrow gamma sensitivity.' }

  python -X utf8 -m Numerical.plotting.BuildT3ZeroExtraReferenceReports `
    --run-id $RunId --weinberg-output-id $WeinbergOutputId --scan-id $ScanId --mode $Mode
  if ($LASTEXITCODE -ne 0) { throw 'Failed rendering reference interactive reports.' }

  $index = Join-Path $root "Reports/output/display_zero_extra_all16/$RunId/index.html"
  if (-not (Test-Path -LiteralPath $index -PathType Leaf)) {
    throw "Missing report index: $index"
  }
  Write-Host "Created reference reports: $index"
  if ($Mode -eq 'boundary_zero') {
    $open = Join-Path $root "Reports/output/display_zero_extra_all16/$RunId/boundary_zero/interactive_rge.html"
  } else {
    $open = $index
  }
  if (-not (Test-Path -LiteralPath $open -PathType Leaf)) {
    throw "Missing expected report: $open"
  }
  Start-Process -FilePath $open
}
finally { Pop-Location }
