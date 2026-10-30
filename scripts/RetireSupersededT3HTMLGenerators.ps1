<#!
.SYNOPSIS
  Preview or archive HTML-only general-scalar reporting entry points AFTER the
  new report has verified 64 scalar and sensitivity cases.
.DESCRIPTION
  - No PDF generator is included.
  - No numerical output, saved fits, or historical Lambda template is touched.
  - No action by default; -Apply moves sources to local .cleanup_quarantine.
  - All current .py imports must be checked before moving these files.
#>
param(
  [string]$RunId = 'frozen_general64_v1',
  [switch]$Apply
)
Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$coverage = Join-Path $project "Reports\output\display_general_scalar_frozen_64\$RunId\coverage.csv"
$targets = @(
  'Numerical\plotting\BuildFrozenGeneralScalarLambdaReports.py',
  'scripts\BuildT3FrozenLambdaNative.py',
  'scripts\BuildT3FrozenDisplaySummary.py'
)
Write-Host 'Candidate HTML-only legacy entry points:'
foreach ($rel in $targets) { Write-Host ('  ' + $rel) }
if (-not $Apply) {
  Write-Host 'PREVIEW ONLY. No files modified.'
  Write-Host 'A move requires validated 64/64 general-scalar scalar and gamma coverage.'
  return
}
if (-not (Test-Path -LiteralPath $coverage -PathType Leaf)) {
  throw "Missing validated new coverage: $coverage"
}
$rows = @(Import-Csv -LiteralPath $coverage)
if ($rows.Count -ne 64) {throw "Expected 64 coverage rows; got $($rows.Count)"}
foreach ($row in $rows) {
  if ($row.scalar_rge -ne 'True' -or $row.complete_gamma_scan -ne 'True') {
    throw "Do not retire generators until all 64 general-scalar RGE and gamma scans are validated."
  }
}
$names = @('BuildFrozenGeneralScalarLambdaReports','BuildT3FrozenLambdaNative','BuildT3FrozenDisplaySummary')
$outside = @()
$sourceFiles = @(Get-ChildItem -LiteralPath (Join-Path $project 'scripts'),(Join-Path $project 'Numerical') -Recurse -File -Filter '*.py')
foreach ($file in $sourceFiles) {
  $relative = [IO.Path]::GetRelativePath($project,$file.FullName)
  if (($targets -contains $relative)) {continue}
  $matches = @(Select-String -LiteralPath $file.FullName -Pattern $names -SimpleMatch)
  if ($matches.Count -gt 0) {$outside += $relative}
}
if ($outside.Count -gt 0) {
  throw ('Other Python scripts still reference these generators: ' + (($outside | Sort-Object -Unique) -join ', '))
}
$stamp = (Get-Date -Format 'yyyyMMdd_HHmmss')
$archive = Join-Path $project ".cleanup_quarantine\html_generators\$stamp"
foreach ($rel in $targets) {
  $source = Join-Path $project $rel
  if (-not (Test-Path -LiteralPath $source -PathType Leaf)) {continue}
  $destination = Join-Path $archive $rel
  New-Item -ItemType Directory -Force -Path (Split-Path -Parent $destination) | Out-Null
  Move-Item -LiteralPath $source -Destination $destination -ErrorAction Stop
  Write-Host "ARCHIVED: $rel"
}
Write-Host "Local archive: $archive"
