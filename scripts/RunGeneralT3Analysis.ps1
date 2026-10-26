# One-command analysis of PREVIOUSLY COMPLETED baseline/extended Matchete runs.
# This does not rerun expensive matching, change input parameters or change
# any production numerical/RGBeta files. Do not treat this as a full UV theory.
# Run from anywhere on Windows PowerShell 5.1+:
# powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3Analysis.ps1
param(
    [string]$OutputRoot = 'output\general_t3_trial',
    [switch]$RequireProof
)
$ErrorActionPreference = 'Stop'
$repo = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$python = (Get-Command 'python' -ErrorAction Stop).Source
$wolframscript = (Get-Command 'wolframscript' -ErrorAction Stop).Source
$dest = if ([System.IO.Path]::IsPathRooted($OutputRoot)) {
    $OutputRoot
} else { Join-Path $repo $OutputRoot }
if (-not (Test-Path -LiteralPath $dest -PathType Container)) {
    throw "Matching results not found at $dest. Run the Batch 4 16-model matching first."
}

function Invoke-AnalysisStage {
    param([string]$StageName, [string]$Executable, [string[]]$Arguments)
    Write-Host "`n=== $StageName ===" -ForegroundColor Cyan
    & $Executable @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw ("{0} returned exit code {1}" -f $StageName, $LASTEXITCODE)
    }
}

# Basic accepted-operator and text-level one-loop result reports.
Invoke-AnalysisStage -StageName 'Matchete operator counts' -Executable $python -Arguments @(
    (Join-Path $repo 'scripts\ReportGeneralScalarTrial.py'),
    '--root', $dest, '--all', '--require-complete',
    '--output', (Join-Path $repo 'output\general_scalar_audit_report.json')
)
Invoke-AnalysisStage -StageName 'Textual C5 signatures' -Executable $python -Arguments @(
    (Join-Path $repo 'scripts\InspectGeneralScalarMatching.py'),
    '--root', $dest, '--require-complete',
    '--output', (Join-Path $repo 'output\general_matching_symbolic_report.json')
)
Invoke-AnalysisStage -StageName 'RGBeta structural coverage' -Executable $python -Arguments @(
    (Join-Path $repo 'scripts\CheckGeneralT3RGBetaCoverage.py'),
    '--output', (Join-Path $repo 'output\general_rgbeta_coverage.json')
)

# FORMAL Mathematica comparison: zero new quartic couplings and compare C5.
# This uses only the existing exported symbolic expressions; no FrontEnd.
$proofPath = Join-Path $dest 'c5_algebraic_report.json'
Write-Host "`n=== Formal C5 zero-extra limit ===" -ForegroundColor Cyan
& $wolframscript '-file' (Join-Path $repo 'scripts\CompareGeneralT3C5.wl') $dest $proofPath
$proofCode = $LASTEXITCODE
Write-Host ("Formal Wolfram comparison exit code: {0}" -f $proofCode)

# Always emit the combined summary even if Wolfram cannot prove identity.
$args = @(
    (Join-Path $repo 'scripts\SummarizeGeneralT3Physics.py'),
    '--root', $dest,
    '--output', (Join-Path $repo 'output\general_t3_physics_summary.json')
)
if ($RequireProof) { $args += '--require-proof' }
Invoke-AnalysisStage -StageName 'Combined physics report' -Executable $python -Arguments $args
if ($proofCode -ne 0) { exit $proofCode }
exit 0
