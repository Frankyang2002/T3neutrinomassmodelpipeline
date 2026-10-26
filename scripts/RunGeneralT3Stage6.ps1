<#
  Formal analysis after Batch 5 with existing Matchete + RGBeta artifacts.
  No production physics file is changed or numerical C5 value fabricated.

  powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3Stage6.ps1
  powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3Stage6.ps1 -SkipTensor
  Note: refreshes the formal C5 comparison from existing matching artifacts.
  -SkipFormalRefresh is only for a rerun after the NEW proof is generated.
  -TensorOnly runs only the standalone Matchete tensor exporter.
#>
param(
    [string]$TrialRoot = 'output\general_t3_trial',
    [string]$RGBetaRoot = 'output',
    [switch]$SkipTensor,
    [switch]$SkipRGE,
    [switch]$SkipFormalRefresh,
    [switch]$TensorOnly
)
$ErrorActionPreference = 'Stop'
$repo = Split-Path -Parent $PSScriptRoot
$python = (Get-Command python -ErrorAction Stop).Source
$wolframscript = (Get-Command wolframscript -ErrorAction Stop).Source
$trial = if ([System.IO.Path]::IsPathRooted($TrialRoot)) { $TrialRoot } else { Join-Path $repo $TrialRoot }
$rgroot = if ([System.IO.Path]::IsPathRooted($RGBetaRoot)) { $RGBetaRoot } else { Join-Path $repo $RGBetaRoot }
$out = Join-Path $repo 'output'
if (-not (Test-Path (Join-Path $trial 'c5_algebraic_report.json'))) {
    throw "Missing Batch 5 formal C5 proof. Run scripts\RunGeneralT3Analysis.ps1 -RequireProof first."
}
function Invoke-Required {
    param([string]$Stage, [string]$Executable, [string[]]$Arguments)
    Write-Host "`n=== $Stage ===" -ForegroundColor Cyan
    & $Executable @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw ("{0} failed with exit code {1}. Inspect its JSON output before continuing." -f $Stage, $LASTEXITCODE)
    }
}
# Fast targeted retry after a tensor-export-only failure: do not re-run
# already proven formal C5 or RGBeta source checks.
if ($TensorOnly) {
    Invoke-Required 'Matchete doublet tensor convention' $wolframscript @(
        '-file', (Join-Path $repo 'scripts\ExportGeneralT3TensorConvention.wl'),
        (Join-Path $out 'general_t3_trial\matchete_active_quartic_tensors.json')
    )
    Write-Host "`nStandalone tensor validation completed. RGBeta normalisation remains unverified." -ForegroundColor Green
    exit 0
}
Invoke-Required 'Wolfram extra-coupling atom self-test' $wolframscript @(
    '-file', (Join-Path $repo 'scripts\CompareGeneralT3C5.wl'), 'SELFTEST'
)
# Refresh the formal response; a previously generated JSON from Batch 5 or
# Batch 6 does not contain the corrected atom extraction and proof fields.
if (-not $SkipFormalRefresh) {
    Invoke-Required 'Refresh formal C5 difference and response atoms' $wolframscript @(
        '-file', (Join-Path $repo 'scripts\CompareGeneralT3C5.wl'),
        $trial, (Join-Path $trial 'c5_algebraic_report.json')
    )
}
Invoke-Required 'Formal delta-C5 response and Yukawa families' $python @(
    (Join-Path $repo 'scripts\AssessGeneralT3C5Response.py'),
    '--root', $trial, '--output', (Join-Path $out 'general_c5_response.json'), '--strict'
)
if (-not $SkipRGE) {
    Invoke-Required 'One-loop RGBeta source term proof for B/C(-1)' $wolframscript @(
        '-file', (Join-Path $repo 'scripts\AnalyzeGeneralT3RGBetaSources.wl'),
        $rgroot, (Join-Path $out 'general_rgbeta_source_proof.json')
    )
}
if (-not $SkipTensor) {
    Invoke-Required 'Matchete doublet tensor convention' $wolframscript @(
        '-file', (Join-Path $repo 'scripts\ExportGeneralT3TensorConvention.wl'),
        (Join-Path $out 'general_t3_trial\matchete_active_quartic_tensors.json')
    )
}
Write-Host "`nStage 6 diagnostics completed. The RGBeta coupling NORMALISATIONS are not yet verified." -ForegroundColor Green
