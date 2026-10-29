<#
  Experimental two-branch numerical probes after Stage 9.
  Uses only existing Matchete C5 derivative and RGBeta source reports.
  NOT a combined end-to-end numerical result: Matchete/RGBeta coupling map
  and the full scalar mass-mixing calculation remain unverified.
#>
param(
    [string]$TrialRoot = 'output\general_t3_trial',
    [string]$Proof = 'output\general_rgbeta_source_proof.json',
    [string]$Config = 'docs\general_c5_toy_benchmark.json',
    [string]$OutputRoot = 'output\general_t3_numerical_pilot'
)
$ErrorActionPreference = 'Stop'
$repo = Split-Path -Parent $PSScriptRoot
$python = (Get-Command python -ErrorAction Stop).Source
function Resolve-RepoPath([string]$p) {
    if ([System.IO.Path]::IsPathRooted($p)) { return $p }
    return Join-Path $repo $p
}
$trial = Resolve-RepoPath $TrialRoot
$proof = Resolve-RepoPath $Proof
$configFile = Resolve-RepoPath $Config
$out = Resolve-RepoPath $OutputRoot
New-Item -ItemType Directory -Force -Path $out | Out-Null
$failures = @()
foreach ($model in @('B_m1', 'C_m1')) {
    Write-Host "`n=== $model Matchete C5 pilot ==="
    & $python (Join-Path $repo 'scripts\RunGeneralT3NumericalPilot.py') --root $trial --model $model --input $configFile --output (Join-Path $out ($model + '_matchete_hard_c5.json'))
    if ($LASTEXITCODE -ne 0) {
        $failures += ('{0}: Matchete C5 pilot exit {1}' -f $model, $LASTEXITCODE)
        Write-Warning ('{0}: C5 blocked; continuing with independent RGBeta source' -f $model)
    }
    Write-Host "`n=== $model RGBeta source pilot ==="
    & $python (Join-Path $repo 'scripts\RunGeneralT3RGSourcePilot.py') --proof $proof --model $model --input $configFile --output (Join-Path $out ($model + '_rgbeta_source.json'))
    if ($LASTEXITCODE -ne 0) {
        $failures += ('{0}: RGBeta source pilot exit {1}' -f $model, $LASTEXITCODE)
        Write-Warning ('{0}: RGBeta source blocked; continuing with other models' -f $model)
    }
}
if ($failures.Count -gt 0) {
    Write-Host "`nStage 11 completed with errors; inspect output JSON files:" -ForegroundColor Yellow
    foreach ($f in $failures) { Write-Host ('  ' + $f) }
    exit 2
}
Write-Host "`nStage 11 experimental probes completed; no production files changed." -ForegroundColor Green
