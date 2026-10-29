<#
  Stage 12 extra-only neutrino-mass and phase scan.
  Reuses the successfully generated Stage 11 Matchete hard C5 JSONs.
  Independent of Matchete/RGBeta coupling normalisation and full RGE.
#>
param(
    [string]$InputRoot = 'output\general_t3_numerical_pilot',
    [string]$OutputRoot = 'output\general_t3_neutrino_impact',
    [double]$VevGeV = 246.22
)
$ErrorActionPreference = 'Stop'
$repo = Split-Path -Parent $PSScriptRoot
$python = (Get-Command python -ErrorAction Stop).Source
function Resolve-RepoPath([string]$p) {
    if ([System.IO.Path]::IsPathRooted($p)) { return $p }
    return Join-Path $repo $p
}
& $python (Join-Path $repo 'scripts\RunGeneralT3ExtraC5PhaseResponse.py') --input-root (Resolve-RepoPath $InputRoot) --output-root (Resolve-RepoPath $OutputRoot) --vev-gev $VevGeV
if ($LASTEXITCODE -ne 0) {
    throw ('Stage 12 mass-impact analysis failed with exit code {0}' -f $LASTEXITCODE)
}
Write-Host 'Stage 12 completed: extra-only masses and toy coupling-phase scan.' -ForegroundColor Green
Write-Host 'The result is not a full physical neutrino spectrum or an RG-improved result.'
