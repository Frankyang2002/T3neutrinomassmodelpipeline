<#
Stage 34: export native projected scalar-only EFT1 betas and integrate Stage 33
renormalisable boundaries. No modification of production runner or GitHub.
Run through: powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3Stage34.ps1
Optional single model: ... -Models 'A_+0'
#>
[CmdletBinding()]
param(
    [string[]]$Models = @(
        'A_+0','A_-2','A_-4','B_+1','B_-1','B_-3',
        'C_+1','C_-1','C_-3','D_+0','D_+2','D_-2',
        'E_+0','E_+2','E_-2','E_-4'
    )
)
$ErrorActionPreference = 'Stop'
$repo = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
Set-Location -LiteralPath $repo
$outdir = Join-Path $repo 'output\general_t3_stage22'
New-Item -ItemType Directory -Path $outdir -Force | Out-Null

# Expected representation and hypercharge parameter for the 16 ordinary models.
$spec = @{
    'A_+0' = @(1,3,2,0);  'A_-2' = @(1,3,2,-2); 'A_-4' = @(1,3,2,-4)
    'B_+1' = @(2,2,1,1);  'B_-1' = @(2,2,1,-1); 'B_-3' = @(2,2,1,-3)
    'C_+1' = @(2,2,3,1);  'C_-1' = @(2,2,3,-1); 'C_-3' = @(2,2,3,-3)
    'D_+0' = @(3,1,2,0);  'D_+2' = @(3,1,2,2);  'D_-2' = @(3,1,2,-2)
    'E_+0' = @(3,3,2,0);  'E_+2' = @(3,3,2,2);  'E_-2' = @(3,3,2,-2)
    'E_-4' = @(3,3,2,-4)
}
foreach ($model in $Models) {
    if (-not $spec.ContainsKey($model)) {
        throw "Stage34 invalid model: $model"
    }
    $rep = $spec[$model]
    $alpha = if ($rep[3] -lt 0) { 'm' + [string](-$rep[3]) } else { [string]$rep[3] }
    $safe = $model.Replace('+','p').Replace('-','m')
    $beta = Join-Path $outdir "stage34_${safe}_eft1_full_betas.json"
    Write-Host "Stage34 exporting $model EFT1 betas ..."
    & wolframscript -file .\scripts\ExportGeneralT3Stage34EFT1OneLoop.wl ([string]$rep[0]) ([string]$rep[1]) ([string]$rep[2]) $alpha $beta
    if ($LASTEXITCODE -ne 0) { throw "Stage34 RGBeta EFT1 export failed: $model" }
    & python -m scripts.Stage34IntermediateEFTRunning --model $model
    if ($LASTEXITCODE -ne 0) { throw "Stage34 supported EFT1 numerical integration failed: $model" }
}
if (-not $PSBoundParameters.ContainsKey('Models')) {
    & python -m scripts.Stage34IntermediateEFTRunning --summary
    if ($LASTEXITCODE -ne 0) { throw 'Stage34 combined 16-model summary failed' }
}
