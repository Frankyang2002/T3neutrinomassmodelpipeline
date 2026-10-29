# Stage31: isolated RGBeta export and numerical UV running for 15 models.
# Run in repository root: .\scripts\RunGeneralT3CoupledUVRunning.ps1
param(
    [string[]]$Models = @('A_-2','A_-4','B_+1','B_-1','B_-3','C_+1',
                          'C_-1','C_-3','D_+0','D_+2','D_-2',
                          'E_+0','E_+2','E_-2','E_-4')
)
$ErrorActionPreference = 'Stop'
$root = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
Push-Location $root
try {
    foreach ($model in $Models) {
        if ($model -notin @('A_-2','A_-4','B_+1','B_-1','B_-3',
                           'C_+1','C_-1','C_-3','D_+0','D_+2','D_-2',
                           'E_+0','E_+2','E_-2','E_-4')) {
            throw "Invalid Stage31 model: $model"
        }
        $token = $model.Replace('+','p').Replace('-','m')
        $path = "output/general_t3_stage22/stage31_${token}_full_uv_betas.json"
        Write-Host "Stage31 exporting $model ..."
        & wolframscript -file .\scripts\ExportGeneralT3AllModelsUVOneLoop.wl $model $path
        if ($LASTEXITCODE -ne 0) { throw "Stage31 RGBeta export failed: $model" }
        & python -m scripts.CoupledUVRunning --model $model
        if ($LASTEXITCODE -ne 0) { throw "Stage31 UV integration failed: $model" }
    }
    if ($Models.Count -eq 15) {
        & python -m scripts.CoupledUVRunning --summary
        if ($LASTEXITCODE -ne 0) { throw 'Stage31 16-model summary failed' }
    } else {
        Write-Host 'Partial model selection: summary skipped; rerun full default list to summarize.'
    }
} finally { Pop-Location }
