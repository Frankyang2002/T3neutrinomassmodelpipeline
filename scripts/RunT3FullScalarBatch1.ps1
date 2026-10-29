<#
Bulk non-destructive continuation of the *existing* T3 general-scalar pipeline.
Project path: scripts/RunT3FullScalarBatch1.ps1
Requires local existing project sources, Python and (for -ExportBetas) wolframscript.
This does not claim complete physical matching or modify the production Lagrangian.
#>
[CmdletBinding()]
param(
    [switch]$ExportBetas,
    [switch]$RunTests,
    [switch]$SkipQuadraticPrototype
)
$ErrorActionPreference = 'Stop'
$repo = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
Push-Location $repo
try {
    $python = (Get-Command python -ErrorAction Stop).Source
    $stages = [System.Collections.Generic.List[object]]::new()
    $resultsDir = Join-Path $repo 'output\audits\full_scalar_bulk_batch1'
    New-Item -ItemType Directory -Path $resultsDir -Force | Out-Null
    $stamp = (Get-Date).ToUniversalTime().ToString('o')
    function Invoke-Stage {
        param([string]$Name, [string[]]$CommandArgs)
        Write-Host "`n=== $Name ===" -ForegroundColor Cyan
        $log = Join-Path $resultsDir (($Name -replace '[^a-zA-Z0-9_-]', '_') + '.log')
        & $python @CommandArgs 2>&1 | Tee-Object -FilePath $log
        $code = $LASTEXITCODE
        $script:stages.Add([pscustomobject]@{stage=$Name; exit_code=$code; log=$log})
        if ($code -ne 0) { throw "$Name failed (exit $code); see $log" }
    }
    Invoke-Stage 'independent_quadratic_basis' @('-m', 'scripts.AuditGeneralT3QuadraticBasis')
    if ($ExportBetas) {
        if (-not (Get-Command wolframscript -ErrorAction SilentlyContinue)) {
            throw 'ExportBetas requested but wolframscript is not on PATH'
        }
        Invoke-Stage 'quartic_uv_and_eft1_betas_all_16' @('-m', 'scripts.RunGeneralT3GeneralQuarticRGBeta', '--check-gauge', '--export-betas', '--stage', 'both')
    } else {
        Invoke-Stage 'quartic_basis_all_16' @('-m', 'scripts.RunGeneralT3GeneralQuarticRGBeta', '--stage', 'both')
    }
    if (-not $SkipQuadraticPrototype) {
        Invoke-Stage 'quadratic_mass_hessian_all_16' @('-m', 'scripts.RunGeneralT3QuadraticAudit')
    }
    if ($RunTests) {
        Invoke-Stage 'python_tests' @('-X','utf8','-m','pytest','-q')
    }
    $registryPath='RGE/general/generated/general_rgbeta_quartic_registry.json'
    if (-not (Test-Path $registryPath)) { throw "Missing registry $registryPath" }
    $registry = Get-Content -LiteralPath $registryPath -Raw -Encoding UTF8 | ConvertFrom-Json
    $modelCount = @($registry.models.PSObject.Properties).Count
    $quarticCount = 0; $quadraticCount = 0
    foreach ($prop in $registry.models.PSObject.Properties) {
        $quarticCount += [int]$prop.Value.total_additional_quartics
        $quadraticCount += [int]$prop.Value.quadratics_not_implemented
    }
    if ($modelCount -ne 16 -or $quarticCount -ne 50 -or $quadraticCount -ne 8) {
        throw "Unexpected scalar basis inventory $modelCount/$quarticCount/$quadraticCount"
    }
    $report = [ordered]@{
        status = 'BulkStructuralAndCandidateRGEExportOnly'
        generated_utc = $stamp
        models = $modelCount
        extra_quartics_across_models = $quarticCount
        extra_quadratics_across_models = $quadraticCount
        candidate_uv_eft1_betas_requested = [bool]$ExportBetas
        candidate_uv_eft1_betas_exported = [bool]$ExportBetas
        coupling_normalisations_fully_certified = $false
        general_quadratic_rg_and_threshold_matching_implemented = $false
        general_sequential_c5_matching_certified = $false
        physical_all_16_neutrino_predictions_certified = $false
        stages = @($stages.ToArray())
    }
    $dest=Join-Path $resultsDir 'batch1_status.json'
    ($report | ConvertTo-Json -Depth 8) | Set-Content -LiteralPath $dest -Encoding UTF8
    Write-Host "`nInventory: $modelCount models / $quarticCount extra quartics / $quadraticCount extra quadratics"
    Write-Host "Saved status: $dest"
    Write-Host 'Physical full-potential sequential matching remains UNVERIFIED.' -ForegroundColor Yellow
} finally {
    Pop-Location
}
