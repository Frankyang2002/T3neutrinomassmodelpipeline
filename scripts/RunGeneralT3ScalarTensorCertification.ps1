# Stage 22 local-only workflow. Production files and GitHub are untouched.
param(
    [switch]$StaticOnly,
    [switch]$ExportMatchete,
    [switch]$ExportRGBetaDeclarations,
    [switch]$ExportRGBetaInvariants,
    [switch]$NativeDiagnosticsOnly,
    [string]$MatcheteTensors = "output/general_t3_stage22/matchete_tensors.json",
    [string]$RGBetaNativeTensors = "output/general_t3_stage22/rgbeta_native_tensors.json"
)
$ErrorActionPreference = "Stop"
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
Push-Location $project
try {
    python -m scripts.RunGeneralT3ScalarTensorCertification --static-only
    if ($LASTEXITCODE -ne 0) { throw "Stage22 static source audit FAILED" }
    if ($StaticOnly) { return }

    if ($ExportMatchete) {
        wolframscript -file .\scripts\ExportGeneralT3ScalarTensorMatchete.wl $MatcheteTensors
        if ($LASTEXITCODE -ne 0) {
            throw "Stage22 Matchete native tensor export FAILED; no certificate can be issued"
        }
    }
    if ($ExportRGBetaDeclarations) {
        wolframscript -file .\scripts\ExportGeneralT3ScalarTensorRGBetaDeclarations.wl `
            "output/general_t3_stage22/rgbeta_native_declarations.json"
        if ($LASTEXITCODE -ne 0) { throw "Stage22 RGBeta declaration export FAILED" }
    }
    # Treat -ExportRGBetaDeclarations as the legacy 'export all RGBeta evidence'
    # option, adding the independently evaluated native invariant expressions.
    if ($ExportRGBetaInvariants -or $ExportRGBetaDeclarations) {
        wolframscript -file .\scripts\ExportGeneralT3ScalarTensorRGBetaInvariants.wl `
            "output/general_t3_stage22/rgbeta_native_invariants.json"
        if ($LASTEXITCODE -ne 0) { throw "Stage22 RGBeta native invariant expression export FAILED" }
    }

    if ($NativeDiagnosticsOnly) {
        python -m scripts.RunGeneralT3ScalarTensorCertification --inspect-native `
            --matchete-tensors $MatcheteTensors
        if ($LASTEXITCODE -ne 0) { throw "Stage22 native evidence integrity FAILED" }
        Write-Host "Stage22 native diagnostic exports passed integrity checks; physical normalisation remains UNCERTIFIED."
        return
    }

    python -m scripts.RunGeneralT3ScalarTensorCertification --matchete-tensors $MatcheteTensors `
        --rgbeta-tensors $RGBetaNativeTensors
    if ($LASTEXITCODE -ne 0) {
        throw "Stage22 physically UNCERTIFIED: inspect output/general_t3_stage22/stage22_summary.json"
    }
}
finally { Pop-Location }
