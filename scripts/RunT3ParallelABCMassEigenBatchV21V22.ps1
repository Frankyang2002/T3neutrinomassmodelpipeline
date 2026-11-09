# V21/V22 read-only bridges to EXISTING B/C mass-basis source reports.
param([switch]$VerifyOnly)
$ErrorActionPreference='Stop'
$root=Split-Path -Parent $PSScriptRoot
Set-Location -LiteralPath $root
$py=Join-Path $root '.venv\Scripts\python.exe'
if(-not (Test-Path -LiteralPath $py)){throw "Missing Python: $py"}
& $py -m unittest discover -s tests -p 'test_t3_parallel_a_bc_mass_bridge_v21_v22.py' -v
if($LASTEXITCODE -ne 0){throw "v21/v22 unit tests failed: $LASTEXITCODE"}
$base='output\parallel_t3\wilson'
if(-not (Test-Path -LiteralPath (Join-Path $base 'full_scalar_census_conditional_v20.json'))){
    throw 'Missing verified v20 scalar census. Run the v18-v20 batch first.'
}
foreach($model in @('B_m1','C_m1')){
    $path="output\general_t3_mass_basis\${model}_final_mass_basis.json"
    if(-not (Test-Path -LiteralPath $path)){
        throw "Missing existing experimental scalar mass basis: $path. Run the repository's RunGeneralT3MixedScalarMassBasis first."
    }
}
if($VerifyOnly){
    foreach($p in @('bc_native_mass_eigen_complex_conditional_v21.json','bc_quartic_mass_readiness_conditional_v22.json')){
        if(-not (Test-Path (Join-Path $base $p))){throw "VerifyOnly: missing $p"}
    }
}
& $py -m scripts.RunT3ParallelABCMassProvenanceV22 --project-root .
if($LASTEXITCODE -ne 0){throw "v21/v22 source provenance blocked: $LASTEXITCODE"}
# v19 uses historical Stage22 CG ordering; v18 exports by structure order.
# Preserve the input witnesses. Verify the named column permutation first.
& $py -m scripts.RunT3ParallelABCMassSubspaceOrderV23 --project-root .
if($LASTEXITCODE -ne 0){throw "v23 B/C seven-source-axis order blocked: $LASTEXITCODE"}
& $py -m scripts.RunT3ParallelABCMassEigenNativeV21 --project-root . --resume
if($LASTEXITCODE -ne 0){throw "Native B/C mass-eigen tensor transport v21 blocked: $LASTEXITCODE"}
& $py -m scripts.RunT3ParallelABCMassBridgeGateV22 --project-root . --resume
if($LASTEXITCODE -ne 0){throw "Source-conditional readiness v22 blocked: $LASTEXITCODE"}
Write-Host 'v21/v22 completed conditionally. No physical Wilson evolution or Match[].'
