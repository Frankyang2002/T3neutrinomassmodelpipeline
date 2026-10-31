# Run from the local T3 repository root. No GitHub writes or Wolfram sessions.
param(
    [string]$NumericManifest = '',
    [string]$Output = 'output\audits\t3_bc_sequential_batch_v14\summary.json'
)
$ErrorActionPreference = 'Stop'
$root = (Get-Location).Path
$required = @(
    'scripts\BuildT3BCYFiniteBubbleInterface.experimental.py',
    'scripts\BuildT3BCFHardBoundaryInterface.experimental.py',
    'scripts\BuildT3BCSequentialMatchingBatch.experimental.py',
    'output\audits\t3_bc_y_warsaw_direct_identity_v11\summary_live.json',
    'output\audits\t3_bc_scalar_finite_term_extraction_v1\summary.json',
    'output\audits\t3_bc_physical_y_matching_gate_v1\summary.json'
)
foreach ($p in $required) {
    if (-not (Test-Path -LiteralPath $p -PathType Leaf)) {
        throw "Required local source/evidence is missing: $p"
    }
}
function Invoke-PythonStage([string[]]$ScriptArgs) {
    & python @ScriptArgs
    if ($LASTEXITCODE -ne 0) { throw "Batch Python stage failed ($LASTEXITCODE): $($ScriptArgs -join ' ')" }
}
Write-Host 'T3_C_BATCH: self-testing F, Y and combined F/Y/X stages'
Invoke-PythonStage -ScriptArgs @('scripts\BuildT3BCYFiniteBubbleInterface.experimental.py','--self-test')
Invoke-PythonStage -ScriptArgs @('scripts\BuildT3BCFHardBoundaryInterface.experimental.py','--self-test')
Invoke-PythonStage -ScriptArgs @('scripts\BuildT3BCSequentialMatchingBatch.experimental.py','--self-test')
Write-Host 'T3_C_BATCH: regenerating Y restricted and F formal ledgers'
Invoke-PythonStage -ScriptArgs @('scripts\BuildT3BCYFiniteBubbleInterface.experimental.py',
    '--output','output\audits\t3_bc_y_conditional_finite_interface_v1\summary.json')
Invoke-PythonStage -ScriptArgs @('scripts\BuildT3BCFHardBoundaryInterface.experimental.py',
    '--output','output\audits\t3_bc_f_hard_boundary_v1\summary.json')
Write-Host 'T3_C_BATCH: auditing Y LLXX and X LLHH; assembling fixed-order ledger'
$argsList = @('scripts\BuildT3BCSequentialMatchingBatch.experimental.py','--output',$Output)
if ($NumericManifest -ne '') {
    if (-not (Test-Path -LiteralPath $NumericManifest -PathType Leaf)) {
        throw "Numeric manifest not found: $NumericManifest"
    }
    $argsList += @('--numeric-manifest',$NumericManifest)
}
Invoke-PythonStage -ScriptArgs $argsList
Write-Host "T3_C_BATCH: done; $Output"
Write-Host 'T3_C_BATCH: physical C5 remains NOT CERTIFIED'
