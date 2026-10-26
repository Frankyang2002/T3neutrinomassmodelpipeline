# Run the experimental T3 scalar-quartic validation across one or many models.
# No GitHub access and no production output files are overwritten.
#
# Examples from the repository root:
# powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralScalarTrial.ps1 -Scope BC -Phase Audit
# powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralScalarTrial.ps1 -Scope All -Phase Audit
# powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralScalarTrial.ps1 -Scope B -Phase Matching
#
# Windows PowerShell 5.1 compatibility: expand names immediately before : as ${name}.
# WARNING: Matching phase is symbolic-only. The RGBeta/numerical state does NOT
# include every additional coupling. Do not treat matching output as a verified
# end-to-end general-T3 prediction.

param(
    [ValidateSet('B','BC','Problems','All')]
    [string]$Scope = 'BC',
    [ValidateSet('Audit','Matching','Both')]
    [string]$Phase = 'Audit',
    [string]$OutputRoot = 'output\general_t3_trial',
    [int]$EFTOrder = 5,
    [int]$LoopOrder = 1
)

$ErrorActionPreference = 'Stop'
$repo = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$wolframscript = (Get-Command 'wolframscript' -ErrorAction Stop).Source
$dest = if ([System.IO.Path]::IsPathRooted($OutputRoot)) {
    $OutputRoot
} else { Join-Path $repo $OutputRoot }
New-Item -ItemType Directory -Path $dest -Force | Out-Null

$points = @(
    @{Class='A'; Alpha='m4'}, @{Class='A'; Alpha='m2'}, @{Class='A'; Alpha='p0'},
    @{Class='B'; Alpha='m3'}, @{Class='B'; Alpha='m1'}, @{Class='B'; Alpha='p1'},
    @{Class='C'; Alpha='m3'}, @{Class='C'; Alpha='m1'}, @{Class='C'; Alpha='p1'},
    @{Class='D'; Alpha='m2'}, @{Class='D'; Alpha='p0'}, @{Class='D'; Alpha='p2'},
    @{Class='E'; Alpha='m4'}, @{Class='E'; Alpha='m2'}, @{Class='E'; Alpha='p0'}, @{Class='E'; Alpha='p2'}
)
$points = @($points | Where-Object {
    if ($Scope -eq 'All') { return $true }
    if ($Scope -eq 'BC') { return $_.Class -in @('B','C') -and $_.Alpha -eq 'm1' }
    if ($Scope -eq 'Problems') { return $_.Class -in @('B','C') -and $_.Alpha -in @('m3','p1') }
    return $_.Class -eq 'B' -and $_.Alpha -eq 'm1'
})

# Expected quartic counts are a gate, not merely an informational report.
# This prevents successful CheckLagrangian calls for a silently incomplete set.
$expectedQuartics = @{
    A_m4=0; A_m2=5; A_p0=5
    B_m3=1; B_m1=7; B_p1=1
    C_m3=1; C_m1=7; C_p1=1
    D_m2=5; D_p0=5; D_p2=0
    E_m4=0; E_m2=6; E_p0=6; E_p2=0
}

$results = [System.Collections.Generic.List[object]]::new()
$originalFlag = [Environment]::GetEnvironmentVariable('T3_EXPERIMENTAL_GENERAL_QUARTICS', 'Process')

function Invoke-WolframTrial {
    param([string[]]$Arguments, [string]$LogFile)
    $ErrorActionPreference = 'Continue'
    & $wolframscript @Arguments 2>&1 | ForEach-Object {
        $_ | Out-File -FilePath $LogFile -Append -Encoding utf8
        Write-Host $_
    }
    $status = $LASTEXITCODE
    if ($null -eq $status) { $status = 999 }
    return [int]$status
}

try {
    foreach ($point in $points) {
        $name = $point.Class + '_' + $point.Alpha
        $modelDir = Join-Path $dest $name
        New-Item -ItemType Directory -Path $modelDir -Force | Out-Null

        $auditSuccess = $true
        if ($Phase -in @('Audit','Both')) {
            $auditJson = Join-Path $modelDir 'audit.json'
            if (Test-Path $auditJson) { Remove-Item -LiteralPath $auditJson }
            $log = Join-Path $modelDir 'audit.log'
            if (Test-Path $log) { Remove-Item $log }
            Write-Host "`n=== Audit $name ===" -ForegroundColor Cyan
            $rc = Invoke-WolframTrial -Arguments @('-file', (Join-Path $repo 'Lagrangian\RunGeneralScalarAudit.wl'), $auditJson, $point.Class, $point.Alpha) -LogFile $log
            $status = 'NoJSON'
            $candidateCount = $null
            $accepted = $null
            if (Test-Path $auditJson) {
                try {
                    $body = Get-Content -Raw -LiteralPath $auditJson | ConvertFrom-Json
                    $status = $body.Status
                    $candidateCount = $body.ExtraQuarticCandidateCount
                    $accepted = $body.AcceptedExtraQuarticCount
                } catch { $status = 'UnparseableJSON' }
            }
            $expected = $expectedQuartics[$name]
            $auditSuccess = ($rc -eq 0 -and $status -eq 'Success' -and
                $null -ne $candidateCount -and $null -ne $accepted -and
                [int]$candidateCount -eq $expected -and [int]$accepted -eq $expected)
            if ($rc -eq 0 -and $status -eq 'Success' -and -not $auditSuccess) {
                $status = 'IncompleteOperatorBasis'
                Write-Warning ("{0}: expected {1} quartics, observed {2}, accepted {3}" -f $name, $expected, $candidateCount, $accepted)
            }
            $results.Add([PSCustomObject]@{Model=$name; Phase='Audit'; ExitCode=$rc; Status=$status; Candidates=$candidateCount; Accepted=$accepted; File=$auditJson})
        }

        if ($Phase -eq 'Both' -and -not $auditSuccess) {
            Write-Warning "Skipping matching for ${name}: symbolic validation failed."
        }
        if ($Phase -eq 'Matching' -or ($Phase -eq 'Both' -and $auditSuccess)) {
            foreach ($variant in @('baseline','extended')) {
                $runDir = Join-Path $modelDir $variant
                New-Item -ItemType Directory -Path $runDir -Force | Out-Null
                # Prevent a failed rerun from leaving a stale Weinberg artifact
                # that the symbolic comparison could mistake for new output.
                foreach ($artifact in @('c5_coefficient.txt','weinberg_holomorphic_terms.txt')) {
                    $oldArtifact = Join-Path $runDir $artifact
                    if (Test-Path $oldArtifact) { Remove-Item -LiteralPath $oldArtifact }
                }
                $log = Join-Path $modelDir ($variant + '_matching.log')
                if (Test-Path $log) { Remove-Item $log }
                if ($variant -eq 'extended') {
                    [Environment]::SetEnvironmentVariable('T3_EXPERIMENTAL_GENERAL_QUARTICS', '1', 'Process')
                } else {
                    [Environment]::SetEnvironmentVariable('T3_EXPERIMENTAL_GENERAL_QUARTICS', $null, 'Process')
                }
                Write-Host "`n=== Matching $name ($variant) ===" -ForegroundColor Cyan
                $rc = Invoke-WolframTrial -Arguments @('-file', (Join-Path $repo 'Lagrangian\RunModel.wl'), $runDir, [string]$EFTOrder, [string]$LoopOrder, $point.Class, $point.Alpha, 'DEBUG') -LogFile $log
                $results.Add([PSCustomObject]@{Model=$name; Phase="Matching-$variant"; ExitCode=$rc; Status=$(if ($rc -eq 0) {'Completed'} else {'Failed'}); Candidates=$null; Accepted=$null; File=$runDir})
            }
        }
    }
}
finally {
    [Environment]::SetEnvironmentVariable('T3_EXPERIMENTAL_GENERAL_QUARTICS', $originalFlag, 'Process')
}

$summary = Join-Path $dest 'trial_runs.json'
@($results.ToArray()) | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath $summary -Encoding utf8
$results | Format-Table Model, Phase, ExitCode, Status, Candidates, Accepted -AutoSize
Write-Host "Run summary: $summary"
if (@($results | Where-Object { $_.ExitCode -ne 0 -or $_.Status -in @('IncompleteOperatorBasis', 'Failed', 'NoJSON', 'UnparseableJSON') }).Count -gt 0) { exit 1 }
exit 0
