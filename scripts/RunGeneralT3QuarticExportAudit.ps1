param(
    [switch] $NumericProbe
)
$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $PSScriptRoot
Push-Location $root
try {
    $argv = @('-m', 'scripts.RunGeneralT3QuarticExportAudit')
    if ($NumericProbe) { $argv += '--numeric-probe' }
    & python @argv
    if ($LASTEXITCODE -ne 0) { throw "Stage 21 failed with exit code $LASTEXITCODE" }
}
finally { Pop-Location }
