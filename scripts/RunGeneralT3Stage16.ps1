param(
    [string]$Config = 'docs\general_t3_threshold_benchmark.json',
    [switch]$RenormalisableOnly
)
$ErrorActionPreference = 'Stop'
$repo = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
Push-Location $repo
try {
    $command = @('-B', (Join-Path $repo 'scripts\RunGeneralT3Stage16.py'), '--config', $Config)
    if ($RenormalisableOnly) { $command += '--renormalisable-only' }
    & python @command
    if ($LASTEXITCODE -ne 0) { throw ('Stage 16 failed with exit code {0}; inspect output\general_t3_stage16\stage16_summary.json' -f $LASTEXITCODE) }
} finally {
    Pop-Location
}
