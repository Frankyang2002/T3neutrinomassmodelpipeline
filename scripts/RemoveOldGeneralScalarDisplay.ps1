# Removes ONLY the superseded general-scalar display entry point and its old tests.
# Run from any PowerShell working directory. No physics data or original Lambda
# report generators are modified. Never invokes GitHub or git.
$ErrorActionPreference = 'Stop'
$repo = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$obsolete = @(
    'scripts\BuildGeneralScalarInteractiveDisplay.py',
    'tests\test_general_scalar_interactive_display.py'
)
foreach ($relative in $obsolete) {
    $path = Join-Path $repo $relative
    if (Test-Path -LiteralPath $path -PathType Leaf) {
        Remove-Item -LiteralPath $path -Force
        Write-Host "Removed $relative"
    } else {
        Write-Host "Already absent: $relative"
    }
}
Write-Host 'Original Display Lambda generators and numerical data left unchanged.'
