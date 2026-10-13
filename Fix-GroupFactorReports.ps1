param([string]$ProjectRoot = (Get-Location).Path)
$ErrorActionPreference = 'Stop'
$target = Join-Path $ProjectRoot 'Reports\GroupFactorReports.py'
$uri = 'https://raw.githubusercontent.com/Frankyang2002/T3neutrinomassmodelpipeline/main/Reports/GroupFactorReports.py'
$source = (Invoke-WebRequest -Uri $uri -UseBasicParsing).Content
if ($source -is [byte[]]) { $source = [Text.Encoding]::UTF8.GetString($source) }
if (-not $source.Contains('MAX_RUN_COLUMNS = 5')) { throw 'Unexpected GitHub source: MAX_RUN_COLUMNS marker missing.' }
if (-not $source.Contains('def _comparison_table(')) { throw 'Unexpected GitHub source: comparison table missing.' }
# Only modify the report layout; no beta functions or coefficients are changed.
$source = $source.Replace('MAX_RUN_COLUMNS = 5', 'MAX_RUN_COLUMNS = 3')
# Three equal bounded columns + one wider structure column.  \linewidth remains the upper bound.
$source = $source.Replace('"@{}p{0.31\\linewidth}" + "c" * len(chunk) + "@{}"', '"@{}p{0.37\\linewidth}" + "p{0.19\\linewidth}" * len(chunk) + "@{}"')
$source = $source.Replace('"@{}p{0.38\\linewidth}" + "c" * len(chunk) + "@{}"', '"@{}p{0.37\\linewidth}" + "p{0.19\\linewidth}" * len(chunk) + "@{}"')
$source = $source.Replace('"@{}l" + "c" * len(chunk) + "@{}"', '"@{}p{0.37\\linewidth}" + "p{0.19\\linewidth}" * len(chunk) + "@{}"')
# Fit overlong mathematical cells individually, retaining the full content.
# The helper is applied only inside comparison-table generators.
$anchor = 'def _comparison_table('
$helper = @'
def _fit_group_factor_math(latex: str) -> str:
    """Fit an inline mathematical table entry within its fixed-width column."""
    return rf"\adjustbox{{max width=\linewidth}}{{{latex}}}"


'@
$source = $source.Replace($anchor, $helper + $anchor)
# Limit the substitutions to the four comparison table functions.
$functions = @('_comparison_table', '_fermion_mass_comparison_table', '_saved_rge_term_table')
foreach ($name in $functions) {
    $pattern = '(?ms)(^def ' + [regex]::Escape($name) + '\(.*?)(?=^def |\z)'
    $match = [regex]::Match($source, $pattern)
    if (-not $match.Success) { throw "Could not locate function $name" }
    $block = $match.Value
    # Wrap every table cell after row construction, including mathematical headings.
    $block = $block.Replace('lines.append(" & ".join(header) + r" \\")', 'lines.append(" & ".join(_fit_group_factor_math(cell) for cell in header) + r" \\")')
    $block = $block.Replace('lines.append(" & ".join(row) + r" \\")', 'lines.append(" & ".join(_fit_group_factor_math(cell) for cell in row) + r" \\")')
    $block = $block.Replace('" & ".join(', '" & ".join(') # keep output syntax stable
    $source = $source.Substring(0, $match.Index) + $block + $source.Substring($match.Index + $match.Length)
}
# Ensure adjustbox is loaded in the report preamble.
$oldPackages = '"amsmath,amssymb,booktabs,longtable,array"'
if (-not $source.Contains($oldPackages)) { throw 'Expected LaTeX package list not found.' }
$source = $source.Replace($oldPackages, '"amsmath,amssymb,booktabs,longtable,array,adjustbox"')
$folder = Split-Path $target
New-Item -ItemType Directory -Path $folder -Force | Out-Null
if (Test-Path $target) { Copy-Item $target ($target + '.bak') -Force }
[IO.File]::WriteAllText($target, $source, [Text.UTF8Encoding]::new($false))
Write-Host "Generated complete replacement: $target"
Write-Host 'Source: GitHub main. No remote files modified.'
python -m py_compile $target
if ($LASTEXITCODE -ne 0) { throw 'Python syntax validation failed.' }
Write-Host 'Python syntax validation passed.'
