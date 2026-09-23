$ErrorActionPreference = "Stop"

$RepoRoot = (Get-Location).Path

$ReportGeneration = Join-Path $RepoRoot "Reports\ReportGeneration.py"
$GroupFactors = Join-Path $RepoRoot "Reports\GroupFactorReports.py"

$ReportGenerationOut = Join-Path $RepoRoot "Reports\ReportGeneration.py.generated"
$GroupFactorsOut = Join-Path $RepoRoot "Reports\GroupFactorReports.py.generated"

if (-not (Test-Path $ReportGeneration)) {
    throw "Missing file: $ReportGeneration"
}

if (-not (Test-Path $GroupFactors)) {
    throw "Missing file: $GroupFactors"
}


function Replace-Required {
    param(
        [string]$Text,
        [string]$Old,
        [string]$New,
        [string]$Description
    )

    if (-not $Text.Contains($Old)) {
        throw "Could not find expected text for: $Description"
    }

    Write-Host "[OK] $Description"
    return $Text.Replace($Old, $New)
}


function Replace-Optional {
    param(
        [string]$Text,
        [string]$Old,
        [string]$New,
        [string]$Description
    )

    if ($Text.Contains($Old)) {
        Write-Host "[OK] $Description"
        return $Text.Replace($Old, $New)
    }

    Write-Host "[INFO] Already changed or not present: $Description"
    return $Text
}


# ============================================================================
# Reports\ReportGeneration.py
# ============================================================================

Write-Host ""
Write-Host "Updating ReportGeneration.py ..."

$rg = Get-Content $ReportGeneration -Raw


# ----------------------------------------------------------------------------
# 1. Human-facing CG notation used by Lagrangian/report normalization
# ----------------------------------------------------------------------------

$oldCgMap = @'
    cg_map = {
        "T3Y1CG": r"\mathcal{I}_{y_1}",
        "T3Y2CG": r"\mathcal{I}_{y_2}",
        "T3MixCG": r"\mathcal{I}_{T3}",
    }
'@

$newCgMap = @'
    cg_map = {
        "T3Y1CG": r"\mathcal{C}_{LS_1F}",
        "T3Y2CG": r"\mathcal{C}_{LS_2F}",
        "T3MixCG": r"\mathcal{I}_{LLS_1S_2}",
    }
'@

$rg = Replace-Required `
    -Text $rg `
    -Old $oldCgMap `
    -New $newCgMap `
    -Description "CG notation: T3Y1CG/T3Y2CG/T3MixCG"


# ----------------------------------------------------------------------------
# 2. Quartic interaction notation
# ----------------------------------------------------------------------------

$oldCross = '        r"$\lambda_{12}^{(\times)}$ & crossed independent $S_1$--$S_2$ contraction \\",'
$newCross = '        r"$\lambda_{12}^{(\times)}$ & $\lambda_{12}^{(\times)}(S_1^\dagger S_2)(S_2^\dagger S_1)$ \\",'

$oldS1 = '        r"$\lambda_{S_1}^{(A)}$ & additional independent $S_1$ self-contraction \\",'
$newS1 = '        r"$\lambda_{S_1}^{(A)}$ & $\lambda_{S_1}^{(A)}(S_1^\dagger T^A S_1)(S_1^\dagger T^A S_1)$ \\",'

$oldS2 = '        r"$\lambda_{S_2}^{(A)}$ & additional independent $S_2$ self-contraction \\",'
$newS2 = '        r"$\lambda_{S_2}^{(A)}$ & $\lambda_{S_2}^{(A)}(S_2^\dagger T^A S_2)(S_2^\dagger T^A S_2)$ \\",'

$rg = Replace-Required `
    -Text $rg `
    -Old $oldCross `
    -New $newCross `
    -Description "lambda12^(x) explicit crossed interaction"

$rg = Replace-Required `
    -Text $rg `
    -Old $oldS1 `
    -New $newS1 `
    -Description "lambdaS1^(A) explicit self-interaction"

$rg = Replace-Required `
    -Text $rg `
    -Old $oldS2 `
    -New $newS2 `
    -Description "lambdaS2^(A) explicit self-interaction"


Set-Content `
    -Path $ReportGenerationOut `
    -Value $rg `
    -Encoding UTF8

Write-Host ""
Write-Host "Generated full replacement:"
Write-Host "  $ReportGenerationOut"


# ============================================================================
# Reports\GroupFactorReports.py
# ============================================================================

Write-Host ""
Write-Host "Updating GroupFactorReports.py ..."

$gf = Get-Content $GroupFactors -Raw


# ----------------------------------------------------------------------------
# 3. Matched-operator SU(2) invariant rendering
# ----------------------------------------------------------------------------

$oldMatchingCg = @'
    rendered: list[str] = []

    for name in names:
        if name == "eps[SU2L]":
            rendered.append(r"\epsilon_{SU(2)_L}")
        else:
            rendered.append(
                rf"\mathrm{{{latex_escape_text(name)}}}"
            )

    return r"\,".join(rendered)
'@

$newMatchingCg = @'
    rendered: list[str] = []

    invariant_tex = {
        "eps[SU2L]": r"\epsilon",
        "T3Y1CG": r"\mathcal{C}_{LS_1F}",
        "T3Y2CG": r"\mathcal{C}_{LS_2F}",
        "T3MixCG": r"\mathcal{I}_{LLS_1S_2}",
    }

    for name in names:
        rendered.append(
            invariant_tex.get(
                name,
                rf"\mathrm{{{latex_escape_text(name)}}}",
            )
        )

    return r"\,".join(rendered)
'@

$gf = Replace-Required `
    -Text $gf `
    -Old $oldMatchingCg `
    -New $newMatchingCg `
    -Description "matched-operator CG names -> physics notation"


# ----------------------------------------------------------------------------
# 4. GroupFactorReports' own notation tables
# ----------------------------------------------------------------------------

$oldGfCross = '        (r"$\lambda_{12}^{(\times)}$", r"crossed independent $S_1$--$S_2$ contraction"),'
$newGfCross = '        (r"$\lambda_{12}^{(\times)}$", r"$\lambda_{12}^{(\times)}(S_1^\dagger S_2)(S_2^\dagger S_1)$"),'

$gf = Replace-Optional `
    -Text $gf `
    -Old $oldGfCross `
    -New $newGfCross `
    -Description "GroupFactorReports lambda12^(x) notation"


# If these rows already exist in a local notation table, replace the prose.
$oldGfS1 = '        (r"$\lambda_{S_1}^{(A)}$", r"additional independent $S_1$ self-contraction"),'
$newGfS1 = '        (r"$\lambda_{S_1}^{(A)}$", r"$\lambda_{S_1}^{(A)}(S_1^\dagger T^AS_1)(S_1^\dagger T^AS_1)$"),'

$oldGfS2 = '        (r"$\lambda_{S_2}^{(A)}$", r"additional independent $S_2$ self-contraction"),'
$newGfS2 = '        (r"$\lambda_{S_2}^{(A)}$", r"$\lambda_{S_2}^{(A)}(S_2^\dagger T^AS_2)(S_2^\dagger T^AS_2)$"),'

$gf = Replace-Optional `
    -Text $gf `
    -Old $oldGfS1 `
    -New $newGfS1 `
    -Description "GroupFactorReports lambdaS1^(A) notation"

$gf = Replace-Optional `
    -Text $gf `
    -Old $oldGfS2 `
    -New $newGfS2 `
    -Description "GroupFactorReports lambdaS2^(A) notation"


Set-Content `
    -Path $GroupFactorsOut `
    -Value $gf `
    -Encoding UTF8

Write-Host ""
Write-Host "Generated full replacement:"
Write-Host "  $GroupFactorsOut"


# ============================================================================
# Final checks
# ============================================================================

Write-Host ""
Write-Host "============================================================"
Write-Host "FULL REPLACEMENT FILES GENERATED"
Write-Host "============================================================"
Write-Host ""
Write-Host "Generated:"
Write-Host "  Reports\ReportGeneration.py.generated"
Write-Host "  Reports\GroupFactorReports.py.generated"
Write-Host ""
Write-Host "Intended repository paths:"
Write-Host "  Reports\ReportGeneration.py"
Write-Host "  Reports\GroupFactorReports.py"
Write-Host ""
Write-Host "The live files have NOT been overwritten."
Write-Host ""