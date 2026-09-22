$ErrorActionPreference = "Stop"

# ============================================================================
# Paths
# ============================================================================

$RepoRoot = (Get-Location).Path

$Target = Join-Path $RepoRoot "Reports\GroupFactorReports.py"
$Backup = Join-Path $RepoRoot "Reports\GroupFactorReports.py.before_complex_matching"
$Generated = Join-Path $RepoRoot "Reports\GroupFactorReports.py.generated"

if (-not (Test-Path $Target)) {
    throw "Cannot find $Target. Run this script from the repository root."
}


# ============================================================================
# Helpers
# ============================================================================

function Replace-Block {
    param(
        [Parameter(Mandatory = $true)]
        [string]$Text,

        [Parameter(Mandatory = $true)]
        [string]$StartMarker,

        [Parameter(Mandatory = $true)]
        [string]$EndMarker,

        [Parameter(Mandatory = $true)]
        [string]$Replacement
    )

    $start = $Text.IndexOf($StartMarker)

    if ($start -lt 0) {
        throw "Could not find start marker: $StartMarker"
    }

    $end = $Text.IndexOf($EndMarker, $start)

    if ($end -lt 0) {
        throw "Could not find end marker: $EndMarker"
    }

    if ($end -le $start) {
        throw "End marker occurs before start marker: $EndMarker"
    }

    $before = $Text.Substring(0, $start)
    $after = $Text.Substring($end)

    return ($before, $Replacement, $after) -join ""
}


# ============================================================================
# Load source
# ============================================================================

Write-Host ""
Write-Host "Source:"
Write-Host "  $Target"

Copy-Item $Target $Backup -Force

$text = Get-Content $Target -Raw

Write-Host ""
Write-Host "Backup:"
Write-Host "  $Backup"


# ============================================================================
# 1. Optional module documentation update
# ============================================================================

$oldDoc = @'
* the F-threshold Wilson tensor is reconstructed from the exported seed rather
  than hard-coded from A--E benchmark values.
'@

$newDoc = @'
* the F-threshold matched dimension-five operators are read directly from the
  exported complex Matchete seed; the real C_ijab expansion remains an internal
  representation used by the generic tensor-RGE machinery.
'@

if ($text.Contains($oldDoc)) {
    $text = $text.Replace($oldDoc, $newDoc)
    Write-Host "[OK] Updated module documentation."
}
else {
    Write-Host "[INFO] Module documentation differs locally; skipped wording update."
}


# ============================================================================
# 2. Replace report-only tensor builder import
# ============================================================================

$oldImport = "from RGE.running.eft1.EFT1TensorAdapters import build_eft1_wilson_tensor"

$newImport = @'
from RGE.running.eft1.EFT1TensorAdapters import (
    _term_prefactor as _eft1_wilson_term_prefactor,
)
'@

if ($text.Contains($oldImport)) {
    $text = $text.Replace(
        $oldImport,
        $newImport.TrimEnd()
    )

    Write-Host "[OK] Replaced real-tensor report import."
}
elseif ($text.Contains("_eft1_wilson_term_prefactor")) {
    Write-Host "[INFO] Complex matching prefactor import is already present."
}
else {
    throw @"
Could not find the expected EFT1TensorAdapters import.

Expected:
$oldImport
"@
}


# ============================================================================
# 3. Replace _matching_tensor_rows
# ============================================================================

if ($text.Contains("def _matching_tensor_rows(")) {

    $newMatchingFunctions = @'
def _matching_scalar_product_tex(term: Mapping[str, object]) -> str:
    """Render the two scalar legs of a matched term in the complex basis.

    The Matchete seed records both the total number of S1/S2 legs and the
    number of barred legs.  A barred complex scalar is displayed as a daggered
    field.  This presentation deliberately occurs before the conversion to
    global real scalar coordinates.
    """
    pieces: list[str] = []

    for index in (1, 2):
        count = int(term.get(f"Scalar{index}Count", 0) or 0)
        barred = int(term.get(f"BarredScalar{index}Count", 0) or 0)

        if barred < 0 or barred > count:
            raise ValueError(
                f"Invalid barred-scalar count for S{index}: "
                f"{barred} barred out of {count} total."
            )

        unbarred = count - barred

        pieces.extend(
            [rf"S_{index}^\dagger"] * barred
        )
        pieces.extend(
            [rf"S_{index}"] * unbarred
        )

    if len(pieces) != 2:
        raise ValueError(
            "Expected exactly two scalar legs in an EFT1 Wilson term, "
            f"found {len(pieces)}."
        )

    return r"\,".join(pieces)


def _matching_cg_tex(term: Mapping[str, object]) -> str:
    """Render the invariant contractions in one matched complex term."""
    names = [
        str(name)
        for name in (term.get("CGNames", []) or [])
        if str(name).strip()
    ]

    if not names:
        return r"\mathbf{1}"

    rendered: list[str] = []

    for name in names:
        if name == "eps[SU2L]":
            rendered.append(r"\epsilon_{SU(2)_L}")
        else:
            rendered.append(
                rf"\mathrm{{{latex_escape_text(name)}}}"
            )

    return r"\,".join(rendered)


def _matching_operator_rows(
    record: RunRecord,
) -> tuple[list[str], str]:
    """Summarise F-threshold operators in the original complex scalar basis.

    Matchete normally exports PL and PR terms which are conjugate descriptions
    of one Hermitian interaction.  The report therefore uses the PL term as
    the canonical representative and appends ``+ h.c.``.

    Distinct scalar channels and distinct invariant contractions remain
    separate.  Genuine S1-S1 and S2-S2 operators are therefore retained rather
    than being collapsed into the mixed S1-S2 channel.

    The real C_ijab representation is still constructed by the EFT1 tensor
    adapter when required by the Wilson-tensor RGE machinery; it is simply not
    expanded component-by-component in this report.
    """
    seed_path = (
        record.output_dir
        / "data"
        / "eft1_after_F_wilson_seed.json"
    )

    if not seed_path.is_file():
        return [], "No exported EFT1 Wilson seed was available for this run."

    import json

    seed = json.loads(
        seed_path.read_text(encoding="utf-8")
    )

    terms = list(
        seed.get("TreeWilsonTerms", []) or []
    )

    representatives = [
        term
        for term in terms
        if str(term.get("Chirality", "")).strip() == "PL"
    ]

    # Defensive fallback for a future Matchete export containing only one
    # chirality rather than an explicit PL/PR pair.
    if not representatives:
        representatives = terms

    rows: list[str] = []

    seen: set[
        tuple[
            str,
            str,
            tuple[str, ...],
        ]
    ] = set()

    for term in representatives:
        scalar_product = _matching_scalar_product_tex(term)

        prefactor = sp.simplify(
            _eft1_wilson_term_prefactor(term)
        )

        prefactor_tex = (
            sp.latex(prefactor)
            .replace("MF", r"M_F")
        )

        cg_names = tuple(
            str(name)
            for name in (term.get("CGNames", []) or [])
            if str(name).strip()
        )

        cg_tex = _matching_cg_tex(term)

        # Remove only exact duplicate descriptions.  Different scalar
        # structures, coefficients or CG contractions remain independent.
        key = (
            scalar_product,
            sp.srepr(prefactor),
            cg_names,
        )

        if key in seen:
            continue

        seen.add(key)

        operator_tex = (
            r"(L_i^T C L_j)\,"
            + scalar_product
            + r"+\mathrm{h.c.}"
        )

        rows.append(
            rf"$ {operator_tex} $"
            rf" & $ {prefactor_tex} $"
            rf" & $ {cg_tex} $ \\"
        )

    count = len(rows)

    noun = (
        "structure"
        if count == 1
        else "structures"
    )

    note = (
        f"{count} independent complex matched {noun} displayed. "
        "The conjugate chirality is represented by + h.c.; "
        "the internal real C_ijab expansion is not printed."
    )

    return rows, note


'@

    $text = Replace-Block `
        -Text $text `
        -StartMarker "def _matching_tensor_rows(" `
        -EndMarker "def _direct_weinberg_comparison(" `
        -Replacement $newMatchingFunctions

    Write-Host "[OK] Replaced _matching_tensor_rows."
}
elseif ($text.Contains("def _matching_operator_rows(")) {
    Write-Host "[INFO] Complex matching functions are already present."
}
else {
    throw @"
Neither the old nor new matching-row function was found.

Expected either:
    def _matching_tensor_rows(

or:
    def _matching_operator_rows(
"@
}


# ============================================================================
# 4. Replace _matching_tensor_appendix
# ============================================================================

if ($text.Contains("def _matching_tensor_appendix(")) {

    $newMatchingSection = @'
def _matching_operator_section(
    records: Sequence[RunRecord],
) -> list[str]:
    """Display F-threshold matching in the original complex scalar basis."""
    lines = [
        r"\section*{F-threshold matched operators}",
        r"The dimension-five matching is displayed in the original complex "
        r"scalar basis exported by Matchete.  A PL term and its PR conjugate "
        r"are shown once as one operator plus its Hermitian conjugate.",
        r"Distinct $S_1S_1$, $S_1S_2$ and $S_2S_2$ channels, when gauge "
        r"allowed, remain separate.  Their exact Yukawa/mass prefactors and "
        r"$SU(2)_L$ invariant contractions are retained.",
        r"The conversion to the global real tensor $C_{ijab}$ is still used "
        r"internally by the generic Wilson-tensor RGE calculation, but its "
        r"many basis-dependent real components are not printed here.",
    ]

    for record in records:
        rows, note = _matching_operator_rows(record)

        lines.append(
            rf"\subsection*{{$ {_model_column_label(record)} $}}"
        )

        lines.append(
            rf"\textit{{{latex_escape_text(note)}}}"
        )

        if not rows:
            continue

        lines.extend(
            [
                r"\begin{longtable}{@{}"
                r"p{0.36\linewidth}"
                r"p{0.25\linewidth}"
                r"p{0.29\linewidth}"
                r"@{}}",
                r"\toprule",
                r"complex matched operator"
                r" & matched prefactor"
                r" & $SU(2)_L$ invariant \\",
                r"\midrule",
                *rows,
                r"\bottomrule",
                r"\end{longtable}",
            ]
        )

    return lines


'@

    $text = Replace-Block `
        -Text $text `
        -StartMarker "def _matching_tensor_appendix(" `
        -EndMarker "def _generated_non_singlet_section(" `
        -Replacement $newMatchingSection

    Write-Host "[OK] Replaced _matching_tensor_appendix."
}
elseif ($text.Contains("def _matching_operator_section(")) {
    Write-Host "[INFO] Complex matching section is already present."
}
else {
    throw @"
Neither the old nor new matching section was found.

Expected either:
    def _matching_tensor_appendix(

or:
    def _matching_operator_section(
"@
}


# ============================================================================
# 5. Optional active-theory wording
# ============================================================================

$oldTheory = @'
            r"Active theory: $\mathrm{SM}+S_1+S_2+C_{LLS_1S_2}$.  "
            r"The fermion $F$ is absent, so $y_1$, $y_2$ and $M_F$ are not "
            r"running EFT couplings."
'@

$newTheory = @'
            r"Active theory: $\mathrm{SM}+S_1+S_2$ with matched "
            r"dimension-five $LLSS$ operators.  The fermion $F$ is absent, "
            r"so $y_1$, $y_2$ and $M_F$ are matching parameters rather than "
            r"running EFT couplings."
'@

if ($text.Contains($oldTheory)) {
    $text = $text.Replace($oldTheory, $newTheory)
    Write-Host "[OK] Updated F-first theory description."
}
else {
    Write-Host "[INFO] F-first theory wording differs locally; skipped wording update."
}


# ============================================================================
# 6. Change report call
# ============================================================================

$oldCall = "    lines.extend(_matching_tensor_appendix(records))"
$newCall = "    lines.extend(_matching_operator_section(records))"

if ($text.Contains($oldCall)) {
    $text = $text.Replace(
        $oldCall,
        $newCall
    )

    Write-Host "[OK] Updated F-first report call."
}
elseif ($text.Contains($newCall)) {
    Write-Host "[INFO] F-first report already calls complex matching section."
}
else {
    throw @"
Could not find the matching-section call.

Expected either:
$oldCall

or:
$newCall
"@
}


# ============================================================================
# 7. Sanity checks
# ============================================================================

$failures = @()

if ($text.Contains("def _matching_tensor_rows(")) {
    $failures += "_matching_tensor_rows still exists"
}

if ($text.Contains("def _matching_tensor_appendix(")) {
    $failures += "_matching_tensor_appendix still exists"
}

if (-not $text.Contains("def _matching_scalar_product_tex(")) {
    $failures += "_matching_scalar_product_tex is missing"
}

if (-not $text.Contains("def _matching_operator_rows(")) {
    $failures += "_matching_operator_rows is missing"
}

if (-not $text.Contains("def _matching_operator_section(")) {
    $failures += "_matching_operator_section is missing"
}

if (-not $text.Contains("_eft1_wilson_term_prefactor")) {
    $failures += "_eft1_wilson_term_prefactor import is missing"
}

if (-not $text.Contains("_matching_operator_section(records)")) {
    $failures += "F-first report does not call _matching_operator_section"
}

if ($failures.Count -gt 0) {
    Write-Host ""
    Write-Host "Sanity-check failures:"

    foreach ($failure in $failures) {
        Write-Host "  - $failure"
    }

    throw "Generated source failed sanity checks."
}


# ============================================================================
# 8. Write COMPLETE generated replacement
# ============================================================================

Set-Content `
    -Path $Generated `
    -Value $text `
    -Encoding UTF8

Write-Host ""
Write-Host "============================================================"
Write-Host "FULL REPLACEMENT GENERATED"
Write-Host "============================================================"
Write-Host ""
Write-Host "Generated file:"
Write-Host "  $Generated"
Write-Host ""
Write-Host "Intended repository path:"
Write-Host "  Reports\GroupFactorReports.py"
Write-Host ""
Write-Host "Backup:"
Write-Host "  $Backup"
Write-Host ""
Write-Host "The live GroupFactorReports.py has NOT been overwritten."
Write-Host ""