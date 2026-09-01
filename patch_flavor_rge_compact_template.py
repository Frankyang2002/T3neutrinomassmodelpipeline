from pathlib import Path

stage_path = Path("FlavorMatchedRGEStage.py")
beta_path = Path("SMEFTWeinbergFlavorRGE.py")

if not stage_path.exists():
    raise FileNotFoundError(stage_path)

if not beta_path.exists():
    raise FileNotFoundError(beta_path)

# ---------------------------------------------------------------------------
# 1. Make beta_weinberg_matrix optionally stay compact.
# ---------------------------------------------------------------------------

beta_text = beta_path.read_text(encoding="utf-8")

old_signature = """def beta_weinberg_matrix(
    K: sp.MatrixBase,
    Ye: sp.MatrixBase,
    Yu: sp.MatrixBase,
    Yd: sp.MatrixBase,
    *,
    higgs_quartic: sp.Expr = lambdaH,
    weak_coupling: sp.Expr = g2,
) -> sp.Matrix:
"""

new_signature = """def beta_weinberg_matrix(
    K: sp.MatrixBase,
    Ye: sp.MatrixBase,
    Yu: sp.MatrixBase,
    Yd: sp.MatrixBase,
    *,
    higgs_quartic: sp.Expr = lambdaH,
    weak_coupling: sp.Expr = g2,
    expand_result: bool = True,
) -> sp.Matrix:
"""

if "expand_result: bool" not in beta_text:
    if old_signature not in beta_text:
        raise RuntimeError(
            "Could not locate beta_weinberg_matrix signature."
        )

    beta_text = beta_text.replace(
        old_signature,
        new_signature,
        1,
    )

old_return = "    return beta.applyfunc(sp.expand)\n"
new_return = """    if expand_result:
        return beta.applyfunc(sp.expand)

    return beta
"""

if old_return in beta_text:
    beta_text = beta_text.replace(
        old_return,
        new_return,
        1,
    )

beta_path.write_text(
    beta_text,
    encoding="utf-8",
)

# ---------------------------------------------------------------------------
# 2. Keep the matched C5 matrix as its own output, but write the SMEFT beta
#    as a compact universal matrix template in symbolic K_ij entries.
#
#    This avoids substituting the potentially very large T3 loop expression
#    into every term of the 3x3 flavor RGE.  No physics is lost:
#
#       matched C5 matrix  +  beta[K]
#
#    completely specifies the boundary condition and its RGE.
# ---------------------------------------------------------------------------

stage_text = stage_path.read_text(encoding="utf-8")

old_import = (
    "from SMEFTWeinbergFlavorRGE import "
    "beta_weinberg_matrix, symbolic_complex_matrix\n"
)

new_import = (
    "from SMEFTWeinbergFlavorRGE import (\n"
    "    beta_weinberg_matrix,\n"
    "    symbolic_complex_matrix,\n"
    "    symbolic_symmetric_matrix,\n"
    ")\n"
)

if "symbolic_symmetric_matrix" not in stage_text:
    if old_import not in stage_text:
        raise RuntimeError(
            "Could not locate SMEFTWeinbergFlavorRGE import."
        )

    stage_text = stage_text.replace(
        old_import,
        new_import,
        1,
    )

old_block = """    K = flavor_result["K"]

    Ye = symbolic_complex_matrix("Ye", 3, 3)
    Yu = symbolic_complex_matrix("Yu", 3, 3)
    Yd = symbolic_complex_matrix("Yd", 3, 3)

    beta_K = beta_weinberg_matrix(
        K,
        Ye,
        Yu,
        Yd,
    )

    symmetry_difference = (
        beta_K - beta_K.T
    ).applyfunc(sp.simplify)

    if any(entry != 0 for entry in symmetry_difference):
        raise RuntimeError(
            "Full-flavor beta_C5 is not symmetric."
        )
"""

new_block = """    # The matched T3 matrix is already written above.  Keep the SMEFT RGE
    # itself in a compact universal K_ij basis instead of inserting the full
    # (potentially very large) loop expression into every beta-function term.
    K = symbolic_symmetric_matrix("K", 3)

    Ye = symbolic_complex_matrix("Ye", 3, 3)
    Yu = symbolic_complex_matrix("Yu", 3, 3)
    Yd = symbolic_complex_matrix("Yd", 3, 3)

    beta_K = beta_weinberg_matrix(
        K,
        Ye,
        Yu,
        Yd,
        expand_result=False,
    )

    # K is symmetric by construction and beta_K has the form
    # A K + K A^T plus a scalar multiple of K, so symmetry is exact.
"""

if old_block not in stage_text:
    # Handle versions where expand_result=False was already added but the
    # expensive symmetry simplification remained.
    start = stage_text.find('    K = flavor_result["K"]\n')
    end = stage_text.find(
        '    beta_path = output_dir / "c5_flavor_beta_matrix.txt"\n'
    )

    if start == -1 or end == -1 or end <= start:
        raise RuntimeError(
            "Could not locate the full-flavor beta construction block."
        )

    stage_text = (
        stage_text[:start]
        + new_block
        + "\n"
        + stage_text[end:]
    )
else:
    stage_text = stage_text.replace(
        old_block,
        new_block,
        1,
    )

# Add explicit metadata so the output's meaning is clear.
old_summary_line = '        "FlavorBetaC5Symmetric": True,\n'
new_summary_lines = (
    '        "FlavorBetaC5Symmetric": True,\n'
    '        "FlavorBetaRepresentation": '
    '"compact universal SMEFT beta[K]; matched K stored separately",\n'
)

if (
    old_summary_line in stage_text
    and "FlavorBetaRepresentation" not in stage_text
):
    stage_text = stage_text.replace(
        old_summary_line,
        new_summary_lines,
        1,
    )

stage_path.write_text(
    stage_text,
    encoding="utf-8",
)

print("Patched:", beta_path.resolve())
print("Patched:", stage_path.resolve())
print("Full-flavor RGE now uses compact symbolic K_ij template.")
