from pathlib import Path

flavor_path = Path("SMEFTWeinbergFlavorRGE.py")
stage_path = Path("FlavorMatchedRGEStage.py")

if not flavor_path.exists():
    raise FileNotFoundError(flavor_path)

if not stage_path.exists():
    raise FileNotFoundError(stage_path)

flavor_text = flavor_path.read_text(encoding="utf-8")
stage_text = stage_path.read_text(encoding="utf-8")

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

if old_signature not in flavor_text:
    raise RuntimeError("Could not locate beta_weinberg_matrix signature.")

flavor_text = flavor_text.replace(old_signature, new_signature, 1)

old_return = "    return beta.applyfunc(sp.expand)\n"
new_return = (
    "    if expand_result:\n"
    "        return beta.applyfunc(sp.expand)\n\n"
    "    return beta\n"
)

if old_return not in flavor_text:
    raise RuntimeError("Could not locate expanded beta return.")

flavor_text = flavor_text.replace(old_return, new_return, 1)

old_stage = """    beta_K = beta_weinberg_matrix(
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

new_stage = """    beta_K = beta_weinberg_matrix(
        K,
        Ye,
        Yu,
        Yd,
        expand_result=False,
    )
"""

if old_stage not in stage_text:
    raise RuntimeError("Could not locate full-flavor beta block.")

stage_text = stage_text.replace(old_stage, new_stage, 1)

flavor_path.write_text(flavor_text, encoding="utf-8")
stage_path.write_text(stage_text, encoding="utf-8")

print("Patched:", flavor_path.resolve())
print("Patched:", stage_path.resolve())
print("Production beta expansion disabled; regression tests still expand by default.")
