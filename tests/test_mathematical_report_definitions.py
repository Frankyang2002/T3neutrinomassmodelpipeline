from pathlib import Path
from scripts.AddMathematicalDefinitions import insert_definitions, MATH_START

def test_injects_and_reinjects_without_duplication():
    src=r"\documentclass{article}\begin{document}\section*{RGE}\end{document}"
    p=Path("RGE/UV.tex")
    out=insert_definitions(src,p)
    assert out.count(MATH_START)==1
    assert insert_definitions(out,p).count(MATH_START)==1
    assert r"C_2(d)&=\frac{d^2-1}{4}" in out

def test_c5_warns_about_integral_convention():
    out=insert_definitions(r"\begin{document}\end{document}",Path("Lagrangian/C5.tex"))
    assert "Unresolved matching convention" in out
    assert "I_3" in out

def test_group_factor_has_generator_algebra():
    out=insert_definitions(r"\begin{document}\end{document}",Path("GroupFactors/GF_UV.tex"))
    assert r"[T^a,T^b]" in out
