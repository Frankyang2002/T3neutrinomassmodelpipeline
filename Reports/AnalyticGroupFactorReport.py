from __future__ import annotations

"""Generate the clean analytic T3 group-factor RGE report.

This is a presentation-layer report: it does not parse RGBeta LaTeX and it
does not use coefficient-insensitive term signatures.  The physics formulas
are written explicitly in terms of representation data and validated group
factors.

Default output:
    output/reports/RGE/AnalyticGroupFactors.tex

Optional:
    --compile
runs pdflatex twice in the output directory.
"""

import argparse
from fractions import Fraction
from pathlib import Path
import subprocess
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from common.Paths import REPORT_OUTPUT_DIR
from common.T3Model import T3_CLASSES
from Reports.ReportGeneration import paper_notation_key_lines

from RGE.group_factors.ScalarMassGroupFactors import (
    scalar_mass_group_factors,
)
from RGE.group_factors.DirectWeinbergGroupFactors import (
    direct_weinberg_group_factor,
)



def c2(d: int) -> Fraction:
    return Fraction(d * d - 1, 4)


def leg_factors(dS: int, dF: int):
    dgt = max(dS, dF)
    return (
        Fraction(dgt, dS),
        Fraction(dgt, 2),
        Fraction(dgt, dF),
    )


def ftex(x: Fraction | str) -> str:
    """Return a fraction/integer as LaTeX suitable inside math mode."""
    if isinstance(x, str):
        x = Fraction(x)
    if x.denominator == 1:
        return str(x.numerator)
    return rf"\frac{{{x.numerator}}}{{{x.denominator}}}"


def mathcell(value: str) -> str:
    """Render arbitrary LaTeX safely in a table cell without relying on $...$."""
    return rf"\ensuremath{{{value}}}"


def sympy_to_latex_text(text: str) -> str:
    """Small exact formatter for the reduced Wigner-6j factors."""
    import sympy as sp
    expr = sp.sympify(text, locals={"sqrt": sp.sqrt})
    return sp.latex(expr)


def model_table_rows() -> list[str]:
    rows = []
    for model, (d1, d2, dF) in T3_CLASSES.items():
        GS1, GL1, GF1 = leg_factors(d1, dF)
        GS2, GL2, GF2 = leg_factors(d2, dF)

        y1_self = Fraction(1, 2) * (GL1 + GF1)
        y1_cross = Fraction(1, 2) * GL2
        y2_self = Fraction(1, 2) * (GL2 + GF2)
        y2_cross = Fraction(1, 2) * GL1

        rw = direct_weinberg_group_factor(d1, d2, dF)

        rows.append(
            " & ".join(
                [
                    rf"T3-{model}",
                    str(d1),
                    str(d2),
                    str(dF),
                    mathcell(ftex(GS1)),
                    mathcell(ftex(y1_self)),
                    mathcell(ftex(y1_cross)),
                    mathcell(ftex(GS2)),
                    mathcell(ftex(y2_self)),
                    mathcell(ftex(y2_cross)),
                    mathcell(sympy_to_latex_text(rw.reduced_factor)),
                ]
            )
            + r" \\"
        )
    return rows


def scalar_mass_rows(alpha: int = 0) -> list[str]:
    rows = []
    for model, (d1, d2, dF) in T3_CLASSES.items():
        gf = scalar_mass_group_factors(d1, d2, dF, alpha)
        b1 = gf.beta_mS1Sq
        b2 = gf.beta_mS2Sq
        rows.append(
            " & ".join(
                [
                    rf"T3-{model}",
                    mathcell(b1["mS1Sq*Tr_y1"]),
                    mathcell(b1["Tr_MF_y1"]),
                    mathcell(b1["lambdaS1*mS1Sq"]),
                    mathcell(b1["lambda12*mS2Sq"]),
                    mathcell(b2["mS2Sq*Tr_y2"]),
                    mathcell(b2["Tr_MF_y2"]),
                    mathcell(b2["lambdaS2*mS2Sq"]),
                    mathcell(b2["lambda12*mS1Sq"]),
                ]
            )
            + r" \\"
        )
    return rows



def document() -> str:
    lines = [
        r"\documentclass[11pt]{article}",
        r"\usepackage[a4paper,margin=1.8cm]{geometry}",
        r"\usepackage{amsmath,amssymb,booktabs,longtable,array}",
        r"\usepackage[T1]{fontenc}",
        r"\setlength{\parindent}{0pt}",
        r"\setlength{\parskip}{0.6em}",
        r"\begin{document}",
        r"\section*{T3 analytic one-loop group factors}",
        *paper_notation_key_lines(),
        (
            r"This report presents the validated analytic representation-dependent "
            r"group factors used by the T3 RGE pipeline.  Model dimensions are "
            r"read from \texttt{common.T3Model.T3\_CLASSES} and model-specific coefficients "
            r"are evaluated from the analytic group-factor implementations, not "
            r"copied from RGBeta comparison tables."
        ),
        r"\subsection*{Representation conventions}",
        r"""\[
j=\frac{d-1}{2},\qquad
C_2(d)=\frac{d^2-1}{4},\qquad
T(d)=\frac{d(d^2-1)}{12}.
\]""",
        r"For each Yukawa vertex define",
        r"""\[
G_S=\frac{d_>}{d_S},\qquad
G_L=\frac{d_>}{2},\qquad
G_F=\frac{d_>}{d_F},\qquad
d_>=\max(d_S,d_F).
\]""",
        r"\subsection*{Yukawa beta functions}",
        r"""\[
\begin{aligned}
16\pi^2\beta_{y_1}={}&
G_{S_1}\operatorname{Tr}(y_1y_1^\dagger)y_1
+\frac12(G_{L_1}+G_{F_1})y_1y_1^\dagger y_1\\
&+\frac12G_{L_2}y_2y_2^\dagger y_1
+\frac12Y_eY_e^\dagger y_1\\
&-3g_2^2[C_2(L)+C_2(F)]y_1
-3g_1^2[Y_L^2+Y_F^2]y_1 .
\end{aligned}
\]""",
        r"""\[
\begin{aligned}
16\pi^2\beta_{y_2}={}&
G_{S_2}\operatorname{Tr}(y_2y_2^\dagger)y_2
+\frac12(G_{L_2}+G_{F_2})y_2y_2^\dagger y_2\\
&+\frac12G_{L_1}y_1y_1^\dagger y_2
+\frac12Y_eY_e^\dagger y_2\\
&-3g_2^2[C_2(L)+C_2(F)]y_2
-3g_1^2[Y_L^2+Y_F^2]y_2\\
&+\delta_{\rm SC}\left[
G_{S_2}\operatorname{Tr}(y_2y_1^\dagger)y_1
+\frac12G_{F_2}y_2y_1^\dagger y_1
\right].
\end{aligned}
\]""",
        r"Here $\delta_{\rm SC}=1$ only for the physical self-conjugate branch "
        r"$Y_F=0$ with odd $d_F$; otherwise it is zero.",
        r"\subsection*{Fermion mass}",
        r"""\[
16\pi^2\beta_{M_F}
=
\frac12G_{F_1}(y_1^Ty_1^*)M_F
+\frac12G_{F_2}M_F(y_2^\dagger y_2)
-6g_2^2C_2(F)M_F
-6g_1^2Y_F^2M_F .
\]""",
        r"\subsection*{Heavy scalar masses}",
        r"""\[
\begin{aligned}
16\pi^2\beta_{m_1^2}={}&
2G_{S_1}m_1^2\operatorname{Tr}(y_1y_1^\dagger)
+\Xi_F\,G_{S_1}\operatorname{Tr}(M_FM_F^\dagger y_1^Ty_1^*)\\
&+2(d_{S_1}+1)\lambda_{S_1}^{(1)}m_1^2
-6C_2(S_1)g_2^2m_1^2
-6Y_{S_1}^2g_1^2m_1^2\\
&+2d_{S_2}\lambda_{12}^{(1)}m_2^2+\cdots ,
\end{aligned}
\]""",
        r"with $\Xi_F=-4$ on the vector-like branch and $\Xi_F=-16$ on the "
        r"physical self-conjugate branch.  The $S_2$ equation follows by "
        r"$1\leftrightarrow2$.",
        r"\subsection*{Direct Weinberg mixing}",
        r"""\[
16\pi^2\beta_\kappa=R_W\,\lambda_5\,C_{12},
\]""",
        r"""\[
R_W=
\frac43\,\eta_F\sqrt{3d_{S_1}d_{S_2}d_F}
\begin{Bmatrix}
\frac12&\frac12&1\\
j_{S_2}&j_{S_1}&j_F
\end{Bmatrix},
\qquad
\eta_F=(+1,+1,-1)\;\text{for}\;d_F=(1,2,3).
\]""",
        r"\subsection*{Representation summary}",
        r"\begin{center}",
        r"\small",
        r"\begin{longtable}{@{}lcccccccccc@{}}",
        r"\toprule",
        r"Model & $d_{S_1}$ & $d_{S_2}$ & $d_F$ "
        r"& $G_{S_1}$ & $c_{11}$ & $c_{21}$ "
        r"& $G_{S_2}$ & $c_{22}$ & $c_{12}$ & $R_W$ \\",
        r"\midrule",
        r"\endfirsthead",
        r"\toprule",
        r"Model & $d_{S_1}$ & $d_{S_2}$ & $d_F$ "
        r"& $G_{S_1}$ & $c_{11}$ & $c_{21}$ "
        r"& $G_{S_2}$ & $c_{22}$ & $c_{12}$ & $R_W$ \\",
        r"\midrule",
        r"\endhead",
        *model_table_rows(),
        r"\bottomrule",
        r"\end{longtable}",
        r"\end{center}",
        (
            r"$c_{11}=\frac12(G_{L_1}+G_{F_1})$, "
            r"$c_{21}=\frac12G_{L_2}$, "
            r"$c_{22}=\frac12(G_{L_2}+G_{F_2})$, "
            r"$c_{12}=\frac12G_{L_1}$."
        ),
        r"\subsection*{Scalar-mass coefficient check at $\alpha=0$}",
        r"\begin{center}",
        r"\small",
        r"\begin{longtable}{@{}lrrrrrrrr@{}}",
        r"\toprule",
        r"Model & $2G_{S_1}$ & $-4G_{S_1}$ & $2(d_{S_1}+1)$ & $2d_{S_2}$ "
        r"& $2G_{S_2}$ & $-4G_{S_2}$ & $2(d_{S_2}+1)$ & $2d_{S_1}$ \\",
        r"\midrule",
        r"\endfirsthead",
        r"\toprule",
        r"Model & $2G_{S_1}$ & $-4G_{S_1}$ & $2(d_{S_1}+1)$ & $2d_{S_2}$ "
        r"& $2G_{S_2}$ & $-4G_{S_2}$ & $2(d_{S_2}+1)$ & $2d_{S_1}$ \\",
        r"\midrule",
        r"\endhead",
        *scalar_mass_rows(alpha=0),
        r"\bottomrule",
        r"\end{longtable}",
        r"\end{center}",
        r"\end{document}",
        "",
    ]
    return "\n".join(lines)


def compile_tex(tex_path: Path) -> None:
    for _ in range(2):
        subprocess.run(
            [
                "pdflatex",
                "-interaction=nonstopmode",
                "-halt-on-error",
                tex_path.name,
            ],
            cwd=tex_path.parent,
            check=True,
        )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=REPORT_OUTPUT_DIR / "RGE" / "AnalyticGroupFactors.tex",
    )
    parser.add_argument("--compile", action="store_true")
    args = parser.parse_args()

    print(f"Report output root: {REPORT_OUTPUT_DIR}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(document(), encoding="utf-8")
    print(f"Wrote {args.output}")

    if args.compile:
        compile_tex(args.output)
        print(f"Wrote {args.output.with_suffix('.pdf')}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
