from __future__ import annotations

"""Report the renormalisable SM RGEs used below the common T3 threshold."""

from pathlib import Path

from common.Paths import REPORT_OUTPUT_DIR
from common.Records import RunRecord
from Reports.ReportGeneration import compile_latex_document, latex_escape_text


def _model_universality_table(records: list[RunRecord]) -> list[str]:
    lines = [
        r"\begin{center}",
        r"\begin{tabular}{lc}",
        r"\toprule",
        r"UV model & Renormalisable EFT below the common threshold \\",
        r"\midrule",
    ]
    for record in records:
        lines.append(
            rf"{latex_escape_text(record.name)} & Standard Model renormalisable sector \\" 
        )
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{center}"])
    return lines


def write_eft_renormalisable_rge(
    records: list[RunRecord],
    output_path: Path | None = None,
) -> Path:
    """Write the common one-loop renormalisable EFT RGE report below M."""
    if output_path is None:
        output_path = REPORT_OUTPUT_DIR / "eft_renormalisable_rge.tex"

    lines: list[str] = [
        r"\documentclass[10pt]{article}",
        r"\usepackage[margin=1.7cm]{geometry}",
        r"\usepackage{amsmath,amssymb,booktabs,array}",
        r"\usepackage[T1]{fontenc}",
        r"\allowdisplaybreaks[4]",
        r"\setlength{\emergencystretch}{3em}",
        r"\begin{document}",
        r"\section*{T3 EFT renormalisable one-loop RGEs}",
        (
            r"The heavy T3 fields $S_1,S_2,F$ have been integrated out at the common matching scale $M$. "
            r"At the dimension-four level the active theory is the Standard Model, so these beta functions seems to be the same for T3-A--E."
        ),
        r"\section{Conventions}",
        r"All equations are written as",
        r"\begin{equation}",
        r"16\pi^2\,\mu\frac{dX}{d\mu}=\beta_X^{(1)}.",
        r"\end{equation}",
        r"The hypercharge convention is $Q=T_3+Y$, and the Higgs potential convention is",
        r"\begin{equation}",
        r"V(H)=\frac{\lambda_H}{2}(H^\dagger H)^2.",
        r"\end{equation}",
        r"Define",
        r"\begin{equation}",
        r"T=\operatorname{Tr}\!\left(Y_eY_e^\dagger+3Y_uY_u^\dagger+3Y_dY_d^\dagger\right).",
        r"\end{equation}",
        r"\section{Gauge couplings}",
        r"\begin{align}",
        r"\beta_{g_Y}^{(1)}&=\frac{41}{6}g_Y^3,\\",
        r"\beta_{g_2}^{(1)}&=-\frac{19}{6}g_2^3,\\",
        r"\beta_{g_3}^{(1)}&=-7g_3^3.",
        r"\end{align}",
        r"\section{Standard Model Yukawa matrices}",
        r"\begin{align}",
        r"\beta_{Y_u}^{(1)}={}&T Y_u"
        r"-\frac{3}{2}\left(Y_dY_d^\dagger Y_u-Y_uY_u^\dagger Y_u\right)"
        r"-\left(\frac{17}{12}g_Y^2+\frac94g_2^2+8g_3^2\right)Y_u,\\[1ex]",
        r"\beta_{Y_d}^{(1)}={}&T Y_d"
        r"+\frac{3}{2}\left(Y_dY_d^\dagger Y_d-Y_uY_u^\dagger Y_d\right)"
        r"-\left(\frac{5}{12}g_Y^2+\frac94g_2^2+8g_3^2\right)Y_d,\\[1ex]",
        r"\beta_{Y_e}^{(1)}={}&T Y_e"
        r"+\frac{3}{2}Y_eY_e^\dagger Y_e"
        r"-\left(\frac{15}{4}g_Y^2+\frac94g_2^2\right)Y_e.",
        r"\end{align}",
        (
            r"No more Yukawa couplings from coupling with heavy fields"
        ),
        r"\section{Higgs quartic}",
        r"With the project normalisation $V(H)=\lambda_H(H^\dagger H)^2/2$:",
        r"\begin{align}",
        r"\beta_{\lambda_H}^{(1)}={}&12\lambda_H^2"
        r"+\lambda_H\left[12\operatorname{Tr}(Y_dY_d^\dagger)"
        r"+4\operatorname{Tr}(Y_eY_e^\dagger)"
        r"+12\operatorname{Tr}(Y_uY_u^\dagger)-9g_2^2-3g_Y^2\right]\\",
        r"&-12\operatorname{Tr}(Y_dY_d^\dagger Y_dY_d^\dagger)"
        r"-4\operatorname{Tr}(Y_eY_e^\dagger Y_eY_e^\dagger)"
        r"-12\operatorname{Tr}(Y_uY_u^\dagger Y_uY_u^\dagger)\\",
        r"&+\frac94 g_2^4+\frac32 g_Y^2g_2^2+\frac34 g_Y^4.",
        r"\end{align}",
        (
            r"No more BSM scalar couplings, so less terms"
        ),
        r"\end{document}",
        "",
    ]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines), encoding="utf-8", newline="\n")
    print(f"\nEFT renormalisable RGE report:\n{output_path}")
    return output_path


def write_and_compile_eft_renormalisable_rge(records: list[RunRecord]) -> Path:
    report_tex = write_eft_renormalisable_rge(records)
    compile_latex_document(report_tex)
    return report_tex
