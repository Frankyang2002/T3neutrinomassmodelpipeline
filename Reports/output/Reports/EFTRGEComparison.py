from __future__ import annotations

"""Across-model report for Weinberg matching and its SMEFT RGE."""

from pathlib import Path

from common.Paths import REPORT_OUTPUT_DIR
from common.Records import RunRecord
from Reports.ReportGeneration import (
    compile_latex_document,
    latex_escape_text,
    latex_fraction,
    latex_model_heading,
    normalise_physics_latex,
    paper_notation_key_lines,
    record_quantum_numbers,
    split_latex_terms,
)


def _status(value: object) -> str:
    if value is True or value == "Success":
        return r"\textbf{Success}"
    if value in (False, "Failed"):
        return r"\textbf{Failed}"
    if value in (None, "", "NotRun"):
        return r"Not run"
    return latex_escape_text(value)


def _coefficient_block(latex: str) -> list[str]:
    """Render C5 as additive lines so long matching expressions remain readable."""
    latex = normalise_physics_latex(latex)
    terms = split_latex_terms(latex.strip()) if latex else []
    if not terms:
        return [r"\textit{Matched $C_5$ LaTeX was not available for this model.}"]

    lines: list[str] = []
    for index, term in enumerate(terms):
        prefix = r"C_5(M)=" if index == 0 else r"\phantom{C_5(M)=}"
        lines.extend(
            [
                r"\noindent\adjustbox{max width=\linewidth}{$\displaystyle "
                + prefix
                + term
                + r"$}\par",
                r"\smallskip",
            ]
        )
    return lines


def _model_status_table(records: list[RunRecord]) -> list[str]:
    lines = [
        r"\begin{center}",
        r"\begin{tabular}{lccc}",
        r"\toprule",
        r"Model & matched $C_5$ & one-generation RGE check & full-flavor RGE \\",
        r"\midrule",
    ]

    for record in records:
        summary = record.summary
        c5_ok = summary.get("WeinbergExtractionStatus")
        lines.append(
            " & ".join(
                [
                    latex_escape_text(record.name),
                    _status(c5_ok),
                    _status(summary.get("RGEStatus")),
                    _status(summary.get("FlavorRGEStatus")),
                ]
            )
            + r" \\" 
        )

    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{center}"])
    return lines


def _matching_summary_table(records: list[RunRecord]) -> list[str]:
    lines = [
        r"\begin{longtable}{@{}lcccccc@{}}",
        r"\toprule",
        r"Model & $d_{S_1}$ & $Y_{S_1}$ & $d_{S_2}$ & $Y_{S_2}$ & $d_F$ & $Y_F$ \\",
        r"\midrule",
        r"\endfirsthead",
        r"\toprule",
        r"Model & $d_{S_1}$ & $Y_{S_1}$ & $d_{S_2}$ & $Y_{S_2}$ & $d_F$ & $Y_F$ \\",
        r"\midrule",
        r"\endhead",
    ]
    for record in records:
        d_s1, y_s1, d_s2, y_s2, d_f, y_f = record_quantum_numbers(record)
        lines.append(
            " & ".join(
                [
                    latex_escape_text(record.name),
                    rf"${d_s1}$",
                    rf"${latex_fraction(y_s1)}$",
                    rf"${d_s2}$",
                    rf"${latex_fraction(y_s2)}$",
                    rf"${d_f}$",
                    rf"${latex_fraction(y_f)}$",
                ]
            )
            + r" \\" 
        )
        lines.append(r"\midrule")
    lines.extend([r"\bottomrule", r"\end{longtable}"])
    return lines


def write_weinberg_rge_comparison(
    records: list[RunRecord],
    output_path: Path | None = None,
) -> Path:
    """Write the across-model Weinberg matching/RGE comparison report."""
    if output_path is None:
        output_path = REPORT_OUTPUT_DIR / "weinberg_rge_comparison.tex"

    lines: list[str] = [
        r"\documentclass[10pt]{article}",
        r"\usepackage[margin=1.6cm]{geometry}",
        r"\usepackage{amsmath,amssymb,adjustbox,longtable,array,booktabs}",
        r"\usepackage[T1]{fontenc}",
        r"\allowdisplaybreaks[4]",
        r"\setlength{\emergencystretch}{3em}",
        r"\begin{document}",
        r"\section*{T3 Weinberg-operator RGE comparison}",
        *paper_notation_key_lines(),
        (
            r"We have Weinberg operator in SMEFT after integrating out our 2 fields."
        ),
        r"\begin{equation}",
        r"\mathcal L_{\rm EFT}\supset \frac12(C_5)_{ij}(L_i\!\cdot H)(L_j\!\cdot H)+\mathrm{h.c.}",
        r"\end{equation}",
        r"\section{Run status}",
        *_model_status_table(records),
        r"\section{Model assignments}",
        *_matching_summary_table(records),
        r"\section{Matched Weinberg coefficients}",
        (
            r"These are the model-dependent boundary conditions for the EFT.  They are not "
            r"five different EFT beta functions: they are five different UV matching results "
            r"feeding the same low-energy SMEFT evolution."
        ),
    ]

    for index, record in enumerate(records):
        if index:
            lines.append(r"\medskip")
        lines.extend(
            [
                rf"\subsection*{{{latex_escape_text(record.name)}}}",
                rf"\noindent ${latex_model_heading(record)}$\par\medskip",
                *_coefficient_block(record.summary.get("WeinbergCoefficientLaTeX", "")),
            ]
        )

    lines.extend(
        [
            r"\clearpage",
            r"\section{Three-generation flavor RGE}",
            r"Writing $K\equiv C_5$,",
            r"\begin{align}",
            r"16\pi^2\frac{dK}{d\ln\mu}",
            r"&=(2\lambda_1-3g_2^2+2T)K",
            r"-\frac32\left[Y_eY_e^\dagger K+K(Y_eY_e^\dagger)^T\right],\\",
            r"T&=\operatorname{Tr}\!\left(Y_eY_e^\dagger+3Y_uY_u^\dagger+3Y_dY_d^\dagger\right).",
            r"\end{align}",
            (
                r"This matrix equation is common to T3-A--E once the heavy T3 particles are "
                r"removed.  Representation dependence survives through the matched boundary "
                r"condition and, in a threshold treatment with non-degenerate masses, through "
                r"where individual heavy particles are integrated out."
            ),
            r"\end{document}",
            "",
        ]
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines), encoding="utf-8", newline="\n")
    print(f"\nWeinberg RGE comparison report:\n{output_path}")
    return output_path


def write_and_compile_weinberg_rge_comparison(records: list[RunRecord]) -> Path:
    report_tex = write_weinberg_rge_comparison(records)
    compile_latex_document(report_tex)
    return report_tex
