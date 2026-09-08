from __future__ import annotations

"""RGE report rendering and compilation utilities."""

import json
import math
from pathlib import Path
from typing import Any

from common.Records import RunRecord
from Reports.ReportGeneration import (
    compile_latex_document,
    latex_model_heading,
    report_output_dir_for,
)


def _latex_escape(value: object) -> str:
    """Escape ordinary text before placing it in a LaTeX document."""

    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }

    return "".join(
        replacements.get(character, character)
        for character in str(value)
    )


def _latex_status(value: object) -> str:
    """Return a compact colored status for the report tables."""

    if value is True or value == "Success":
        return r"\textcolor{green!45!black}{\textbf{Success}}"

    if value in (False, "Failed"):
        return r"\textcolor{red!70!black}{\textbf{Failed}}"

    if value in (None, "", "NotRun"):
        return r"\textcolor{black!55}{Not run}"

    return _latex_escape(value)


def _split_latex_terms(latex: str) -> list[str]:
    """Split a long LaTeX sum only at top-level plus and minus signs."""

    text = latex.strip()

    if not text:
        return []

    terms: list[str] = []
    current: list[str] = []
    brace_depth = 0
    parenthesis_depth = 0
    bracket_depth = 0
    escaped = False

    for index, character in enumerate(text):
        if escaped:
            current.append(character)
            escaped = False
            continue

        if character == "\\":
            current.append(character)
            escaped = True
            continue

        if character == "{":
            brace_depth += 1
        elif character == "}":
            brace_depth = max(0, brace_depth - 1)
        elif character == "(":
            parenthesis_depth += 1
        elif character == ")":
            parenthesis_depth = max(0, parenthesis_depth - 1)
        elif character == "[":
            bracket_depth += 1
        elif character == "]":
            bracket_depth = max(0, bracket_depth - 1)

        top_level = (
            brace_depth == 0
            and parenthesis_depth == 0
            and bracket_depth == 0
        )
        previous = text[index - 1] if index else ""

        if (
            character in "+-"
            and top_level
            and current
            and previous not in "^_eE"
        ):
            term = "".join(current).strip()

            if term:
                terms.append(term)

            current = [character]
        else:
            current.append(character)

    final_term = "".join(current).strip()

    if final_term:
        terms.append(final_term)

    return terms


def _coefficient_block(coefficient_latex: str) -> list[str]:
    """Render a long matched coefficient as readable additive lines."""

    terms = _split_latex_terms(coefficient_latex)

    if not terms:
        return [r"\textit{The matched coefficient was not available as LaTeX.}"]

    lines: list[str] = []

    for index, term in enumerate(terms):
        prefix = r"C_5=" if index == 0 else r"\phantom{C_5=}"
        lines.extend(
            [
                r"\noindent\adjustbox{max width=\linewidth}"
                r"{$\displaystyle "
                + prefix
                + term
                + r"$}\par",
                r"\smallskip",
            ]
        )

    return lines


def _scientific_latex(value: float, *, digits: int = 5) -> str:
    """Format one real number as compact LaTeX scientific notation."""

    value = float(value)

    if not math.isfinite(value):
        return _latex_escape(value)

    if value == 0.0:
        return "0"

    exponent = int(math.floor(math.log10(abs(value))))
    mantissa = value / (10.0**exponent)

    if -2 <= exponent <= 2:
        return f"{value:.{digits}g}"

    return rf"{mantissa:.{digits}g}\times 10^{{{exponent}}}"


def _complex_latex(real: float, imag: float) -> str:
    """Format a complex matrix entry without hiding its phase."""

    real = float(real)
    imag = float(imag)
    scale = max(abs(real), abs(imag), 1.0e-300)
    tolerance = 1.0e-12 * scale

    if abs(imag) <= tolerance:
        return _scientific_latex(real)

    if abs(real) <= tolerance:
        return _scientific_latex(imag) + r"\,i"

    sign = "+" if imag >= 0 else "-"
    return (
        _scientific_latex(real)
        + rf"\;{sign}\;"
        + _scientific_latex(abs(imag))
        + r"\,i"
    )


def _read_complex_matrix(
    output_dir: Path,
    relative_path: object,
) -> tuple[list[list[float]], list[list[float]]] | None:
    """Read a numerical complex-matrix JSON file when it is available."""

    if not relative_path:
        return None

    path = output_dir / str(relative_path)

    if not path.is_file():
        return None

    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
        real = payload["real"]
        imag = payload["imag"]
    except (OSError, KeyError, TypeError, json.JSONDecodeError):
        return None

    if (
        not isinstance(real, list)
        or not isinstance(imag, list)
        or len(real) != len(imag)
    ):
        return None

    return real, imag


def _matrix_block(
    title: str,
    matrix: tuple[list[list[float]], list[list[float]]] | None,
    *,
    symbol: str,
    scale: float = 1.0,
    units: str = "",
) -> list[str]:
    """Create a readable LaTeX block for one numerical complex matrix."""

    lines = [rf"\paragraph{{{title}}}"]

    if matrix is None:
        lines.append(r"\textit{Numerical matrix not available.}")
        return lines

    real, imag = matrix

    if not real or any(len(row) != len(real) for row in real):
        lines.append(r"\textit{Numerical matrix has an unsupported shape.}")
        return lines

    entries: list[str] = []

    for real_row, imag_row in zip(real, imag):
        if len(real_row) != len(imag_row):
            lines.append(r"\textit{Numerical matrix parts have different shapes.}")
            return lines

        entries.append(
            " & ".join(
                _complex_latex(scale * real_value, scale * imag_value)
                for real_value, imag_value in zip(real_row, imag_row)
            )
        )

    matrix_body = r" \\ ".join(entries)
    unit_suffix = rf"\;{units}" if units else ""

    lines.extend(
        [
            r"\begin{equation*}",
            r"\adjustbox{max width=\linewidth}{$\displaystyle",
            rf"{symbol}=\begin{{pmatrix}}{matrix_body}\end{{pmatrix}}{unit_suffix}",
            r"$}",
            r"\end{equation*}",
        ]
    )
    return lines


def _status_table(summary: dict[str, Any]) -> list[str]:
    """Build the stage-status table."""

    return [
        r"\begin{center}",
        r"\begin{tabular}{>{\bfseries}l l >{\bfseries}l l}",
        r"\toprule",
        (
            r"Matched EFT RGE & "
            + _latex_status(summary.get("RGEStatus"))
            + r" & Flavor RGE & "
            + _latex_status(summary.get("FlavorRGEStatus"))
            + r" \\"
        ),
        (
            r"Neutrino mass & "
            + _latex_status(summary.get("NeutrinoMassStatus"))
            + r" & Numerical running & "
            + _latex_status(summary.get("NumericalRGEStatus"))
            + r" \\"
        ),
        (
            r"Observables & "
            + _latex_status(summary.get("NeutrinoObservableStatus"))
            + r" & Internal benchmark & "
            + (
                _latex_status("Success")
                if summary.get("RGEStatus") == "Success"
                else _latex_status(None)
            )
            + r" \\"
        ),
        r"\bottomrule",
        r"\end{tabular}",
        r"\end{center}",
    ]


def _numerical_summary(
    output_dir: Path,
    summary: dict[str, Any],
) -> list[str]:
    """Build numerical running, matrix, and observable tables."""

    lines = [r"\section{Numerical running and observables}"]

    if summary.get("NumericalRGEStatus") != "Success":
        if summary.get("NumericalRGEStatus") == "Failed":
            error = summary.get("NumericalRGEError", "Unknown numerical error")
            lines.append(
                r"\textit{Numerical running failed: }"
                + _latex_escape(error)
            )
        else:
            lines.append(
                r"\textit{Numerical running was not requested. Supply "
                r"\texttt{--numerical CONFIG.json} to include this section.}"
            )
        return lines

    lines.extend(
        [
            r"\begin{center}",
            r"\begin{tabular}{l l}",
            r"\toprule",
            (
                r"Matching scale & $"
                + _scientific_latex(summary["NumericalMatchingScaleGeV"])
                + r"\ \mathrm{GeV}$ \\"
            ),
            (
                r"Low scale & $"
                + _scientific_latex(summary["NumericalLowScaleGeV"])
                + r"\ \mathrm{GeV}$ \\"
            ),
            (
                r"ODE evaluations & "
                + _latex_escape(summary.get("NumericalSolverEvaluations", "N/A"))
                + r" \\"
            ),
            (
                r"Mass ordering & "
                + _latex_escape(summary.get("Ordering", "N/A"))
                + r" \\"
            ),
            r"\bottomrule",
            r"\end{tabular}",
            r"\end{center}",
        ]
    )

    initial_c5 = _read_complex_matrix(
        output_dir,
        summary.get("C5FlavorMatrixNumericFile"),
    )
    low_c5 = _read_complex_matrix(
        output_dir,
        summary.get("C5FlavorMatrixLowScaleFile"),
    )
    low_mass = _read_complex_matrix(
        output_dir,
        summary.get("NeutrinoMassMatrixLowScaleFile"),
    )

    lines.extend(
        _matrix_block(
            "Matched coefficient",
            initial_c5,
            symbol=r"C_5(M)",
            units=r"\mathrm{GeV}^{-1}",
        )
    )
    lines.extend(
        _matrix_block(
            "Low-scale coefficient",
            low_c5,
            symbol=r"C_5(\mu)",
            units=r"\mathrm{GeV}^{-1}",
        )
    )
    lines.extend(
        _matrix_block(
            "Low-scale neutrino-mass matrix",
            low_mass,
            symbol=r"m_\nu(\mu)",
            scale=1.0e9,
            units=r"\mathrm{eV}",
        )
    )

    if summary.get("NeutrinoObservableStatus") != "Success":
        lines.append(r"\textit{Neutrino observables were not available.}")
        return lines

    masses = summary.get("MassesEV", [])
    pmns = summary.get("PMNSAbs", [])

    if len(masses) == 3:
        lines.extend(
            [
                r"\subsection*{Mass spectrum}",
                r"\begin{center}",
                r"\begin{tabular}{c c}",
                r"\toprule",
                r"Quantity & Value \\",
                r"\midrule",
                rf"$m_1$ & ${_scientific_latex(masses[0])}\ \mathrm{{eV}}$ \\",
                rf"$m_2$ & ${_scientific_latex(masses[1])}\ \mathrm{{eV}}$ \\",
                rf"$m_3$ & ${_scientific_latex(masses[2])}\ \mathrm{{eV}}$ \\",
                (
                    r"$\Delta m_{21}^2$ & $"
                    + _scientific_latex(summary.get("DeltaM21SqEV2", 0.0))
                    + r"\ \mathrm{eV}^2$ \\"
                ),
                (
                    r"$\Delta m_{31}^2$ & $"
                    + _scientific_latex(summary.get("DeltaM31SqEV2", 0.0))
                    + r"\ \mathrm{eV}^2$ \\"
                ),
                (
                    r"Takagi residual & $"
                    + _scientific_latex(summary.get("TakagiResidual", 0.0))
                    + r"$ \\"
                ),
                r"\bottomrule",
                r"\end{tabular}",
                r"\end{center}",
            ]
        )

    if (
        isinstance(pmns, list)
        and len(pmns) == 3
        and all(isinstance(row, list) and len(row) == 3 for row in pmns)
    ):
        rows = r" \\ ".join(
            " & ".join(f"{float(value):.5f}" for value in row)
            for row in pmns
        )
        lines.extend(
            [
                r"\subsection*{Absolute PMNS matrix}",
                r"\begin{equation*}",
                rf"|U_{{\rm PMNS}}|=\begin{{pmatrix}}{rows}\end{{pmatrix}}.",
                r"\end{equation*}",
            ]
        )

    return lines


def _output_table(summary: dict[str, Any]) -> list[str]:
    """List machine-readable outputs without exposing empty debug entries."""

    output_keys = (
        ("Matched coefficient", "WeinbergCoefficientFile"),
        ("Rendered matched coefficient", "WeinbergCoefficientPDFFile"),
        ("One-generation beta function", "C5BetaFile"),
        ("Full flavor beta matrix", "C5FlavorBetaMatrixFile"),
        ("Symbolic neutrino-mass matrix", "NeutrinoMassMatrixFile"),
        ("Numerical observables", "NeutrinoObservablesFile"),
        ("Debug beta ratio", "C5BetaOverC5File"),
        ("Debug flavor matrix", "C5FlavorMatrixFile"),
        ("Debug loop kernel", "C5LoopKernelFile"),
    )

    rows = [
        (label, summary.get(key))
        for label, key in output_keys
        if summary.get(key)
    ]

    if not rows:
        return [r"\textit{No machine-readable outputs were recorded.}"]

    lines = [
        r"\begin{center}",
        r"\begin{tabular}{l l}",
        r"\toprule",
        r"Result & File \\",
        r"\midrule",
    ]

    for label, path in rows:
        lines.append(
            _latex_escape(label)
            + r" & \texttt{"
            + _latex_escape(path)
            + r"} \\"
        )

    lines.extend(
        [
            r"\bottomrule",
            r"\end{tabular}",
            r"\end{center}",
        ]
    )
    return lines


def write_rge_report(
    *,
    output_dir: Path,
    model_heading_latex: str,
    summary: dict[str, Any],
    coefficient_latex: str = "",
) -> Path:
    """Write one standalone LaTeX report for a model's complete RGE result."""

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    report_path = output_dir / "rge_report.tex"

    lines: list[str] = [
        r"\documentclass[10pt]{article}",
        r"\usepackage[margin=1.8cm]{geometry}",
        r"\usepackage{amsmath,amssymb,adjustbox,booktabs,array,xcolor}",
        r"\usepackage[T1]{fontenc}",
        r"\allowdisplaybreaks[4]",
        r"\setlength{\emergencystretch}{3em}",
        r"\begin{document}",
        r"\section*{T3 Weinberg-coefficient RGE report}",
        rf"\noindent ${model_heading_latex}$",
        r"\par\medskip",
        *_status_table(summary),
        r"\section{Conventions and scope}",
        r"Below the common heavy-particle threshold, the active theory is SMEFT.",
        r"\begin{align*}",
        r"t&=\ln\mu,\\",
        r"V(H)&=\frac{\lambda_H}{2}(H^\dagger H)^2.",
        r"\end{align*}",
        (
            r"The flavor construction assumes three lepton generations, a "
            r"diagonal heavy-fermion mass basis, and common scalar masses."
        ),
        r"\section{Matched Weinberg coefficient}",
    ]

    lines.extend(_coefficient_block(coefficient_latex))

    lines.extend(
        [
            r"\section{One-generation SMEFT check}",
            (
                r"The tensor calculation is required to reproduce the compact "
                r"one-generation benchmark"
            ),
            r"\begin{equation}",
            r"16\pi^2\beta_{C_5}=C_5\left[",
            r"-3g_2^2+2\lambda_H+6|y_u|^2+6|y_d|^2-|y_e|^2",
            r"\right].",
            r"\end{equation}",
            (
                r"The RGE stage stops with an error if the calculated expression "
                r"does not reduce to this benchmark."
            ),
            r"\section{Three-generation flavor result}",
            r"The matched one-generation result is factorized as",
            r"\begin{equation*}",
            r"C_5^{(1)}=F_{\rm loop}\,\bar y_1\bar y_2,",
            r"\end{equation*}",
            r"and lifted to the symmetric flavor matrix",
            r"\begin{equation}",
            r"(C_5)_{pq}=\frac12\sum_r F_r\left[",
            r"y_{1,pr}^*y_{2,qr}^*+y_{2,pr}^*y_{1,qr}^*",
            r"\right].",
            r"\end{equation}",
            r"Writing this matrix as $K\equiv C_5$, the implemented RGE is",
            r"\begin{align}",
            r"16\pi^2\frac{dK}{d\ln\mu}",
            r"&=(2\lambda_H-3g_2^2+2T)K",
            r"-\frac32\left[Y_eY_e^\dagger K+K(Y_eY_e^\dagger)^T\right],\\",
            r"T&=\operatorname{Tr}\!\left(",
            r"Y_eY_e^\dagger+3Y_uY_u^\dagger+3Y_dY_d^\dagger",
            r"\right).",
            r"\end{align}",
            (
                r"The fully expanded symbolic matrix is retained as a data file "
                r"rather than printed here, because it is much less readable than "
                r"the matrix equation above."
            ),
            r"\section{Neutrino-mass conversion}",
            r"The current symbolic and numerical implementation uses",
            r"\begin{equation}",
            r"m_\nu=-v^2C_5.",
            r"\end{equation}",
            *_numerical_summary(output_dir, summary),
            r"\section{Saved data}",
            *_output_table(summary),
            r"\end{document}",
            "",
        ]
    )

    report_path.write_text(
        "\n".join(lines),
        encoding="utf-8",
        newline="\n",
    )
    return report_path


def write_and_compile_rge_report(record: RunRecord) -> Path:
    """Create and compile the concise per-model RGE report."""

    report_tex = write_rge_report(
        output_dir=report_output_dir_for(record),
        model_heading_latex=latex_model_heading(record),
        summary=record.summary,
        coefficient_latex=record.summary.get(
            "WeinbergCoefficientLaTeX",
            "",
        ),
    )

    record.summary["RGEReportTeXFile"] = report_tex.name
    compile_latex_document(report_tex)

    report_pdf = report_tex.with_suffix(".pdf")
    if report_pdf.exists():
        record.summary["RGEReportStatus"] = "Success"
        record.summary["RGEReportPDFFile"] = report_pdf.name
        print(f"  {record.name}: RGE report -> {report_pdf}")
    else:
        record.summary["RGEReportStatus"] = "TeXOnly"
        record.summary["RGEReportPDFFile"] = ""
        print(f"  {record.name}: RGE report source -> {report_tex}")

    return report_tex
