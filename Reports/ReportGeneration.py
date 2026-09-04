from __future__ import annotations

import re
import shutil
import subprocess
import sys
from fractions import Fraction
from pathlib import Path

from common.Paths import EFT_ORDER, LOOP_ORDER, REPORT_OUTPUT_DIR
from common.Records import RunRecord


def report_output_dir_for(record: RunRecord) -> Path:
    # Human-readable per-model reports are kept away from raw calculation data.
    output_dir = REPORT_OUTPUT_DIR / record.output_dir.name
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


# ---------------------------------------------------------------------------
# LaTeX/report helpers
# ---------------------------------------------------------------------------

def summary_fraction(summary: dict, key: str) -> Fraction | None:
    """Read an exact number such as -1/2 from the Wolfram summary."""

    raw_value = summary.get(key)

    if raw_value in (None, ""):
        return None

    try:
        return Fraction(str(raw_value))
    except (ValueError, ZeroDivisionError):
        return None


def latex_fraction(value: Fraction | None) -> str:
    """Convert an exact Fraction to LaTeX."""

    if value is None:
        return r"\text{N/A}"

    if value.denominator == 1:
        return str(value.numerator)

    sign = "-" if value < 0 else ""
    magnitude = abs(value)

    return (
        rf"{sign}\frac{{{magnitude.numerator}}}"
        rf"{{{magnitude.denominator}}}"
    )


def record_quantum_numbers(
    record: RunRecord,
) -> tuple[int, Fraction | None, int, Fraction | None, int, Fraction | None]:
    """Return dimensions and actual hypercharges reported by Wolfram."""

    summary = record.summary

    return (
        int(summary.get("Scalar1SU2", record.d_s1)),
        summary_fraction(summary, "Scalar1Hypercharge"),
        int(summary.get("Scalar2SU2", record.d_s2)),
        summary_fraction(summary, "Scalar2Hypercharge"),
        int(summary.get("FermionSU2", record.d_f)),
        summary_fraction(summary, "FermionHypercharge"),
    )


def latex_model_heading(record: RunRecord) -> str:
    """Create the representation heading used in generated reports."""

    d_s1, y_s1, d_s2, y_s2, d_f, y_f = record_quantum_numbers(record)

    return (
        rf"{record.name},\ \alpha={record.alpha}:\quad "
        rf"S_1=({d_s1},{latex_fraction(y_s1)}),\quad "
        rf"S_2=({d_s2},{latex_fraction(y_s2)}),\quad "
        rf"F=({d_f},{latex_fraction(y_f)})"
    )


def split_latex_equation(latex: str, target_length: int = 92) -> list[str]:
    """Split a LaTeX sum at top-level + and - operators."""

    if not latex:
        return [r"\text{Not available}"]

    terms = split_latex_terms(latex)

    if not terms:
        return [latex.strip()]

    lines: list[str] = []
    current_line = ""

    for term in terms:
        candidate = f"{current_line} {term}".strip()

        if current_line and len(candidate) > target_length:
            lines.append(current_line)
            current_line = term
        else:
            current_line = candidate

    if current_line:
        lines.append(current_line)

    return lines


def latex_aligned_block(latex: str) -> list[str]:
    """Create page-breakable LaTeX displays for a long expression."""

    output: list[str] = []

    for equation_line in split_latex_equation(latex):
        output.extend(
            [
                (
                    r"\noindent\adjustbox{max width=\linewidth}"
                    r"{$\displaystyle "
                    + equation_line
                    + r"$}\par"
                ),
                r"\smallskip",
            ]
        )

    return output


def split_latex_terms(latex: str) -> list[str]:
    """Split a LaTeX expression into signed top-level additive terms."""

    if not latex or not latex.strip():
        return []

    text = latex.strip()

    terms: list[str] = []
    current: list[str] = []

    brace_depth = 0
    paren_depth = 0
    bracket_depth = 0
    escaped = False

    for index, char in enumerate(text):
        if escaped:
            current.append(char)
            escaped = False
            continue

        if char == "\\":
            current.append(char)
            escaped = True
            continue

        if char == "{":
            brace_depth += 1
        elif char == "}":
            brace_depth = max(0, brace_depth - 1)
        elif char == "(":
            paren_depth += 1
        elif char == ")":
            paren_depth = max(0, paren_depth - 1)
        elif char == "[":
            bracket_depth += 1
        elif char == "]":
            bracket_depth = max(0, bracket_depth - 1)

        top_level = (
            brace_depth == 0
            and paren_depth == 0
            and bracket_depth == 0
        )

        previous = text[index - 1] if index else ""

        if (
            char in "+-"
            and top_level
            and current
            and previous not in "^_eE"
        ):
            term = "".join(current).strip()

            if term:
                terms.append(term)

            current = [char]
        else:
            current.append(char)

    term = "".join(current).strip()

    if term:
        terms.append(term)

    return terms


# The Lagrangian builder currently prints the generic T3 heavy fields as
# NewScalar1, NewScalar2 and N. These names are translated here into the
# physical labels S1, S2 and F used by the pipeline.
FIELD_PATTERNS: tuple[tuple[str, tuple[str, ...]], ...] = (
    (r"S_1", (r"\\text\{NewScalar1\}",)),
    (r"S_2", (r"\\text\{NewScalar2\}",)),
    (r"F", (r"N(?:_|\^)",)),
    (r"H", (r"H(?:_|\^)",)),
    (r"\ell", (r"\\ell", r"\\mathcal\{l\}", r"\\mathscr\{l\}")),
    (r"e", (r"(?<![A-Za-z\\])e(?:_|\^)",)),
    (r"q", (r"(?<![A-Za-z\\])q(?:_|\^)",)),
    (r"u", (r"(?<![A-Za-z\\])u(?:_|\^)",)),
    (r"d", (r"(?<![A-Za-z\\])d(?:_|\^)",)),
    (r"B_{\mu\nu}", (r"B(?:_|\^)",)),
    (r"W_{\mu\nu}", (r"W(?:_|\^)",)),
    (r"G_{\mu\nu}", (r"G(?:X)?(?:_|\^)",)),
)


def latex_field_signature(term: str) -> tuple[str, ...]:
    """Return the physical fields occurring in one LaTeX term."""

    searchable_term = term

    # Remove Lorentz-index/scale commands before looking for fields.
    ignored_commands = (
        r"\mubar",
        r"\overline\mu",
        r"\bar\mu",
        r"\mu",
        r"\nu",
        r"\rho",
        r"\sigma",
        r"\alpha",
        r"\beta",
        r"\gamma",
        r"\delta",
        r"\epsilon",
        r"\varepsilon",
        r"\hbar",
    )

    for command in ignored_commands:
        searchable_term = re.sub(
            re.escape(command),
            " ",
            searchable_term,
        )

    fields: list[str] = []

    for label, patterns in FIELD_PATTERNS:
        if any(
            re.search(pattern, searchable_term)
            for pattern in patterns
        ):
            fields.append(label)

    # Keep covariant derivatives as a useful distinction for kinetic terms.
    if re.search(
        r"(?:\\mathrm\{D\}|\\mathcal\{D\}|D)(?:_|\^)",
        term,
    ):
        fields.append(r"D")

    return tuple(fields) or (r"\text{constant}",)


def grouped_lagrangian_terms(
    latex: str,
) -> dict[tuple[str, ...], list[str]]:
    """Group additive Lagrangian terms by field content."""

    grouped: dict[tuple[str, ...], list[str]] = {}

    for term in split_latex_terms(latex):
        signature = latex_field_signature(term)

        grouped.setdefault(signature, []).append(term)

    return grouped


def field_signature_label(signature: tuple[str, ...]) -> str:
    """Render a field combination as a LaTeX table heading."""

    return "$" + r"\,".join(signature) + "$"


def matrix_cell(terms: list[str], empty_value: str = "") -> str:
    """Render all terms belonging to one field combination."""

    if not terms:
        return empty_value

    joined = r" \\ ".join(terms)

    return (
        r"\scriptsize\adjustbox{max width=\linewidth}"
        r"{$\begin{gathered} "
        + joined
        + r"\end{gathered}$}"
    )


# ---------------------------------------------------------------------------
# Primary Lagrangian report
# ---------------------------------------------------------------------------

def latex_escape_text(value: object) -> str:
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


def latex_status(value: object) -> str:
    """Render one pipeline status as compact LaTeX text."""

    if value is True or value == "Success":
        return r"\textbf{Success}"

    if value is False:
        return r"\textbf{Failed}"

    return latex_escape_text(value if value not in (None, "") else "Unknown")


def append_report_expression(
    lines: list[str],
    title: str,
    latex: str,
    available: bool,
) -> None:
    """Append one readable expression or an explicit unavailable notice."""

    lines.append(rf"\subsubsection*{{{title}}}")

    if available and latex.strip():
        lines.extend(latex_aligned_block(latex))
    else:
        lines.append(r"\textit{Not available for this model.}")


def write_lagrangian_report(records: list[RunRecord]) -> Path:
    """Create the concise human-facing Lagrangian and Weinberg report."""

    output_path = REPORT_OUTPUT_DIR / "lagrangian_report.tex"

    lines: list[str] = [
        r"\documentclass[10pt]{article}",
        r"\usepackage[margin=1.7cm]{geometry}",
        r"\usepackage{amsmath,amssymb,adjustbox,booktabs,array}",
        r"\usepackage[T1]{fontenc}",
        r"\allowdisplaybreaks[4]",
        r"\setlength{\emergencystretch}{3em}",
        r"\begin{document}",
        r"\section*{T3 Lagrangian and Weinberg-operator report}",
        rf"EFT order: ${EFT_ORDER}$; loop order: ${LOOP_ORDER}$.",
    ]

    if not records:
        lines.append(r"No models were run.")

    for model_number, record in enumerate(records, start=1):
        summary = record.summary

        if model_number > 1:
            lines.append(r"\clearpage")

        lines.extend(
            [
                rf"\subsection*{{$" + latex_model_heading(record) + r"$}",
                r"\begin{center}",
                r"\begin{tabular}{>{\bfseries}l l >{\bfseries}l l}",
                r"\toprule",
                (
                    r"Build & "
                    + latex_status(summary.get("BuildStatus"))
                    + r" & Matching & "
                    + latex_status(summary.get("MatchingStatus"))
                    + r" \\"
                ),
                (
                    r"T3 ingredients & "
                    + latex_status(summary.get("T3IngredientsPresent"))
                    + r" & Weinberg extraction & "
                    + latex_status(summary.get("WeinbergExtractionStatus"))
                    + r" \\"
                ),
                r"\bottomrule",
                r"\end{tabular}",
                r"\end{center}",
            ]
        )

        accepted = summary.get("AcceptedInteractions", [])
        lines.append(r"\subsubsection*{Accepted interactions}")

        if accepted:
            lines.append(
                ", ".join(
                    rf"\texttt{{{latex_escape_text(name)}}}"
                    for name in accepted
                )
                + "."
            )
        else:
            lines.append(r"\textit{No accepted interactions were reported.}")

        append_report_expression(
            lines,
            r"BSM free Lagrangian $\mathcal{L}_{\mathrm{free}}$",
            summary.get("FreeLagrangianLaTeX", ""),
            summary.get("FreeLagrangianConversionSuccess") is True,
        )
        append_report_expression(
            lines,
            r"Interaction Lagrangian $\mathcal{L}_{\mathrm{int}}$",
            summary.get("InteractionLagrangianLaTeX", ""),
            summary.get("InteractionLagrangianConversionSuccess") is True,
        )
        append_report_expression(
            lines,
            (
                r"BSM-induced EFT contribution "
                r"$\Delta\mathcal{L}_{\mathrm{EFT}}$"
            ),
            summary.get("BSMEFTLagrangianLaTeX", ""),
            summary.get("BSMEFTConversionSuccess") is True,
        )
        append_report_expression(
            lines,
            r"Weinberg sector",
            summary.get("WeinbergSectorLaTeX", ""),
            summary.get("WeinbergSectorConversionSuccess") is True,
        )
        append_report_expression(
            lines,
            r"Weinberg coefficient $C_5$",
            summary.get("WeinbergCoefficientLaTeX", ""),
            summary.get("WeinbergCoefficientConversionSuccess") is True,
        )

    lines.extend(
        [
            r"\end{document}",
            "",
        ]
    )

    output_path.write_text(
        "\n".join(lines),
        encoding="utf-8",
    )

    print(f"\nLagrangian report:\n{output_path}")

    return output_path


# ---------------------------------------------------------------------------
# Full Lagrangian reports
# ---------------------------------------------------------------------------

def write_latex_lagrangian_table(records: list[RunRecord]) -> Path:
    """Create the complete UV and matched EFT Lagrangian report."""

    output_path = REPORT_OUTPUT_DIR / "eft_lagrangians.tex"

    lines: list[str] = [
        r"\documentclass[10pt]{article}",
        r"\usepackage[margin=1.5cm]{geometry}",
        r"\usepackage{amsmath,amssymb,adjustbox,pdflscape}",
        r"\usepackage[T1]{fontenc}",
        r"\allowdisplaybreaks[4]",
        r"\setlength{\emergencystretch}{3em}",
        r"\begin{document}",
        r"\section*{T3 UV and EFT Lagrangians}",
        rf"EFT order: ${EFT_ORDER}$; loop order: ${LOOP_ORDER}$.",
        r"\begin{landscape}",
        r"\footnotesize",
    ]

    included_models = 0

    for record in records:
        summary = record.summary

        if (
            summary.get("BuildStatus") != "Success"
            or summary.get("MatchingStatus") != "Success"
        ):
            continue

        if summary.get("UVConversionSuccess") is not True:
            continue

        if summary.get("EFTConversionSuccess") is not True:
            continue

        uv_latex = summary.get("UVLagrangianLaTeX", "").strip()
        eft_latex = summary.get("EFTLagrangianLaTeX", "").strip()

        if not uv_latex or not eft_latex:
            continue

        included_models += 1

        lines.extend(
            [
                rf"\subsection*{{$"
                + latex_model_heading(record)
                + r"$}",
                r"\textbf{Full UV Lagrangian}",
                *latex_aligned_block(uv_latex),
                r"\textbf{Matched EFT Lagrangian}",
                *latex_aligned_block(eft_latex),
                r"\medskip\hrule\medskip",
            ]
        )

    if included_models == 0:
        lines.append(
            r"No fully converted UV/EFT results are available."
        )

    lines.extend(
        [
            r"\end{landscape}",
            r"\end{document}",
            "",
        ]
    )

    output_path.write_text(
        "\n".join(lines),
        encoding="utf-8",
    )

    print(f"\nFull UV/EFT Lagrangian report:\n{output_path}")

    return output_path


def write_bsm_lagrangian_table(records: list[RunRecord]) -> Path:
    """Create the BSM UV and induced EFT Lagrangian report."""

    output_path = REPORT_OUTPUT_DIR / "bsm_lagrangians.tex"

    lines: list[str] = [
        r"\documentclass[10pt]{article}",
        r"\usepackage[margin=1.5cm]{geometry}",
        r"\usepackage{amsmath,amssymb,adjustbox,pdflscape}",
        r"\usepackage[T1]{fontenc}",
        r"\allowdisplaybreaks[4]",
        r"\setlength{\emergencystretch}{3em}",
        r"\begin{document}",
        r"\section*{T3 BSM UV and matched EFT contributions}",
        rf"EFT order: ${EFT_ORDER}$; loop order: ${LOOP_ORDER}$.",
        (
            r"The EFT expressions show "
            r"$\Delta\mathcal{L}_{\mathrm{EFT}}="
            r"\mathcal{L}_{\mathrm{matched}}^{\mathrm{SM+BSM}}"
            r"-\mathcal{L}_{\mathrm{matched}}^{\mathrm{SM}}$."
        ),
        r"\begin{landscape}",
        r"\footnotesize",
    ]

    included_models = 0

    for record in records:
        summary = record.summary

        if (
            summary.get("BuildStatus") != "Success"
            or summary.get("MatchingStatus") != "Success"
        ):
            continue

        uv_ok = summary.get("BSMUVConversionSuccess") is True
        eft_ok = summary.get("BSMEFTConversionSuccess") is True

        uv_latex = summary.get(
            "BSMUVLagrangianLaTeX",
            "",
        ).strip()

        eft_latex = summary.get(
            "BSMEFTLagrangianLaTeX",
            "",
        ).strip()

        if not uv_ok and not eft_ok:
            continue

        included_models += 1

        lines.extend(
            [
                rf"\subsection*{{$"
                + latex_model_heading(record)
                + r"$}",
                r"\textbf{BSM UV Lagrangian}",
            ]
        )

        if uv_ok and uv_latex:
            lines.extend(latex_aligned_block(uv_latex))
        else:
            lines.append(r"\textit{BSM UV Lagrangian unavailable.}")

        lines.append(
            r"\textbf{$\Delta\mathcal{L}_{\mathrm{EFT}}$}"
        )

        if eft_ok and eft_latex:
            lines.extend(latex_aligned_block(eft_latex))
        else:
            lines.append(
                r"\textit{Matched BSM EFT difference unavailable.}"
            )

        lines.append(r"\medskip\hrule\medskip")

    if included_models == 0:
        lines.append(
            r"No fully converted BSM results are available."
        )

    lines.extend(
        [
            r"\end{landscape}",
            r"\end{document}",
            "",
        ]
    )

    output_path.write_text(
        "\n".join(lines),
        encoding="utf-8",
    )

    print(f"\nBSM-only Lagrangian report:\n{output_path}")

    return output_path


# ---------------------------------------------------------------------------
# Lagrangian term tables
# ---------------------------------------------------------------------------

def write_bsm_field_table(
    records: list[RunRecord],
    *,
    summary_key: str,
    conversion_key: str,
    output_stem: str,
    title: str,
    description: str,
    empty_value: str = "",
) -> Path:
    """Create a table grouping terms according to their field content."""

    output_path = REPORT_OUTPUT_DIR / f"{output_stem}.tex"

    rows = []
    signatures: set[tuple[str, ...]] = set()

    for record in records:
        summary = record.summary

        if (
            summary.get("BuildStatus") != "Success"
            or summary.get("MatchingStatus") != "Success"
            or summary.get(conversion_key) is not True
        ):
            continue

        grouped = grouped_lagrangian_terms(
            summary.get(summary_key, "")
        )

        if not grouped:
            continue

        quantum_numbers = record_quantum_numbers(record)

        rows.append(
            (
                record,
                quantum_numbers,
                grouped,
            )
        )

        signatures.update(grouped)

    ordered_signatures = sorted(
        signatures,
        key=lambda signature: (
            len(signature),
            signature,
        ),
    )

    lines: list[str] = [
        r"\documentclass[8pt]{article}",
        r"\usepackage[margin=0.8cm]{geometry}",
        (
            r"\usepackage{amsmath,amssymb,adjustbox,"
            r"pdflscape,longtable,array,booktabs}"
        ),
        r"\usepackage[T1]{fontenc}",
        r"\setlength{\tabcolsep}{2pt}",
        r"\renewcommand{\arraystretch}{1.25}",
        r"\begin{document}",
        r"\begin{landscape}",
        rf"\section*{{{title}}}",
        rf"EFT order: ${EFT_ORDER}$; loop order: ${LOOP_ORDER}$. ",
        description,
        r"\tiny",
    ]

    if not rows or not ordered_signatures:
        lines.append(
            r"No fully converted results are available for this sector."
        )

    else:
        # Limit the number of field-content columns on each page.
        chunk_size = 4

        chunks = [
            ordered_signatures[index:index + chunk_size]
            for index in range(
                0,
                len(ordered_signatures),
                chunk_size,
            )
        ]

        for page_number, signature_chunk in enumerate(
            chunks,
            start=1,
        ):
            if len(chunks) > 1:
                lines.append(
                    rf"\subsection*{{Field combinations "
                    rf"{page_number} of {len(chunks)}}}"
                )

            widths = " ".join(
                r">{\raggedright\arraybackslash}p{0.18\linewidth}"
                for _ in signature_chunk
            )

            column_spec = (
                r"@{}cccccccc "
                + widths
                + r"@{}"
            )

            header = [
                r"Model",
                r"$\alpha$",
                r"$d_{S_1}$",
                r"$Y_{S_1}$",
                r"$d_{S_2}$",
                r"$Y_{S_2}$",
                r"$d_F$",
                r"$Y_F$",
                *[
                    field_signature_label(signature)
                    for signature in signature_chunk
                ],
            ]

            lines.extend(
                [
                    rf"\begin{{longtable}}{{{column_spec}}}",
                    r"\toprule",
                    " & ".join(header) + r" \\",
                    r"\midrule",
                    r"\endfirsthead",
                    r"\toprule",
                    " & ".join(header) + r" \\",
                    r"\midrule",
                    r"\endhead",
                ]
            )

            for record, quantum_numbers, grouped in rows:
                (
                    d_s1,
                    y_s1,
                    d_s2,
                    y_s2,
                    d_f,
                    y_f,
                ) = quantum_numbers

                cells = [
                    record.name,
                    rf"${record.alpha}$",
                    rf"${d_s1}$",
                    rf"${latex_fraction(y_s1)}$",
                    rf"${d_s2}$",
                    rf"${latex_fraction(y_s2)}$",
                    rf"${d_f}$",
                    rf"${latex_fraction(y_f)}$",
                    *[
                        matrix_cell(
                            grouped.get(signature, []),
                            empty_value,
                        )
                        for signature in signature_chunk
                    ],
                ]

                lines.append(
                    " & ".join(cells)
                    + r" \\"
                )

                lines.append(r"\midrule")

            lines.extend(
                [
                    r"\bottomrule",
                    r"\end{longtable}",
                ]
            )

            if page_number != len(chunks):
                lines.append(r"\clearpage")

    lines.extend(
        [
            r"\end{landscape}",
            r"\end{document}",
            "",
        ]
    )

    output_path.write_text(
        "\n".join(lines),
        encoding="utf-8",
    )

    print(f"\nBSM field-combination table:\n{output_path}")

    return output_path


def write_bsm_uv_field_table(
    records: list[RunRecord],
    empty_value: str = "",
) -> Path:
    """Create the BSM UV table using UV-only field-combination columns."""

    return write_bsm_field_table(
        records,
        summary_key="BSMUVLagrangianLaTeX",
        conversion_key="BSMUVConversionSuccess",
        output_stem="bsm_uv_field_table",
        title="T3 BSM UV terms grouped by field content",
        description=(
            r"Columns are determined only from BSM UV terms. "
            r"Empty cells denote no UV term with that field combination."
        ),
        empty_value=empty_value,
    )


def write_bsm_matched_field_table(
    records: list[RunRecord],
    empty_value: str = "",
) -> Path:
    """Create the matched BSM table using matched-only field columns."""

    return write_bsm_field_table(
        records,
        summary_key="BSMEFTLagrangianLaTeX",
        conversion_key="BSMEFTConversionSuccess",
        output_stem="bsm_matched_field_table",
        title=(
            r"T3 matched BSM EFT terms by field content"
        ),
        description=(
            r"Columns are determined only from matched BSM EFT terms. "
            r"Empty cells denote no matched term with that field combination."
        ),
        empty_value=empty_value,
    )


def write_c5_coefficient_report(record: RunRecord) -> Path:
    """Create a standalone PDF-ready document for one C5 coefficient."""

    output_path = report_output_dir_for(record) / "c5_coefficient.tex"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    coefficient_latex = record.summary.get(
        "WeinbergCoefficientLaTeX",
        "",
    ).strip()

    lines: list[str] = [
        r"\documentclass[11pt]{article}",
        r"\usepackage[margin=2cm]{geometry}",
        r"\usepackage{amsmath,amssymb,adjustbox}",
        r"\usepackage[T1]{fontenc}",
        r"\allowdisplaybreaks[4]",
        r"\setlength{\emergencystretch}{3em}",
        r"\begin{document}",
        r"\section*{Weinberg-operator coefficient}",
        rf"\noindent ${latex_model_heading(record)}$",
        r"\par",
        r"\bigskip",
    ]

    if (
        record.summary.get("WeinbergCoefficientConversionSuccess") is True
        and coefficient_latex
    ):
        lines.extend(
            latex_aligned_block(
                r"C_5 = " + coefficient_latex
            )
        )
    else:
        lines.append(
            r"\textit{The Weinberg coefficient is not available.}"
        )

    lines.extend(
        [
            r"\end{document}",
            "",
        ]
    )

    output_path.write_text(
        "\n".join(lines),
        encoding="utf-8",
    )

    print(f"\nC5 coefficient report:\n{output_path}")

    return output_path


def compile_latex_document(tex_path: Path) -> None:
    """Compile a generated LaTeX report using latexmk or pdflatex."""

    latexmk = shutil.which("latexmk")
    pdflatex = shutil.which("pdflatex")

    if latexmk is None and pdflatex is None:
        print(
            "latexmk and pdflatex were not found; the .tex file was generated "
            "but not compiled."
        )
        return

    outputs: list[str] = []
    compilation_succeeded = False

    tex_path.with_suffix(".pdf").unlink(missing_ok=True)

    if latexmk is not None:
        latexmk_result = subprocess.run(
            [
                latexmk,
                "-pdf",
                "-interaction=nonstopmode",
                "-halt-on-error",
                tex_path.name,
            ],
            cwd=tex_path.parent,
            capture_output=True,
            text=True,
            check=False,
        )
        outputs.extend(
            [
                "--- latexmk output ---",
                latexmk_result.stdout,
                latexmk_result.stderr,
            ]
        )
        compilation_succeeded = latexmk_result.returncode == 0

        if not compilation_succeeded and pdflatex is not None:
            print(
                "latexmk failed; falling back to pdflatex."
            )

    if not compilation_succeeded and pdflatex is not None:
        compilation_succeeded = True

        for pass_number in (1, 2):
            pdflatex_result = subprocess.run(
                [
                    pdflatex,
                    "-interaction=nonstopmode",
                    "-halt-on-error",
                    tex_path.name,
                ],
                cwd=tex_path.parent,
                capture_output=True,
                text=True,
                check=False,
            )
            outputs.extend(
                [
                    f"--- pdflatex pass {pass_number} output ---",
                    pdflatex_result.stdout,
                    pdflatex_result.stderr,
                ]
            )

            if pdflatex_result.returncode != 0:
                compilation_succeeded = False
                break

    if not compilation_succeeded:
        log_path = tex_path.with_suffix(".compile.log")

        log_path.write_text(
            "\n".join(outputs),
            encoding="utf-8",
        )

        print(
            f"LaTeX compilation failed. See:\n{log_path}",
            file=sys.stderr,
        )

        return

    tex_path.with_suffix(".compile.log").unlink(
        missing_ok=True
    )

    # Remove auxiliary LaTeX files.
    for suffix in (
        ".aux",
        ".log",
        ".out",
        ".fls",
        ".fdb_latexmk",
        ".synctex.gz",
    ):
        tex_path.with_suffix(suffix).unlink(
            missing_ok=True
        )

    print(
        f"Compiled PDF:\n"
        f"{tex_path.with_suffix('.pdf')}"
    )


def write_reports(
    records: list[RunRecord],
    debug_reports: bool = False,
) -> None:
    """Generate all Lagrangian reports and term tables."""

    report_tex = write_lagrangian_report(records)
    compile_latex_document(report_tex)

    uv_table_tex = write_bsm_uv_field_table(records)
    compile_latex_document(uv_table_tex)

    matched_table_tex = write_bsm_matched_field_table(records)
    compile_latex_document(matched_table_tex)

    for record in records:
        summary = record.summary

        if (
            summary.get("WeinbergCoefficientConversionSuccess") is not True
            or not summary.get("WeinbergCoefficientLaTeX", "").strip()
        ):
            continue

        coefficient_tex = write_c5_coefficient_report(record)
        compile_latex_document(coefficient_tex)

    if debug_reports:
        full_tex = write_latex_lagrangian_table(records)
        compile_latex_document(full_tex)

        bsm_tex = write_bsm_lagrangian_table(records)
        compile_latex_document(bsm_tex)


