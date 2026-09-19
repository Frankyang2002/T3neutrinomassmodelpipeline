from __future__ import annotations

import os

import re
import shutil
import subprocess
import sys
from fractions import Fraction
from pathlib import Path
from typing import Iterable

from common.Paths import EFT_ORDER, LOOP_ORDER, REPORT_OUTPUT_DIR
from common.Records import RunRecord



# ---------------------------------------------------------------------------
# Shared report notation
# Consolidated from the former Reports/Notation.py module.
# Presentation-only: no matching normalization or RGE physics is defined here.
# ---------------------------------------------------------------------------

PAPER_SYMBOL_LATEX: dict[str, str] = {
    "gY": r"g_1",
    "g2": r"g_2",
    "g3": r"g_3",
    "yu": r"Y_u",
    "yd": r"Y_d",
    "ye": r"Y_e",
    "y1": r"y_1",
    "y2": r"y_2",
    "h": r"h",
    "MF": r"M_F",
    "mSSq": r"m_S^2",
    "mS1Sq": r"m_1^2",
    "mS2Sq": r"m_2^2",
    "lambdaH": r"\lambda_1",
    "lambdaS": r"\lambda_2",
    "lambda3": r"\lambda_3",
    "lambda4": r"\lambda_4",
    "lambda5": r"\lambda_5",
    "lambdaS1": r"\lambda_{S_1}^{(1)}",
    "lambdaS2": r"\lambda_{S_2}^{(1)}",
    "lambdaH1": r"\lambda_{HS_1}^{(1)}",
    "lambdaH2": r"\lambda_{HS_2}^{(1)}",
    "lambdaH1Adj": r"\lambda_{HS_1}^{(A)}",
    "lambdaH2Adj": r"\lambda_{HS_2}^{(A)}",
    "lambdaT3": r"\lambda_5",
    "lambda12": r"\lambda_{12}^{(1)}",
    "lambda12Adj": r"\lambda_{12}^{(A)}",
    "lambda12Cross": r"\lambda_{12}^{(\times)}",
    "lambdaS1Adj": r"\lambda_{S_1}^{(A)}",
    "lambdaS2Adj": r"\lambda_{S_2}^{(A)}",
    "lambdaHHdagS1barS1bar": r"\lambda_{7}",
    "lambdaHHdagS2S2": r"\lambda_{8}",
    "lambdaS1bar2S2bar2": r"\lambda_{9}",
    "lambdaS1barS2S2bar2": r"\lambda_{10}",
    "lambdaS1S1bar2S2bar": r"\lambda_{11}",
    "lambdaHHdagS1barS2barCross": r"\lambda_{12}",
}


# RGBeta's raw special-cross symbol historically renders with its explicit
# field-content name before the later paper-notation cleanup. Keep that
# presentation behavior explicit rather than hiding it in RGEComparison.
RGBETA_SYMBOL_OVERRIDES: dict[str, str] = {
    "lambdaHHdagS1barS2barCross": (
        r"\lambda_{H H^\dagger S_1^\dagger S_2^\dagger}^{\rm Cross}"
    ),
}


RENDERED_PAPER_REPLACEMENTS: dict[str, str] = {
    r"\lambda_{H H^\dagger S_1^\dagger S_2^\dagger}^{\rm Cross}": r"\lambda_{12}",
    r"\lambda_{S_1 S_1^{\dagger 2}S_2^\dagger}": r"\lambda_{11}",
    r"\lambda_{S_1^\dagger S_2 S_2^{\dagger 2}}": r"\lambda_{10}",
    r"\lambda_{S_1^{\dagger 2}S_2^{\dagger 2}}": r"\lambda_9",
    r"\lambda_{H^\dagger H^\dagger S_2 S_2}": r"\lambda_8",
    r"\lambda_{H^\dagger H^\dagger S_1^\dagger S_1^\dagger}": r"\lambda_7",
    r"\lambda_{12}^{\rm Cross}": r"\lambda_{12}^{(\times)}",
    r"\lambda_{12}^{\rm Adj}": r"\lambda_{12}^{(A)}",
    r"\lambda_{H2}^{\rm Adj}": r"\lambda_{HS_2}^{(A)}",
    r"\lambda_{H1}^{\rm Adj}": r"\lambda_{HS_1}^{(A)}",
    r"\lambda_{S_2}^{\rm Adj}": r"\lambda_{S_2}^{(A)}",
    r"\lambda_{S_1}^{\rm Adj}": r"\lambda_{S_1}^{(A)}",
    r"\lambda_{T3}": r"\lambda_5",
    r"\lambda_{12}": r"\lambda_{12}^{(1)}",
    r"\lambda_{H2}": r"\lambda_{HS_2}^{(1)}",
    r"\lambda_{H1}": r"\lambda_{HS_1}^{(1)}",
    r"\lambda_{S_2}": r"\lambda_{S_2}^{(1)}",
    r"\lambda_{S_1}": r"\lambda_{S_1}^{(1)}",
    r"\lambda_H": r"\lambda_1",
    r"g_Y": r"g_1",
    r"m_{S_1}^2": r"m_1^2",
    r"m_{S_2}^2": r"m_2^2",
    r"M_{S_1}": r"m_1",
    r"M_{S_2}": r"m_2",
}


def paper_symbol_latex(name: str) -> str:
    """Return the canonical human-facing symbol for one pipeline quantity."""
    return PAPER_SYMBOL_LATEX.get(str(name), str(name))


def rgbeta_symbol_latex(name: str) -> str:
    """Return the report symbol used while normalising raw RGBeta LaTeX."""
    key = str(name)
    return RGBETA_SYMBOL_OVERRIDES.get(key, paper_symbol_latex(key))


BETA_LATEX_OVERRIDES: dict[str, str] = {
    # Preserve exact historical GroupFactorReports spelling.
    "h": r"\beta_h",
    "mSSq": r"\beta_{m_S^{2}}",
    "mS1Sq": r"\beta_{m_1^{2}}",
    "mS2Sq": r"\beta_{m_2^{2}}",
}


def beta_symbol_latex(name: str) -> str:
    """Return the canonical human-facing beta-function label."""
    key = str(name)

    if key in BETA_LATEX_OVERRIDES:
        return BETA_LATEX_OVERRIDES[key]

    if key in PAPER_SYMBOL_LATEX:
        return rf"\beta_{{{PAPER_SYMBOL_LATEX[key]}}}"

    return rf"\beta_{{\mathrm{{{key}}}}}"


def paper_notation_key_lines(*, heading: str = "Notation key") -> list[str]:
    """Return the common symbol/index key used in human-facing reports."""
    return [
        rf"\section*{{{heading}}}",
        r"\begin{longtable}{@{}p{0.23\linewidth}p{0.69\linewidth}@{}}",
        r"\toprule",
        r"symbol & interaction / definition \\",
        r"\midrule",
        r"$g_1$ & $D_\mu\supset i g_1YB_\mu$ \\",
        r"$g_2$ & $D_\mu\supset i g_2T^AW_\mu^A$ \\",
        r"$g_3$ & $D_\mu\supset i g_3t^AG_\mu^A$ \\",
        r"$Y_u$ & $-\bar QY_u\widetilde Hu+\mathrm{h.c.}$ \\",
        r"$Y_d$ & $-\bar QY_dHd+\mathrm{h.c.}$ \\",
        r"$Y_e$ & $-\bar LY_eHe+\mathrm{h.c.}$ \\",
        r"$y_1$ & $LFS_1+\mathrm{h.c.}$ \\",
        r"$y_2$ & $LFS_2+\mathrm{h.c.}$ \\",
        r"$M_F$ & $-\bar FM_FF$ \\",
        r"$T$ & $\operatorname{Tr}(Y_e^\dagger Y_e+3Y_u^\dagger Y_u+3Y_d^\dagger Y_d)$ \\",
        r"$T_\nu^{(1)}$ & $\operatorname{Tr}(y_1^\dagger y_1)$ \\",
        r"$T_\nu^{(2)}$ & $\operatorname{Tr}(y_2^\dagger y_2)$ \\",
        r"$T_\nu$ & $T_\nu^{(1)}+T_\nu^{(2)}$ \\",
        r"$m_1^2$ & $m_1^2S_1^\dagger S_1$ \\",
        r"$m_2^2$ & $m_2^2S_2^\dagger S_2$ \\",
        r"$\lambda_1$ & $\frac12\lambda_1(H^\dagger H)^2$ \\",
        r"$\lambda_{S_1}^{(1)}$ & $\frac12\lambda_{S_1}^{(1)}(S_1^\dagger S_1)^2$ \\",
        r"$\lambda_{S_2}^{(1)}$ & $\frac12\lambda_{S_2}^{(1)}(S_2^\dagger S_2)^2$ \\",
        r"$\lambda_{HS_1}^{(1)}$ & $\lambda_{HS_1}^{(1)}(H^\dagger H)(S_1^\dagger S_1)$ \\",
        r"$\lambda_{HS_2}^{(1)}$ & $\lambda_{HS_2}^{(1)}(H^\dagger H)(S_2^\dagger S_2)$ \\",
        r"$\lambda_{HS_1}^{(A)}$ & $\lambda_{HS_1}^{(A)}(H^\dagger T^AH)(S_1^\dagger T^AS_1)$ \\",
        r"$\lambda_{HS_2}^{(A)}$ & $\lambda_{HS_2}^{(A)}(H^\dagger T^AH)(S_2^\dagger T^AS_2)$ \\",
        r"$\lambda_5$ & $\lambda_5HHS_1S_2^\dagger+\mathrm{h.c.}$ \\",
        r"$\lambda_{12}^{(1)}$ & $\lambda_{12}^{(1)}(S_1^\dagger S_1)(S_2^\dagger S_2)$ \\",
        r"$\lambda_{12}^{(A)}$ & $\lambda_{12}^{(A)}(S_1^\dagger T^AS_1)(S_2^\dagger T^AS_2)$ \\",
        r"$\lambda_{12}^{(\times)}$ & crossed independent $S_1$--$S_2$ contraction \\",
        r"$\lambda_{S_1}^{(A)}$ & additional independent $S_1$ self-contraction \\",
        r"$\lambda_{S_2}^{(A)}$ & additional independent $S_2$ self-contraction \\",
        r"$\alpha$ & $Y_{S_1}=\alpha/2,\ Y_F=(\alpha+1)/2,\ Y_{S_2}=(\alpha+2)/2$ \\",
        r"$C_2(d)$ & $(d^2-1)/4$ \\",
        r"$T(d)$ & $d(d^2-1)/12$ \\",
        r"$C_{LLS_1S_2}^{ij;AB}$ & $L_iL_jS_1^AS_2^B$ \\",
        r"$C_{ijab}$ & $L_iL_j\phi_a\phi_b$ \\",
        r"$C_5^{ij}$ & $(L_i^TCL_j)HH$ \\",
        r"$R_W$ & $16\pi^2\beta_{C_5}\supset R_W\lambda_5C_{LLS_1S_2}$ \\",
        r"\bottomrule",
        r"\end{longtable}",
    ]


def apply_paper_notation(latex: str) -> str:
    """Convert rendered implementation symbols to common report notation."""
    if not latex:
        return latex

    text = str(latex)

    for raw, pretty in RENDERED_PAPER_REPLACEMENTS.items():
        text = text.replace(raw, pretty)

    for raw in sorted(PAPER_SYMBOL_LATEX, key=len, reverse=True):
        pretty = PAPER_SYMBOL_LATEX[raw]
        text = text.replace(r"\text{" + raw + "}", "{" + pretty + "}")
        text = re.sub(
            rf"(?<![A-Za-z0-9]){re.escape(raw)}(?![A-Za-z0-9])",
            lambda _match, replacement=pretty: replacement,
            text,
        )

    return text


# ---------------------------------------------------------------------------
# Stage-aware report paths
# Consolidated from the former Reports/StageReports.py module.
# ---------------------------------------------------------------------------

def final_eft_stage_label(records: Iterable[object]) -> str:
    """Return the final recorded physical EFT-stage label.

    Stage labels come from each run's physical threshold metadata. If no
    stage metadata is available, use the neutral fallback ``EFT`` rather than
    assuming the ordinary-T3 ``S1,S2`` field content.
    """
    for record in records:
        stages = getattr(record, "eft_stages", None)
        if stages:
            return stages[-1].label

    return "EFT"


def lagrangian_report_path(
    stage_label: str,
    report_root: Path | None = None,
) -> Path:
    root = report_root or REPORT_OUTPUT_DIR
    return root / "Lagrangian" / f"{stage_label}.tex"


def rge_report_path(
    stage_label: str,
    report_root: Path | None = None,
) -> Path:
    root = report_root or REPORT_OUTPUT_DIR
    return root / "RGE" / f"{stage_label}.tex"


def group_factor_report_path(
    stage_label: str,
    report_root: Path | None = None,
) -> Path:
    root = report_root or REPORT_OUTPUT_DIR
    return root / "GroupFactors" / f"GF_{stage_label}.tex"


def c5_report_path(report_root: Path | None = None) -> Path:
    root = report_root or REPORT_OUTPUT_DIR
    return root / "Lagrangian" / "C5.tex"


def report_output_dir_for(record: RunRecord) -> Path:
    # Human-readable per-model reports are kept away from raw calculation data.
    output_dir = REPORT_OUTPUT_DIR / record.output_dir.name
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


# ---------------------------------------------------------------------------
# LaTeX/report helpers
# ---------------------------------------------------------------------------


def normalise_physics_latex(latex: str) -> str:
    """Clean implementation/Matchete notation for human-facing reports.

    This is deliberately a presentation-only pass: it does not alter the
    underlying matching output.  It maps Matchete/internal symbol names to the
    notation used in the report and keeps physically meaningful conjugation
    explicit.
    """
    if not latex:
        return latex

    text = latex

    # Renormalisation scale. Mathematica/TeXForm can render the symbol mubar2
    # as the particularly ugly ``\\text{$\\mu $bar2}``.
    text = text.replace(r"\text{$\mu $bar2}", r"\bar{\mu}^{2}")
    text = text.replace(r"\text{mubar2}", r"\bar{\mu}^{2}")
    text = text.replace(r"\mathrm{mubar2}", r"\bar{\mu}^{2}")

    # Matchete coefficients generated for the scalar bilinears after the first
    # threshold.  Name them by the operator they multiply rather than exposing
    # the internal CNewScalar12/CNewScalar22 symbols.  Do this before replacing
    # the field names so we do not accidentally modify the coefficient name.
    text = text.replace(r"\text{CNewScalar12}", r"C_{S_1^\dagger S_1}")
    text = text.replace(r"\text{CNewScalar22}", r"C_{S_2^\dagger S_2}")
    text = re.sub(
        r"(?<![A-Za-z0-9])CNewScalar12(?![A-Za-z0-9])",
        lambda _match: r"C_{S_1^\dagger S_1}",
        text,
    )
    text = re.sub(
        r"(?<![A-Za-z0-9])CNewScalar22(?![A-Za-z0-9])",
        lambda _match: r"C_{S_2^\dagger S_2}",
        text,
    )

    # Generic heavy-field implementation names.
    text = text.replace(r"\text{NewScalar1}", r"S_1")
    text = text.replace(r"\text{NewScalar2}", r"S_2")
    text = re.sub(
        r"(?<![A-Za-z0-9])NewScalar1(?![A-Za-z0-9])",
        r"S_1",
        text,
    )
    text = re.sub(
        r"(?<![A-Za-z0-9])NewScalar2(?![A-Za-z0-9])",
        r"S_2",
        text,
    )

    # T3 Clebsch--Gordan / invariant tensors.  Bar[...] is not discarded: it
    # denotes the complex-conjugate invariant tensor.  Render it as a star.
    cg_map = {
        "T3Y1CG": r"\mathcal{I}_{y_1}",
        "T3Y2CG": r"\mathcal{I}_{y_2}",
        "T3MixCG": r"\mathcal{I}_{T3}",
    }
    for raw_name, pretty in cg_map.items():
        text = text.replace(
            r"\text{Bar}(\text{" + raw_name + "})",
            "{" + pretty + r"}^{*}",
        )
        text = text.replace(r"\text{" + raw_name + "}", pretty)

    # The SU(2) epsilon tensor is real in the convention used here, so the
    # conjugated Matchete form does not need a visible bar.
    text = text.replace(
        r"\text{Bar}(\text{eps}(\text{SU2L}))",
        r"\epsilon",
    )
    text = text.replace(r"\text{eps}(\text{SU2L})", r"\epsilon")

    # Some tensors have already been partially prettified by the Wolfram
    # exporter.  Remove remaining \text{} wrappers from their labels.
    invariant_label_map = {
        "y1": "y_1",
        "y2": "y_2",
        "T3": "T3",
        "H1": "H1",
        "H2": "H2",
        "S1": "S_1",
        "S2": "S_2",
    }
    for raw_label, pretty_label in invariant_label_map.items():
        text = text.replace(
            r"\mathcal{I}_{\text{" + raw_label + "},",
            r"\mathcal{I}_{" + pretty_label + ",",
        )
        text = text.replace(
            r"\mathcal{I}_{\text{" + raw_label + "}}",
            r"\mathcal{I}_{" + pretty_label + "}",
        )

    # Scalar masses.
    text = text.replace(r"M_{\text{S1}}", r"m_1")
    text = text.replace(r"M_{\text{S2}}", r"m_2")

    # Scalar couplings.  The trailing integer labels independent invariant
    # contractions, so display it as (n), not as an algebraic power.
    def _lambda_text_repl(match: re.Match[str]) -> str:
        label = match.group(1)
        invariant = match.group(2)
        label_map = {"S1": "S_1", "S2": "S_2"}
        pretty = label_map.get(label, label)
        result = rf"\lambda_{{{pretty}}}"
        if invariant is not None:
            result += rf"^{{({invariant})}}"
        return result

    text = re.sub(
        r"\\lambda\s*_\{\\text\{([^}]+)\}\}(?:\{\}\^(\d+))?",
        _lambda_text_repl,
        text,
    )

    def _lambda_plain_repl(match: re.Match[str]) -> str:
        label = match.group(1)
        invariant = match.group(2)
        result = rf"\lambda_{{{label}}}"
        if invariant is not None:
            result += rf"^{{({invariant})}}"
        return result

    text = re.sub(
        r"\\lambda\s*_\{(12)\}(?:\{\}\^(\d+))?",
        _lambda_plain_repl,
        text,
    )

    return apply_paper_notation(text)

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

    if record.shared_scalar:
        return (
            rf"{record.name}:\quad "
            rf"S=({d_s2},{latex_fraction(y_s2)}),\quad "
            rf"\widetilde S=({d_s1},{latex_fraction(y_s1)}),\quad "
            rf"F=({d_f},{latex_fraction(y_f)})"
        )

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
    (
        r"S_1",
        (
            r"\\text\{NewScalar1\}",
            r"S_\{1(?:,[^}]*)?\}",
            r"S_1(?![A-Za-z0-9])",
        ),
    ),
    (
        r"S_2",
        (
            r"\\text\{NewScalar2\}",
            r"S_\{2(?:,[^}]*)?\}",
            r"S_2(?![A-Za-z0-9])",
        ),
    ),
    (
        r"F",
        (
            r"N(?:_|\^)",
            r"(?<![A-Za-z\\])F(?:_|\^)",
        ),
    ),
    (r"H", (r"(?<![A-Za-z\\])H(?:_|\^)",)),
    (
        r"\ell",
        (
            r"\\ell",
            r"\\mathcal\{l\}",
            r"\\mathscr\{l\}",
            r"(?<![A-Za-z\\])L(?:_|\^)",
        ),
    ),
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

    latex = normalise_physics_latex(latex)
    grouped: dict[tuple[str, ...], list[str]] = {}

    for term in split_latex_terms(latex):
        signature = latex_field_signature(term)

        grouped.setdefault(signature, []).append(term)

    return grouped


def field_signature_label(signature: tuple[str, ...]) -> str:
    """Render a field combination as a LaTeX table heading."""

    return "$" + r"\,".join(signature) + "$"


FIELD_SIGNATURE_ORDER = {
    r"H": 0,
    r"S_1": 1,
    r"S_2": 2,
    r"F": 3,
    r"\ell": 4,
    r"e": 5,
    r"q": 6,
    r"u": 7,
    r"d": 8,
    r"B_{\mu\nu}": 9,
    r"W_{\mu\nu}": 10,
    r"G_{\mu\nu}": 11,
    r"D": 12,
    r"\text{constant}": 99,
}


def field_signature_sort_key(signature: tuple[str, ...]) -> tuple:
    """Sort columns by physical field content; the synthetic constant bucket is last."""
    if signature == (r"\text{constant}",):
        return (99, 99, ())
    physical_fields = tuple(field for field in signature if field != r"D")
    ranks = tuple(FIELD_SIGNATURE_ORDER.get(field, 50) for field in physical_fields)
    return (len(physical_fields), ranks, int(r"D" in signature), signature)


def _term_chunks(terms: list[str], max_terms: int = 3) -> list[list[str]]:
    """Split long cells into continuation rows so longtable can page-break."""
    return [terms[i:i + max_terms] for i in range(0, len(terms), max_terms)]


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
        lines.extend(latex_aligned_block(normalise_physics_latex(latex)))
    else:
        lines.append(r"\textit{Not available for this model.}")


def write_lagrangian_report(records: list[RunRecord]) -> Path:
    """Create the concise human-facing Lagrangian and Weinberg report."""

    output_path = REPORT_OUTPUT_DIR / "lagrangian_report.tex"

    lines: list[str] = [
        r"\documentclass[10pt]{article}",
        r"\usepackage[margin=1.7cm]{geometry}",
        r"\usepackage{amsmath,amssymb,adjustbox,booktabs,array,longtable}",
        r"\usepackage[T1]{fontenc}",
        r"\allowdisplaybreaks[4]",
        r"\setlength{\emergencystretch}{3em}",
        r"\begin{document}",
        r"\section*{T3 Lagrangian and Weinberg-operator report}",
        *paper_notation_key_lines(),
        rf"EFT order: ${EFT_ORDER}$; loop order: ${LOOP_ORDER}$.",
        (
            r"Implementation names are translated into physics notation: "
            r"$S_1,S_2,F$ denote the new fields and "
            r"$\mathcal{I}^{(n)}$ denotes the $n$th independent gauge-invariant "
            r"tensor contraction when more than one invariant exists."
        ),
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
        r"\usepackage{amsmath,amssymb,adjustbox,pdflscape,longtable}",
        r"\usepackage[T1]{fontenc}",
        r"\allowdisplaybreaks[4]",
        r"\setlength{\emergencystretch}{3em}",
        r"\begin{document}",
        r"\section*{T3 UV and EFT Lagrangians}",
        *paper_notation_key_lines(),
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
                *latex_aligned_block(normalise_physics_latex(uv_latex)),
                r"\textbf{Matched EFT Lagrangian}",
                *latex_aligned_block(normalise_physics_latex(eft_latex)),
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
        r"\usepackage{amsmath,amssymb,adjustbox,pdflscape,longtable}",
        r"\usepackage[T1]{fontenc}",
        r"\allowdisplaybreaks[4]",
        r"\setlength{\emergencystretch}{3em}",
        r"\begin{document}",
        r"\section*{T3 BSM UV and matched EFT contributions}",
        *paper_notation_key_lines(),
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
            lines.extend(latex_aligned_block(normalise_physics_latex(uv_latex)))
        else:
            lines.append(r"\textit{BSM UV Lagrangian unavailable.}")

        lines.append(
            r"\textbf{$\Delta\mathcal{L}_{\mathrm{EFT}}$}"
        )

        if eft_ok and eft_latex:
            lines.extend(latex_aligned_block(normalise_physics_latex(eft_latex)))
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

def _study_lagrangian_report_path(stage_label: str, report_root: Path | None = None) -> Path:
    """Return the Lagrangian report path, optionally rooted in one study directory.

    ``pipeline.finish_runs`` passes ``report_root`` for study-aware output such as
    ``Reports/output/single`` or ``Reports/output/full/hypercharge``.  Keep the
    historical StageReports helper as the fallback for callers that do not pass
    a study root.
    """
    if report_root is None:
        return lagrangian_report_path(stage_label)
    return Path(report_root) / "Lagrangian" / f"{stage_label}.tex"

def write_bsm_field_table(
    records: list[RunRecord],
    *,
    summary_key: str,
    output_path: Path,
    title: str,
    description: str,
    empty_value: str = "",
    fallback_summary_key: str | None = None,
    require_matching: bool = True,
) -> Path:
    """Create a comparison table grouping Lagrangian terms by field content.

    The preferred source is ``summary_key`` (normally the BSM-only expression).
    If that expression is unavailable, ``fallback_summary_key`` is used so that
    report generation does not silently become empty merely because one optional
    conversion flag or BSM subtraction failed.

    UV reports only require a successful build. EFT reports additionally require
    successful matching.
    """

    output_path.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    signatures: set[tuple[str, ...]] = set()

    for record in records:
        summary = record.summary

        if summary.get("BuildStatus") != "Success":
            continue

        if require_matching and summary.get("MatchingStatus") != "Success":
            continue

        latex = str(summary.get(summary_key, "") or "").strip()

        if not latex and fallback_summary_key:
            latex = str(summary.get(fallback_summary_key, "") or "").strip()

        if not latex:
            continue

        grouped = grouped_lagrangian_terms(latex)

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
        key=field_signature_sort_key,
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
        *paper_notation_key_lines(),
        r"\tiny",
    ]

    if not rows or not ordered_signatures:
        lines.append(
            r"No Lagrangian expressions were available for this sector."
        )
    else:
        # Keep only a few field-content columns on each page so the expressions
        # remain readable rather than being compressed into an unusable table.
        chunk_size = 3

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
                r">{\raggedright\arraybackslash}p{0.21\linewidth}"
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

                per_signature_chunks = {
                    signature: _term_chunks(grouped.get(signature, []))
                    for signature in signature_chunk
                }
                continuation_rows = max(
                    [len(parts) for parts in per_signature_chunks.values()] + [1]
                )

                for continuation_index in range(continuation_rows):
                    if continuation_index == 0:
                        metadata_cells = [
                            latex_escape_text(record.name),
                            rf"${record.alpha}$",
                            rf"${d_s1}$",
                            rf"${latex_fraction(y_s1)}$",
                            rf"${d_s2}$",
                            rf"${latex_fraction(y_s2)}$",
                            rf"${d_f}$",
                            rf"${latex_fraction(y_f)}$",
                        ]
                    else:
                        metadata_cells = [r"\textit{cont.}", "", "", "", "", "", "", ""]

                    expression_cells = []
                    for signature in signature_chunk:
                        parts = per_signature_chunks[signature]
                        if continuation_index < len(parts):
                            expression_cells.append(
                                matrix_cell(parts[continuation_index], empty_value)
                            )
                        else:
                            expression_cells.append(empty_value)

                    lines.append(
                        " & ".join(metadata_cells + expression_cells) + r" \\"
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

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        "\n".join(lines),
        encoding="utf-8",
    )

    print(f"\nLagrangian comparison table:\n{output_path}")

    return output_path

def write_bsm_uv_field_table(
    records: list[RunRecord],
    empty_value: str = "",
    *,
    report_root: Path | None = None,
) -> Path:
    """Create the UV Lagrangian comparison table by field content."""

    return write_bsm_field_table(
        records,
        summary_key="BSMUVLagrangianLaTeX",
        fallback_summary_key="UVLagrangianLaTeX",
        require_matching=False,
        output_path=_study_lagrangian_report_path("UV", report_root),
        title="T3 UV Lagrangian terms grouped by field content",
        description=(
            r"Rows are model configurations and columns are field configurations. "
            r"The BSM-only UV expression is preferred; the full UV Lagrangian is "
            r"used as a fallback when the BSM-only conversion is unavailable."
        ),
        empty_value=empty_value,
    )

def write_bsm_matched_field_table(
    records: list[RunRecord],
    empty_value: str = "",
    *,
    report_root: Path | None = None,
) -> Path:
    """Create one Lagrangian comparison report for every matched EFT stage."""

    stage_labels: list[str] = []

    for record in records:
        for stage in record.summary.get("EFTStages", []):
            label = str(stage.get("Label", "")).strip()
            if label and label not in stage_labels:
                stage_labels.append(label)

    if not stage_labels:
        # Historical fallback for old summaries.
        stage_labels = [final_eft_stage_label(records)]

    last_path = _study_lagrangian_report_path(stage_labels[-1], report_root)

    for stage_label in stage_labels:
        stage_records: list[RunRecord] = []

        for record in records:
            stage_payload = next(
                (
                    stage
                    for stage in record.summary.get("EFTStages", [])
                    if stage.get("Label") == stage_label
                ),
                None,
            )

            if stage_payload is None:
                continue

            stage_summary = dict(record.summary)
            stage_summary["BSMEFTLagrangianLaTeX"] = stage_payload.get(
                "BSMEFTLagrangianLaTeX",
                "",
            )
            stage_summary["EFTLagrangianLaTeX"] = stage_payload.get(
                "EFTLagrangianLaTeX",
                "",
            )
            stage_summary["MatchingStatus"] = (
                "Success"
                if record.summary.get("SequentialMatchingStatus") == "Success"
                else record.summary.get("MatchingStatus", "Unknown")
            )

            stage_records.append(
                RunRecord(
                    name=record.name,
                    alpha=record.alpha,
                    d_s1=record.d_s1,
                    d_s2=record.d_s2,
                    d_f=record.d_f,
                    return_code=record.return_code,
                    summary=stage_summary,
                    output_dir=record.output_dir,
                    eft_stages=record.eft_stages,
                )
            )

        if not stage_records:
            continue

        stage_info = next(
            (
                stage
                for stage in stage_records[0].summary.get("EFTStages", [])
                if stage.get("Label") == stage_label
            ),
            {},
        )
        active = stage_info.get("ActiveHeavyFields", [])
        integrated = stage_info.get("IntegratedFields", [])

        description = (
            r"Rows are model configurations and columns are field configurations. "
            + "This stage integrates out "
            + ", ".join(integrated)
            + ". Remaining active T3 fields: "
            + (", ".join(active) if active else "none")
            + "."
        )

        last_path = write_bsm_field_table(
            stage_records,
            summary_key="BSMEFTLagrangianLaTeX",
            fallback_summary_key="EFTLagrangianLaTeX",
            require_matching=True,
            output_path=_study_lagrangian_report_path(stage_label, report_root),
            title=(
                "T3 matched EFT Lagrangian terms: "
                + stage_label.replace("_", r"\_")
            ),
            description=description,
            empty_value=empty_value,
        )

        compile_latex_document(last_path)

    return last_path


def _c5_display_indices(latex: str) -> str:
    """Rename Matchete's generated free/dummy labels for human-facing C5 output."""
    text = str(latex)
    replacements = {
        "i_1": "i",
        "i_2": "j",
        "r_1": "r",
    }
    for raw, pretty in replacements.items():
        text = text.replace(raw, pretty)
    return text


def _swap_c5_flavour_indices(latex: str) -> str:
    r"""Return the same ordered coefficient with only flavour labels i <-> j.

    Do not replace bare letters globally: that would corrupt LaTeX commands such
    as ``\right`` and ``\overline``.  Here the free flavour labels occur only as
    comma-delimited indices inside the displayed Yukawa subscripts.
    """
    text = _c5_display_indices(latex)

    # Swap only index tokens, using a temporary placeholder.
    text = re.sub(r"(?<=,)i(?=[,}])", "__C5_I__", text)
    text = re.sub(r"(?<=,)j(?=[,}])", "i", text)
    text = text.replace("__C5_I__", "j")

    return text


def write_c5_coefficient_report(
    records: list[RunRecord],
    *,
    report_root: Path | None = None,
) -> Path:
    """Create a compact study-level C5 comparison table.

    The report deliberately contains only the ordered coefficient extracted
    from Matchete and the corresponding symmetric physical Wilson coefficient.
    """

    output_path = c5_report_path(report_root=report_root)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    rows: list[tuple[RunRecord, str]] = []

    for record in records:
        summary = record.summary

        coefficient_latex = str(
            summary.get("WeinbergCanonicalCoefficientLaTeX")
            or summary.get("CanonicalCoefficientLaTeX")
            or summary.get("WeinbergCoefficientLaTeX")
            or ""
        ).strip()

        # Presentation-only consistency: use uppercase I_3 in the report.
        coefficient_latex = coefficient_latex.replace(
            r"i_3\left(",
            r"I_3\left(",
        )
        coefficient_latex = _c5_display_indices(coefficient_latex)
        coefficient_latex = normalise_physics_latex(coefficient_latex)

        coefficient_ok = (
            summary.get("WeinbergCanonicalCoefficientConversionSuccess") is True
            or summary.get("CanonicalCoefficientConversionSuccess") is True
            or summary.get("WeinbergCoefficientConversionSuccess") is True
        )

        if coefficient_ok and coefficient_latex:
            rows.append((record, coefficient_latex))

    lines: list[str] = [
        r"\documentclass[8pt]{article}",
        r"\usepackage[margin=0.7cm]{geometry}",
        r"\usepackage{amsmath,amssymb,adjustbox,longtable,array,booktabs,pdflscape}",
        r"\usepackage[T1]{fontenc}",
        r"\setlength{\tabcolsep}{4pt}",
        r"\renewcommand{\arraystretch}{1.35}",
        r"\begin{document}",
        r"\begin{landscape}",
        r"\section*{T3 Weinberg-operator coefficient comparison}",
        *paper_notation_key_lines(),
    ]

    if not rows:
        lines.append(
            r"\textit{No successfully converted Weinberg coefficient was available.}"
        )
    else:
        lines.extend(
            [
                r"\scriptsize",
                r"\begin{longtable}{@{}l "
                r">{\raggedright\arraybackslash}p{0.42\linewidth} "
                r">{\raggedright\arraybackslash}p{0.42\linewidth}@{}}",
                r"\toprule",
                (
                    r"Model & Ordered coefficient extracted from Matchete "
                    r"& Symmetric physical Wilson coefficient \\"
                ),
                r"\midrule",
                r"\endfirsthead",
                r"\toprule",
                (
                    r"Model & Ordered coefficient extracted from Matchete "
                    r"& Symmetric physical Wilson coefficient \\"
                ),
                r"\midrule",
                r"\endhead",
            ]
        )

        for record, coefficient_latex in rows:
            ordered_cell = (
                r"\(\displaystyle A_5^{ij}="
                + coefficient_latex
                + r"\)"
            )
            # Factor out the common mass/loop/coupling prefactor so the
            # symmetric coefficient is compact and easy to compare across models.
            ordered_flavour = (
                r"\overline{y_{1,i,r}} \overline{y_{2,j,r}}"
            )
            swapped_flavour = (
                r"\overline{y_{1,j,r}} \overline{y_{2,i,r}}"
            )

            common_prefactor = coefficient_latex.replace(
                ordered_flavour,
                "",
                1,
            ).strip()

            symmetric_cell = (
                r"\(\displaystyle C_5^{ij}="
                + common_prefactor
                + r"\left("
                + ordered_flavour
                + r"+"
                + swapped_flavour
                + r"\right)\)"
            )

            lines.append(
                " & ".join(
                    [
                        latex_escape_text(record.name),
                        ordered_cell,
                        symmetric_cell,
                    ]
                )
                + r" \\"
            )
            lines.append(r"\midrule")

        lines.extend(
            [
                r"\bottomrule",
                r"\end{longtable}",
            ]
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

    print(f"\nC5 coefficient report:\n{output_path}")
    return output_path

def compile_latex_document(tex_path: Path) -> None:
    """Compile a generated LaTeX report using latexmk or pdflatex."""

    latexmk = shutil.which("latexmk")
    pdflatex = shutil.which("pdflatex")

    # MiKTeX's latexmk wrapper requires Perl.  On Windows, if Perl is not
    # installed, calling latexmk only produces a noisy failure before we fall
    # back to pdflatex anyway.  Skip latexmk in that situation.
    if os.name == "nt" and shutil.which("perl") is None:
        latexmk = None

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
    *,
    report_root: Path | None = None,
) -> None:
    """Generate all Lagrangian reports and term tables."""

    report_tex = write_lagrangian_report(records)
    compile_latex_document(report_tex)

    uv_table_tex = write_bsm_uv_field_table(
        records,
        report_root=report_root,
    )
    compile_latex_document(uv_table_tex)

    matched_table_tex = write_bsm_matched_field_table(
        records,
        report_root=report_root,
    )
    compile_latex_document(matched_table_tex)

    coefficient_tex = write_c5_coefficient_report(
        records,
        report_root=report_root,
    )
    compile_latex_document(coefficient_tex)

    if debug_reports:
        full_tex = write_latex_lagrangian_table(records)
        compile_latex_document(full_tex)

        bsm_tex = write_bsm_lagrangian_table(records)
        compile_latex_document(bsm_tex)


