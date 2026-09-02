from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

from RGE.matching.MatchedEFTRGE import run_matched_eft_rge
from RGE.stages.FlavorMatchedRGEStage import run_flavor_matched_rge
from RGE.stages.NeutrinoMassStage import run_neutrino_mass_stage
from RGE.stages.NumericalPipelineStage import run_numerical_pipeline_stage
from RGE.stages.RGEReportStage import write_rge_report
from RGE.phenomenology.NeutrinoObservables import run_neutrino_observables_stage

PROJECT_ROOT = Path(__file__).resolve().parent

LAGRANGIAN_DIR = PROJECT_ROOT / "Lagrangian"
RGE_DIR = PROJECT_ROOT / "RGE"
OUTPUT_DIR = PROJECT_ROOT / "output"

RUN_MODEL_SCRIPT = LAGRANGIAN_DIR / "RunModel.wl"

# Weinberg 5D at 1 loop
EFT_ORDER = 5
LOOP_ORDER = 1


# Known T3 models from the original classification with specific dimensions
T3_CLASSES = {
    "A": (1, 3, 2),
    "B": (2, 2, 1),
    "C": (2, 2, 3),
    "D": (3, 1, 2),
    "E": (3, 3, 2),
}


# Some interesting models to look into.
# The first entry is the known T3 class and the second is alpha
# Hypercharges are fixed by gauge invariance of the Yukawa interactions
# and the four-scalar interaction:
#   Y(S1) = alpha
#   Y(F)  = alpha + 1
#   Y(S2) = alpha + 2
# The A-E labels only specify the SU(2) representations.
INTERESTING = [("A", 0), ("B", -1), ("C", -1), ("D", -2), ("E", 0)]

SMOKE = [("B", -1), ("C", -1), ("A", 0), ("D", -2), ("E", 0)]

EXTENDED = [
    ("A", 0),
    ("A", -2),
    ("B", -1),
    ("C", -1),
    ("D", -2),
    ("E", 0),
    ("E", -2),
]


# For each completed run we have this object
@dataclass
class RunRecord:
    name: str
    alpha: int
    d_s1: int
    d_s2: int
    d_f: int
    return_code: int
    summary: dict
    output_dir: Path


def encode_alpha(alpha: int) -> str:
    """For Wolfram, recode +5 -> p5 and -5 -> m5."""
    return f"m{abs(alpha)}" if alpha < 0 else f"p{alpha}"


def valid_t3_dimensions(d_s1: int, d_s2: int, d_f: int) -> bool:
    """We check if our SU2 dimensions for our fields are valid for T3"""
    # Only positive dimensions
    if min(d_s1, d_s2, d_f) < 1:
        return False

    # We check for yukawa interaction, knowing that it interacts with a lepton doublet. Thus we have 2⊗dF = 1/2 ⊗ jF -> (1/2+jF)+(1/2-jF)
    # Note that we need j = 0, so when we have j1 ⊗ j2 = j1-j2 +.... j1+j2, we need j1=j2 for us to get a singlet
    # Thus we have dS=1+dF or ds = 1-dF, which is what this part checks
    if abs(d_s1 - d_f) != 1 or abs(d_s2 - d_f) != 1:
        return False

    # d=2j+1; apply the usual SU(2) angular-momentum addition rule for J=1.
    j1 = (d_s1 - 1) / 2
    j2 = (d_s2 - 1) / 2

    # Now we check for scalars with HH part, where we know that our Higgs are doublets, giving us 2x2=3+1 for higgs. Higgs are the same field so our antisymmetric singlet is not included
    # So we only have dimension 3 for our 2 Higgs and isospin charge of J=1. 
    # As isospin is kinda associative and commutative, we just need to match j1+j2 and HH, 
    # This we need to check if j1 x j2 = j1-j2 + ... ,j1+j2 contains a j=1 term to match HH, where we need j1+j2 to be an integer so that it would have j=1 
    return abs(j1 - j2) <= 1 <= j1 + j2 and float(j1 + j2).is_integer()



def identify_t3_class(d_s1: int, d_s2: int, d_f: int) -> str | None:
    """Return the known A-E label if these dimensions match one."""

    dimensions = (d_s1, d_s2, d_f)

    for model_class, known_dimensions in T3_CLASSES.items():
        if dimensions == known_dimensions:
            return model_class

    return None


def run_model(
    name: str,
    alpha: int,
    d_s1: int,
    d_s2: int,
    d_f: int,
    output_dir: Path,
    model_args: list[str],
    debug_reports: bool = False,
) -> RunRecord:
    """What this does is 
    1. Delete previous output directory and recreate for new results
    2. Run Runmodel.wl with out inputs
    3. Get its output and errors into a file
    4. Get debug reports and summaries
    5. Return a RunRecord object with all the data."""

    # Delete the previous output directory and recreate it.
    # This prevents an old successful result being mistaken for a new result
    # if the current Wolfram run fails before producing its summary.
    shutil.rmtree(output_dir, ignore_errors=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Run the model with inputs:
    #   output directory
    #   EFT order
    #   loop order
    #   model arguments
    # We use RunModel.wl+
    command = [
        "wolframscript",
        "-file",
        str(RUN_MODEL_SCRIPT),
        str(output_dir),
        str(EFT_ORDER),
        str(LOOP_ORDER),
        *model_args,
        *(["DEBUG"] if debug_reports else []),
    ]

    # Launch Wolfram and capture both normal output and errors as text.
    process = subprocess.run(
        command,
        cwd=PROJECT_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )

    # Save Wolfram output for debugging.
    if debug_reports:
        debug_dir = output_dir / "debug"
        debug_dir.mkdir(parents=True, exist_ok=True)
        (debug_dir / "wolfram_stdout.log").write_text(
            process.stdout,
            encoding="utf-8",
        )
        (debug_dir / "wolfram_stderr.log").write_text(
            process.stderr,
            encoding="utf-8",
        )
    elif process.returncode != 0:
        (output_dir / "run_error.log").write_text(
            process.stdout + "\n" + process.stderr,
            encoding="utf-8",
        )

    # If Wolfram fails, show only the final part of its output so the
    # terminal remains readable while still giving useful debug information.
    if process.returncode != 0:
        if process.stdout:
            print("\n".join(process.stdout.splitlines()[-40:]))

        if process.stderr:
            print(
                "\n".join(process.stderr.splitlines()[-20:]),
                file=sys.stderr,
            )

    # Read the summary produced by the Wolfram side.
    summary_path = output_dir / "comparison_summary.json"

    summary = (
        json.loads(summary_path.read_text(encoding="utf-8"))
        if summary_path.exists()
        else {
            "BuildStatus": "ProcessFailed",
            "MatchingStatus": "NotRun",
        }
    )

    # Put everything associated with this run into one object.
    return RunRecord(
        name=name,
        alpha=alpha,
        d_s1=d_s1,
        d_s2=d_s2,
        d_f=d_f,
        return_code=process.returncode,
        summary=summary,
        output_dir=output_dir,
    )


def run_dimensions(
    d_s1: int,
    d_s2: int,
    d_f: int,
    alpha: int,
    debug_reports: bool = False,
) -> RunRecord:
    """All it does is 
    1. Check if dimensions are correct, if not then return error
    2. Identify if its an T3-A..E model and name the output folder after it
    3. Use run_model"""

    # First check whether the dimensions can form the required T3
    # Yukawa and scalar interactions.
    if not valid_t3_dimensions(d_s1, d_s2, d_f):
        raise ValueError(
            f"({d_s1}, {d_s2}, {d_f}) is not a valid T3 SU(2) assignment: "
            "each scalar must have dS=dF±1 and S1⊗S2 must contain the triplet."
        )

    # Check whether these dimensions correspond to one of the known
    # T3-A ... T3-E models from the original classification. (Its for output report names)
    model_class = identify_t3_class(d_s1, d_s2, d_f)

    if model_class is not None:
        # Known model: keep its familiar A-E name.
        name = f"T3-{model_class}"
        output_dir = OUTPUT_DIR / (
            f"T3_{model_class}_alpha_{encode_alpha(alpha)}"
        )
    else:
        # New/generalised representation: identify it directly by dimensions.
        name = f"T3-d{d_s1}-d{d_s2}-F{d_f}"
        output_dir = OUTPUT_DIR / (
            f"T3_d{d_s1}_d{d_s2}_F{d_f}_alpha_{encode_alpha(alpha)}"
        )

    # Give arguments for the alpha and the dimensions
    model_args = [
        "DIMS",
        str(d_s1),
        str(d_s2),
        str(d_f),
        encode_alpha(alpha),
    ]

    print(
        f"Running {name}, "
        f"dims=({d_s1}, {d_s2}, {d_f}), "
        f"alpha={alpha} ...",
        flush=True,
    )

    return run_model(
        name,
        alpha,
        d_s1,
        d_s2,
        d_f,
        output_dir,
        model_args,
        debug_reports,
    )


def run_known_class(
    model_class: str,
    alpha: int,
    debug_reports: bool = False,
) -> RunRecord:
    """Convert a known A-E benchmark into dimensions and run normally."""

    if model_class not in T3_CLASSES:
        raise ValueError(f"Unknown T3 model class: {model_class}")

    d_s1, d_s2, d_f = T3_CLASSES[model_class]

    return run_dimensions(
        d_s1,
        d_s2,
        d_f,
        alpha,
        debug_reports,
    )


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

    output_path = OUTPUT_DIR / "lagrangian_report.tex"

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

    output_path = OUTPUT_DIR / "eft_lagrangians.tex"

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

    output_path = OUTPUT_DIR / "bsm_lagrangians.tex"

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

    output_path = OUTPUT_DIR / f"{output_stem}.tex"

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

    output_path = record.output_dir / "c5_coefficient.tex"
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


def print_summary(records: list[RunRecord]) -> int:
    """Print the results of all completed T3 runs."""

    print(
        "\n"
        + "=" * 72
        + "\nT3 MODEL SUMMARY\n"
        + "=" * 72
    )

    successful = 0

    for record in records:
        summary = record.summary

        build_ok = summary.get("BuildStatus") == "Success"
        match_ok = summary.get("MatchingStatus") == "Success"

        successful += int(build_ok and match_ok)

        print(
            f"{record.name} alpha={record.alpha}: "
            f"build={summary.get('BuildStatus')}, "
            f"match={summary.get('MatchingStatus')}, "
            f"T3={summary.get('T3IngredientsPresent')}, "
            f"Weinberg={summary.get('WeinbergOperatorPresent')}"
        )

        # If matching did not generate the Weinberg operator, there is
        # no C5 coefficient to extract.
        if not summary.get("WeinbergOperatorPresent"):
            continue

        extraction = summary.get(
            "WeinbergExtractionStatus",
            "Unknown",
        )

        if extraction == "Success":
            n_holo = summary.get(
                "WeinbergHolomorphicTermCount",
                0,
            )
            n_hc = summary.get(
                "WeinbergConjugateTermCount",
                0,
            )

            print(
                f"  C5 extraction: Success "
                f"({n_holo} holomorphic + {n_hc} HC terms)"
            )

            if coefficient_file := summary.get("WeinbergCoefficientFile"):
                print(
                    f"  C5: "
                    f"{record.output_dir / coefficient_file}"
                )

        else:
            print(
                f"  Weinberg terms: "
                f"{summary.get('WeinbergTermCount', 0)}; "
                f"C5: {extraction}"
            )

    # Save all model summaries together so we can easily compare runs
    # or use them for regression testing.
    aggregate = OUTPUT_DIR / "t3_model_comparison.json"

    aggregate.write_text(
        json.dumps(
            [record.summary for record in records],
            indent=2,
        ),
        encoding="utf-8",
    )

    print(
        f"\n{successful}/{len(records)} completed build+matching."
        f"\nAggregate: {aggregate}"
    )

    # Standard command-line convention:
    #   0 = everything succeeded
    #   1 = at least one model failed
    return 0 if successful == len(records) else 1


def organise_c5_input(record: RunRecord) -> Path | None:
    """Move the machine-readable matched coefficient into the data folder."""

    summary = record.summary
    coefficient_file = summary.get("WeinbergCoefficientFile")

    if not coefficient_file:
        return None

    c5_path = record.output_dir / coefficient_file

    if not c5_path.is_file():
        return None

    data_dir = record.output_dir / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    organised_path = data_dir / c5_path.name

    if c5_path != organised_path:
        if organised_path.exists():
            raise FileExistsError(
                f"Refusing to overwrite existing RGE input: {organised_path}"
            )

        c5_path.replace(organised_path)
        c5_path = organised_path

    summary["WeinbergCoefficientFile"] = (
        c5_path.relative_to(record.output_dir).as_posix()
    )
    return c5_path


def finish_runs(
    records: list[RunRecord],
    numerical_config: Path | None = None,
    debug_reports: bool = False,
) -> int:
    """Print the scan summary and generate all Lagrangian reports."""

    # Organise the matched coefficient before printing paths or writing the
    # aggregate summary, so every reported filename points to its final place.
    for record in records:
        organise_c5_input(record)

    status = print_summary(records)

    write_reports(records, debug_reports)

    for record in records:
        coefficient_pdf = record.output_dir / "c5_coefficient.pdf"
        record.summary["WeinbergCoefficientPDFFile"] = (
            coefficient_pdf.name if coefficient_pdf.exists() else ""
        )

    for record in records:
        summary = record.summary

        if summary.get("BuildStatus") != "Success":
            continue

        if summary.get("MatchingStatus") != "Success":
            continue

        if summary.get("WeinbergExtractionStatus") != "Success":
            continue

        coefficient_file = summary.get("WeinbergCoefficientFile")

        if not coefficient_file:
            continue

        c5_path = record.output_dir / coefficient_file

        if not c5_path.exists():
            continue

        print(f"  {record.name}: starting matched-EFT RGE stage...", flush=True)
        try:
            rge_summary = run_matched_eft_rge(
                c5_path=c5_path,
                output_dir=record.output_dir,
                debug_outputs=debug_reports,
            )
        except Exception as exc:
            summary["RGEStatus"] = "Failed"
            summary["RGEError"] = str(exc)
            status = 1
            print(
                f"  {record.name}: matched-EFT RGE failed: {exc}"
            )
            write_and_compile_rge_report(record)
            continue

        summary.update(rge_summary)

        print(f"  {record.name}: starting symbolic full-flavor RGE stage...", flush=True)
        try:
            flavor_summary = run_flavor_matched_rge(
                c5_path=c5_path,
                output_dir=record.output_dir,
                debug_outputs=debug_reports,
            )
        except Exception as exc:
            summary["FlavorRGEStatus"] = "Failed"
            summary["FlavorRGEError"] = str(exc)
            status = 1
            print(
                f"  {record.name}: full-flavor RGE failed: {exc}"
            )
            write_and_compile_rge_report(record)
            continue

        summary.update(flavor_summary)

        print(f"  {record.name}: starting symbolic neutrino mass stage...", flush=True)
        try:
            mass_summary = run_neutrino_mass_stage(
                c5_path=c5_path,
                output_dir=record.output_dir,
            )
        except Exception as exc:
            summary["NeutrinoMassStatus"] = "Failed"
            summary["NeutrinoMassError"] = str(exc)
            status = 1
            print(
                f"  {record.name}: neutrino mass stage failed: {exc}"
            )
            write_and_compile_rge_report(record)
            continue

        summary.update(mass_summary)

        if numerical_config is not None:
            print(f"  {record.name}: starting numerical RGE stage...", flush=True)
            try:
                numerical_summary = run_numerical_pipeline_stage(
                    c5_path=c5_path,
                    output_dir=record.output_dir,
                    config_path=numerical_config,
                )
            except Exception as exc:
                summary["NumericalRGEStatus"] = "Failed"
                summary["NumericalRGEError"] = str(exc)
                status = 1
                print(
                    f"  {record.name}: numerical RGE failed: {exc}"
                )
                write_and_compile_rge_report(record)
                continue

            summary.update(numerical_summary)

            print(f"  {record.name}: numerical RGE calculation finished.", flush=True)

            mass_matrix_path = (
                record.output_dir
                / numerical_summary["NeutrinoMassMatrixLowScaleFile"]
            )
            numerical_payload = json.loads(
                numerical_config.read_text(encoding="utf-8")
            )
            ordering = numerical_payload.get("ordering", "NO")
            print(f"  {record.name}: starting neutrino observables...", flush=True)
            try:
                observable_summary = run_neutrino_observables_stage(
                    mass_matrix_path=mass_matrix_path,
                    output_dir=record.output_dir,
                     ordering=ordering,
                )
            except Exception as exc:
                summary["NeutrinoObservableStatus"] = "Failed"
                summary["NeutrinoObservableError"] = str(exc)
                status = 1
                print(
                    f"  {record.name}: neutrino observables failed: {exc}"
                )
                write_and_compile_rge_report(record)
                continue

            summary.update(observable_summary)

            print(
                f"  {record.name}: neutrino observables=Success"
                f" -> "
                f"{record.output_dir / observable_summary['NeutrinoObservablesFile']}"
            )

            print(
                f"  {record.name}: numerical RGE=Success"
                f" -> "
                f"{record.output_dir / numerical_summary['NeutrinoMassMatrixLowScaleFile']}"
            )

        print(
            f"  {record.name}: neutrino mass=Success"
            f" -> "
            f"{record.output_dir / mass_summary['NeutrinoMassMatrixFile']}"
        )

        print(
            f"  {record.name}: full-flavor RGE=Success"
            f" -> "
            f"{record.output_dir / flavor_summary['C5FlavorBetaMatrixFile']}"
        )

        print(
            f"  {record.name}: matched-EFT RGE=Success"
            f" -> {record.output_dir / rge_summary['C5BetaFile']}"
        )

        write_and_compile_rge_report(record)

    aggregate = OUTPUT_DIR / "t3_model_comparison.json"

    aggregate.write_text(
        json.dumps(
            [record.summary for record in records],
            indent=2,
        ),
        encoding="utf-8",
    )

    return status


def write_and_compile_rge_report(record: RunRecord) -> Path:
    """Create and compile the concise per-model RGE report."""

    report_tex = write_rge_report(
        output_dir=record.output_dir,
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


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Run T3 matching for known benchmark models "
            "or arbitrary valid SU(2) irreps."
        )
    )

    mode = parser.add_mutually_exclusive_group()

    # Simple quick run, we first fill in the details based on input arguments
    mode.add_argument(
        "--smoke",
        action="store_true",
        help="five T3 regression models",
    )

    # More run
    mode.add_argument(
        "--extended",
        action="store_true",
        help="seven historical benchmark points",
    )

    # Dimension input
    mode.add_argument(
        "--dims",
        nargs=3,
        type=int,
        metavar=("DS1", "DS2", "DF"),
        help=(
            "run one representation assignment, "
            "e.g. --dims 3 5 4"
        ),
    )

    # Gives alpha
    parser.add_argument(
        "--numerical",
        type=Path,
        default=None,
        help=(
            "optional JSON parameter point for numerical "
            "matched-EFT running"
        ),
    )

    parser.add_argument(
        "--alpha",
        type=int,
        default=0,
        help=(
            "T3 hypercharge parameter for --dims mode "
            "(default: 0)"
        ),
    )

    parser.add_argument(
        "--debug-reports",
        action="store_true",
        help=(
            "also keep raw Wolfram logs and generate the full "
            "UV/EFT expression reports"
        ),
    )

    args = parser.parse_args()

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )


    # If the dimensions happen to match T3-A ... T3-E, run_dimensions()
    # automatically recognises and labels the model appropriately.
    if args.dims:
        d_s1, d_s2, d_f = args.dims

        try:
            record = run_dimensions(
                d_s1,
                d_s2,
                d_f,
                args.alpha,
                args.debug_reports,
            )
        except ValueError as exc:
            parser.error(str(exc))

        # This gives us  our reports based on our records
        return finish_runs(
            [record],
            args.numerical,
            args.debug_reports,
        )

    # The benchmark lists still use the familiar A-E notation because
    # it is convenient for regression testing and comparison with the paper.
    if args.smoke:
        mode_name = "smoke"
        points = SMOKE

    elif args.extended:
        mode_name = "extended"
        points = EXTENDED

    else:
        mode_name = "interesting"
        points = INTERESTING

    print(
        f"T3 scan mode: {mode_name}; "
        f"{len(points)} model(s)."
    )

    # This line only exist if we use the special interesting, extended and smoke options where we dont input any dimensions
    # Known A-E models are converted to dimensions first and then sent
    # through exactly the same run_dimensions() path as generalised models.
    records = [
        run_known_class(
            model_class,
            alpha,
            args.debug_reports,
        )
        for model_class, alpha in points
    ]

    # This only gives us our reports created from the records
    return finish_runs(
        records,
        args.numerical,
        args.debug_reports,
    )


if __name__ == "__main__":
    raise SystemExit(main())
