from __future__ import annotations

import re

import sympy as sp

"""Across-model comparison report for the UV one-loop T3 RGEs from RGBeta."""

import json
from pathlib import Path
from typing import Any

from common.Paths import REPORT_OUTPUT_DIR
from common.Records import RunRecord
from Reports.StageReports import final_eft_stage_label, rge_report_path
from Reports.ReportGeneration import (
    compile_latex_document,
    latex_escape_text,
    latex_fraction,
    matrix_cell,
    record_quantum_numbers,
    split_latex_terms,
)


# Put the common couplings first; representation-dependent quartics follow in
# alphabetical order so new RGBeta couplings are picked up automatically.
PREFERRED_COUPLING_ORDER = (
    "gY",
    "g2",
    "g3",
    "yu",
    "yd",
    "ye",
    "y1",
    "y2",
    "MF",
    "mS1Sq",
    "mS2Sq",
    "lambdaH",
    "lambdaS1",
    "lambdaS2",
    "lambdaH1",
    "lambdaH2",
    "lambda12",
    "lambdaT3",
)




GLOSSARY_ROWS = (
    (r"$g_Y$", r"$U(1)_Y$ gauge coupling", r"Contained in $D_\mu=\partial_\mu+i g_Y Y B_\mu+\cdots$.", r"Hypercharge gauge interaction."),
    (r"$g_2$", r"$SU(2)_L$ gauge coupling", r"Contained in $D_\mu=\partial_\mu+i g_2 T^a W^a_\mu+\cdots$.", r"Strongly sensitive to the $SU(2)_L$ representations of the new fields."),
    (r"$g_3$", r"$SU(3)_c$ gauge coupling", r"Contained in the QCD covariant derivative.", r"The T3 fields used here are colour singlets, so their one-loop contribution to $\beta_{g_3}$ vanishes."),
    (r"$Y_u$", r"SM up-type Yukawa matrix", r"$-\bar Q\,Y_u\,\widetilde H\,u+\mathrm{h.c.}$", r"SM Yukawa interaction."),
    (r"$Y_d$", r"SM down-type Yukawa matrix", r"$-\bar Q\,Y_d\,H\,d+\mathrm{h.c.}$", r"SM Yukawa interaction."),
    (r"$Y_e$", r"SM charged-lepton Yukawa matrix", r"$-\bar L\,Y_e\,H\,e+\mathrm{h.c.}$", r"SM charged-lepton Yukawa interaction."),
    (r"$y_1$", r"T3 Yukawa matrix", r"$-\,\bar L\,y_1\,F\,S_1^{(\dagger)}+\mathrm{h.c.}$", r"Schematic only: the required conjugation and $SU(2)$ invariant tensor depend on the T3 representation assignment."),
    (r"$y_2$", r"T3 Yukawa matrix", r"$-\,\bar L\,y_2\,F\,S_2^{(\dagger)}+\mathrm{h.c.}$", r"Schematic only. Together with $y_1$ and $\lambda_{T3}$ it participates in the neutrino-mass loop."),
    (r"$M_F$", r"BSM fermion mass matrix", r"$-\bar F\,M_F F$ (Dirac/vector-like case), or the corresponding allowed mass term.", r"May be matrix-valued when several generations of $F$ are present."),
    (r"$m_{S_1}^2$", r"scalar mass-squared", r"$-m_{S_1}^2 S_1^\dagger S_1$", r"Quadratic scalar-potential parameter."),
    (r"$m_{S_2}^2$", r"scalar mass-squared", r"$-m_{S_2}^2 S_2^\dagger S_2$", r"Quadratic scalar-potential parameter."),
    (r"$\lambda_H$", r"SM Higgs quartic", r"$-\lambda_H(H^\dagger H)^2$", r"Overall sign follows the project's Lagrangian/potential convention."),
    (r"$\lambda_{S_1}$", r"$S_1$ self quartic", r"$-\lambda_{S_1}(S_1^\dagger S_1)^2$", r"Singlet contraction of the scalar bilinears."),
    (r"$\lambda_{S_2}$", r"$S_2$ self quartic", r"$-\lambda_{S_2}(S_2^\dagger S_2)^2$", r"Singlet contraction of the scalar bilinears."),
    (r"$\lambda_{H1}$", r"Higgs--$S_1$ quartic", r"$-\lambda_{H1}(H^\dagger H)(S_1^\dagger S_1)$", r"Singlet--singlet mixed quartic."),
    (r"$\lambda_{H2}$", r"Higgs--$S_2$ quartic", r"$-\lambda_{H2}(H^\dagger H)(S_2^\dagger S_2)$", r"Singlet--singlet mixed quartic."),
    (r"$\lambda_{12}$", r"$S_1$--$S_2$ quartic", r"$-\lambda_{12}(S_1^\dagger S_1)(S_2^\dagger S_2)$", r"Singlet--singlet mixed quartic."),
    (r"$\lambda_{T3}$", r"T3 loop-closing quartic", r"$-\lambda_{T3}\,H H S_1 S_2+\mathrm{h.c.}$", r"Schematic field content. Conjugations and the invariant tensor are fixed by hypercharge and the $SU(2)$ representations; this is the quartic joining the two scalar lines to the two Higgs legs in T3."),
    (r"$\lambda_{H1}^{\rm Adj}$", r"adjoint Higgs--$S_1$ quartic", r"$-\lambda_{H1}^{\rm Adj}(H^\dagger T^a H)(S_1^\dagger T^a S_1)$", r"Exists only when the relevant representations admit the adjoint bilinear."),
    (r"$\lambda_{H2}^{\rm Adj}$", r"adjoint Higgs--$S_2$ quartic", r"$-\lambda_{H2}^{\rm Adj}(H^\dagger T^a H)(S_2^\dagger T^a S_2)$", r"Independent from the singlet contraction $\lambda_{H2}$."),
    (r"$\lambda_{12}^{\rm Adj}$", r"adjoint $S_1$--$S_2$ quartic", r"$-\lambda_{12}^{\rm Adj}(S_1^\dagger T^a S_1)(S_2^\dagger T^a S_2)$", r"Representation-dependent independent invariant."),
    (r"$\lambda_{S_1}^{\rm Adj}$", r"adjoint-type $S_1$ self quartic", r"$-\lambda_{S_1}^{\rm Adj}(S_1^\dagger T^a S_1)(S_1^\dagger T^a S_1)$", r"Only appears where that self-contraction is independent/non-vanishing."),
    (r"$\lambda_{S_2}^{\rm Adj}$", r"adjoint-type $S_2$ self quartic", r"$-\lambda_{S_2}^{\rm Adj}(S_2^\dagger T^a S_2)(S_2^\dagger T^a S_2)$", r"Only appears where that self-contraction is independent/non-vanishing."),
    (r"$\lambda_{12}^{\rm Cross}$", r"additional crossed $S_1$--$S_2$ invariant", r"Model-dependent four-scalar contraction of $S_1,S_1^\dagger,S_2,S_2^\dagger$.", r"The exact Clebsch--Gordan/invariant-tensor contraction is representation dependent, so the report does not guess a universal component formula."),
    (r"$\lambda_{H^\dagger H^\dagger S_1^\dagger S_1^\dagger}$", r"special quartic", r"$-\lambda_{H^\dagger H^\dagger S_1^\dagger S_1^\dagger}\,H^\dagger H^\dagger S_1^\dagger S_1^\dagger+\mathrm{h.c.}$", r"Only registered when gauge invariance allows this field content."),
    (r"$\lambda_{H^\dagger H^\dagger S_2 S_2}$", r"special quartic", r"$-\lambda_{H^\dagger H^\dagger S_2 S_2}\,H^\dagger H^\dagger S_2 S_2+\mathrm{h.c.}$", r"Representation/hypercharge dependent."),
    (r"$\lambda_{S_1^{\dagger 2}S_2^{\dagger 2}}$", r"special $S_1$--$S_2$ quartic", r"$-\lambda_{S_1^{\dagger 2}S_2^{\dagger 2}}(S_1^\dagger)^2(S_2^\dagger)^2+\mathrm{h.c.}$", r"Exact invariant contraction is model dependent."),
    (r"$\lambda_{S_1^\dagger S_2 S_2^{\dagger 2}}$", r"special $S_1$--$S_2$ quartic", r"$-\lambda_{S_1^\dagger S_2 S_2^{\dagger 2}}\,S_1^\dagger S_2(S_2^\dagger)^2+\mathrm{h.c.}$", r"Field-content notation mirrors the RGBeta coupling name."),
    (r"$\lambda_{S_1 S_1^{\dagger 2}S_2^\dagger}$", r"special $S_1$--$S_2$ quartic", r"$-\lambda_{S_1 S_1^{\dagger 2}S_2^\dagger}\,S_1(S_1^\dagger)^2S_2^\dagger+\mathrm{h.c.}$", r"Field-content notation mirrors the RGBeta coupling name."),
    (r"$\lambda_{H H^\dagger S_1^\dagger S_2^\dagger}^{\rm Cross}$", r"special crossed Higgs--scalar quartic", r"A crossed invariant with field content $H,H^\dagger,S_1^\dagger,S_2^\dagger$.", r"Exact $SU(2)$ tensor contraction depends on the representation assignment."),
)


COUPLING_LATEX = {
    "gY": r"g_Y",
    "g2": r"g_2",
    "g3": r"g_3",
    "yu": r"Y_u",
    "yd": r"Y_d",
    "ye": r"Y_e",
    "y1": r"y_1",
    "y2": r"y_2",
    "MF": r"M_F",
    "mS1Sq": r"m_{S_1}^2",
    "mS2Sq": r"m_{S_2}^2",
    "lambdaH": r"\lambda_H",
    "lambdaS1": r"\lambda_{S_1}",
    "lambdaS2": r"\lambda_{S_2}",
    "lambdaH1": r"\lambda_{H1}",
    "lambdaH2": r"\lambda_{H2}",
    "lambda12": r"\lambda_{12}",
    "lambdaT3": r"\lambda_{T3}",
    "lambdaH1Adj": r"\lambda_{H1}^{\rm Adj}",
    "lambdaH2Adj": r"\lambda_{H2}^{\rm Adj}",
    "lambdaS1Adj": r"\lambda_{S_1}^{\rm Adj}",
    "lambdaS2Adj": r"\lambda_{S_2}^{\rm Adj}",
    "lambda12Adj": r"\lambda_{12}^{\rm Adj}",
    "lambda12Cross": r"\lambda_{12}^{\rm Cross}",
    "lambdaHHdagS2S2": r"\lambda_{H^\dagger H^\dagger S_2 S_2}",
    "lambdaHHdagS1barS1bar": r"\lambda_{H^\dagger H^\dagger S_1^\dagger S_1^\dagger}",
    "lambdaS1bar2S2bar2": r"\lambda_{S_1^{\dagger 2}S_2^{\dagger 2}}",
    "lambdaS1barS2S2bar2": r"\lambda_{S_1^\dagger S_2 S_2^{\dagger 2}}",
    "lambdaS1S1bar2S2bar": r"\lambda_{S_1 S_1^{\dagger 2}S_2^\dagger}",
    "lambdaHHdagS1barS2barCross": (
        r"\lambda_{H H^\dagger S_1^\dagger S_2^\dagger}^{\rm Cross}"
    ),
}


def _load_uv_payload(record: RunRecord) -> dict[str, Any] | None:
    """Read one model's saved RGBeta JSON payload."""

    relative_path = record.summary.get("UVRGEFile")
    if not relative_path:
        return None

    path = record.output_dir / str(relative_path)
    if not path.is_file():
        return None

    try:
        payload = json.loads(path.read_text(encoding="utf-8-sig"))
    except (OSError, json.JSONDecodeError):
        return None

    if payload.get("status") != "Success":
        return None

    return payload



def _load_json_payload_from_summary(
    record: RunRecord,
    summary_key: str,
    *,
    require_success_key: str | None = None,
) -> dict[str, Any] | None:
    """Load one JSON payload referenced by a record summary field."""

    if require_success_key is not None:
        if record.summary.get(require_success_key) != "Success":
            return None

    relative_path = record.summary.get(summary_key)
    if not relative_path:
        return None

    path = record.output_dir / str(relative_path)
    if not path.is_file():
        return None

    try:
        payload = json.loads(path.read_text(encoding="utf-8-sig"))
    except (OSError, json.JSONDecodeError):
        return None

    if payload.get("status") != "Success":
        return None

    return payload


def _load_eft1_renormalisable_payload(
    record: RunRecord,
) -> dict[str, Any] | None:
    """Read the saved RGBeta payload for EFT1 = SM + S1 + S2."""

    return _load_json_payload_from_summary(
        record,
        "EFT1RenormalisableRGEFile",
        require_success_key="EFT1RenormalisableRGEStatus",
    )


def _load_eft1_wilson_payload(
    record: RunRecord,
) -> dict[str, Any] | None:
    """Read the component-level dimension-five Wilson RGE payload."""

    return _load_json_payload_from_summary(
        record,
        "EFT1WilsonRGEFile",
        require_success_key="EFT1WilsonRGEStatus",
    )


def _coupling_symbol(name: str) -> str:
    return COUPLING_LATEX.get(name, r"\mathrm{" + latex_escape_text(name) + "}")


def _ordered_couplings(names: set[str]) -> list[str]:
    preferred = [name for name in PREFERRED_COUPLING_ORDER if name in names]
    remaining = sorted(names.difference(preferred))
    return preferred + remaining


def _extract_balanced(text: str, open_index: int, opener: str, closer: str) -> tuple[str, int] | None:
    """Return the contents and closing index of one balanced delimiter group."""

    if open_index >= len(text) or text[open_index] != opener:
        return None

    depth = 0
    for index in range(open_index, len(text)):
        char = text[index]
        if char == opener:
            depth += 1
        elif char == closer:
            depth -= 1
            if depth == 0:
                return text[open_index + 1:index], index
    return None


def _split_top_level_commas(text: str) -> list[str]:
    """Split function arguments without splitting nested function calls."""

    parts: list[str] = []
    start = 0
    paren = bracket = brace = 0
    for index, char in enumerate(text):
        if char == "(":
            paren += 1
        elif char == ")":
            paren -= 1
        elif char == "[":
            bracket += 1
        elif char == "]":
            bracket -= 1
        elif char == "{":
            brace += 1
        elif char == "}":
            brace -= 1
        elif char == "," and paren == bracket == brace == 0:
            parts.append(text[start:index].strip())
            start = index + 1
    parts.append(text[start:].strip())
    return [part for part in parts if part]


def _clean_index(text: str) -> str:
    """Turn RGBeta flavor-index wrappers into ordinary i,j,... indices."""

    value = text.strip()
    value = value.replace(r"\text{$\$$i}", "i")
    value = value.replace(r"\text{$\$$j}", "j")
    value = value.replace(r"\text{$\$$k}", "k")
    value = value.replace(r"\text{$\$$l}", "l")
    value = value.replace(r"\text{$\$$m}", "m")
    value = value.replace(r"\text{$\$$n}", "n")
    return value


def _rewrite_rgbeta_calls(latex: str) -> str:
    """Convert RGBeta's internal TeXForm heads to conventional matrix notation."""

    text = latex

    # Work from the innermost call outward.  TeXForm renders the relevant
    # RGBeta heads as \text{Name}(...) (or Tr[...]).
    call_names = ("Bar", "Trans", "gen", "heavy", "Matrix")

    changed = True
    while changed:
        changed = False
        best_pos = -1
        best_name = None
        for name in call_names:
            pos = text.rfind(r"\text{" + name + "}(")
            if pos > best_pos:
                best_pos = pos
                best_name = name

        if best_name is not None and best_pos >= 0:
            prefix = r"\text{" + best_name + "}"
            open_index = best_pos + len(prefix)
            parsed = _extract_balanced(text, open_index, "(", ")")
            if parsed is not None:
                argument, close_index = parsed
                end_index = close_index

                if best_name in {"gen", "heavy"}:
                    replacement = _clean_index(argument)

                elif best_name == "Bar":
                    replacement = rf"\left({argument}\right)^{{*}}"

                elif best_name == "Trans":
                    stripped = argument.strip()
                    star_suffix = r"\right)^{*}"
                    if stripped.startswith(r"\left(") and stripped.endswith(star_suffix):
                        core = stripped[len(r"\left("):-len(star_suffix)]
                        replacement = rf"\left({core}\right)^{{\dagger}}"
                    else:
                        replacement = rf"\left({argument}\right)^{{T}}"

                else:  # Matrix
                    factors = _split_top_level_commas(argument)
                    product = r"\,".join(factors)
                    replacement = product

                    # RGBeta represents an indexed matrix product as
                    # Matrix[A,B,...](gen[$i],gen[$j]).
                    next_index = close_index + 1
                    if next_index < len(text) and text[next_index] == "(":
                        parsed_indices = _extract_balanced(text, next_index, "(", ")")
                        if parsed_indices is not None:
                            indices_text, second_close = parsed_indices
                            indices = _split_top_level_commas(indices_text)
                            indices = [_clean_index(index) for index in indices]
                            if indices:
                                replacement = (
                                    rf"\left({product}\right)_{{"
                                    + "".join(indices)
                                    + "}"
                                )
                            end_index = second_close

                text = text[:best_pos] + replacement + text[end_index + 1:]
                changed = True
                continue

        # Tr is rendered with square brackets.
        tr_pos = text.rfind(r"\text{Tr}[")
        if tr_pos >= 0:
            open_index = tr_pos + len(r"\text{Tr}")
            parsed = _extract_balanced(text, open_index, "[", "]")
            if parsed is not None:
                argument, close_index = parsed
                argument = argument.replace(".", r"\,")
                replacement = rf"\operatorname{{Tr}}\!\left({argument}\right)"
                text = text[:tr_pos] + replacement + text[close_index + 1:]
                changed = True

    return text


def _normalise_rgbeta_latex(latex: str) -> str:
    """Return human-readable physics LaTeX from RGBeta's TeXForm output."""

    text = _rewrite_rgbeta_calls(latex.strip())

    # RGBeta symbols are ordinary Mathematica symbols, so TeXForm otherwise
    # emits them as \text{symbol}.  Map the physics names explicitly.
    symbol_map = {
        "gY": r"g_Y",
        "g2": r"g_2",
        "g3": r"g_3",
        "yu": r"Y_u",
        "yd": r"Y_d",
        "ye": r"Y_e",
        "y1": r"y_1",
        "y2": r"y_2",
        "MF": r"M_F",
        "mS1Sq": r"m_{S_1}^2",
        "mS2Sq": r"m_{S_2}^2",
        "lambdaH": r"\lambda_H",
        "lambdaS1": r"\lambda_{S_1}",
        "lambdaS2": r"\lambda_{S_2}",
        "lambdaH1": r"\lambda_{H1}",
        "lambdaH2": r"\lambda_{H2}",
        "lambda12": r"\lambda_{12}",
        "lambdaT3": r"\lambda_{T3}",
        "lambdaH1Adj": r"\lambda_{H1}^{\rm Adj}",
        "lambdaH2Adj": r"\lambda_{H2}^{\rm Adj}",
        "lambdaS1Adj": r"\lambda_{S_1}^{\rm Adj}",
        "lambdaS2Adj": r"\lambda_{S_2}^{\rm Adj}",
        "lambda12Adj": r"\lambda_{12}^{\rm Adj}",
        "lambda12Cross": r"\lambda_{12}^{\rm Cross}",
        "lambdaHHdagS2S2": r"\lambda_{H^\dagger H^\dagger S_2 S_2}",
        "lambdaHHdagS1barS1bar": r"\lambda_{H^\dagger H^\dagger S_1^\dagger S_1^\dagger}",
        "lambdaS1bar2S2bar2": r"\lambda_{S_1^{\dagger 2}S_2^{\dagger 2}}",
        "lambdaS1barS2S2bar2": r"\lambda_{S_1^\dagger S_2 S_2^{\dagger 2}}",
        "lambdaS1S1bar2S2bar": r"\lambda_{S_1 S_1^{\dagger 2}S_2^\dagger}",
        "lambdaHHdagS1barS2barCross": r"\lambda_{H H^\dagger S_1^\dagger S_2^\dagger}^{\rm Cross}",
    }

    # Long names must be replaced before shorter names such as lambdaH.
    for name in sorted(symbol_map, key=len, reverse=True):
        # Treat every replacement as one TeX atom.  This is important for
        # couplings that already carry a superscript in their display name,
        # e.g. lambdaH1Adj -> \lambda_{H1}^{\rm Adj}.  RGBeta may then
        # raise that coupling to a power; without grouping, TeX sees
        # \lambda_{H1}^{\rm Adj}^2 and aborts with "Double superscript".
        replacement = "{" + symbol_map[name] + "}"
        text = text.replace(r"\text{" + name + "}", replacement)

    text = _clean_index(text)
    return text


def _beta_cell(latex: str | None, raw: str | None) -> str:
    """Render one beta function as additive lines inside a comparison cell."""

    if latex and latex.strip():
        cleaned = _normalise_rgbeta_latex(latex)
        terms = split_latex_terms(cleaned)
        return matrix_cell(terms or [cleaned], empty_value=r"---")

    if raw and raw.strip():
        # Old output files can still be inspected, but rerunning is required for
        # the physics-aware LaTeX conversion and conventional gauge normalization.
        return r"\scriptsize\texttt{" + latex_escape_text(raw.strip()) + r"}"

    return r"---"


def _comparison_term_signature(term: str) -> str:
    """Return a coefficient-insensitive key used to align RGE terms across models.

    The *cell* always keeps the complete original term.  This key is used only
    to decide which terms belong in the same comparison column.
    """

    text = term.strip()
    text = re.sub(r"^[+-]\s*", "", text)
    text = text.replace(r"\left", "").replace(r"\right", "")
    text = re.sub(r"\s+", "", text)

    # Ignore representation-dependent numerical prefactors while retaining the
    # symbolic tensor/matrix structure.  Powers/subscripts such as g_2^2 remain.
    text = re.sub(
        r"\\frac\{[-+]?\d+\}\{[-+]?\d+\}",
        r"\\mathsf{c}",
        text,
    )
    text = re.sub(
        r"\\sqrt\{\d+\}",
        r"\\mathsf{c}",
        text,
    )

    # Replace ordinary multiplicative integers, but avoid numbers that are
    # explicitly part of a subscript or superscript.
    text = re.sub(
        r"(?<![_^])(?<![_^]\{)(?<![A-Za-z])\d+(?![A-Za-z])",
        "c",
        text,
    )

    return text


def _coupling_term_rows(
    rows: list[tuple[RunRecord, dict[str, Any]]],
    coupling: str,
) -> tuple[
    list[str],
    dict[int, dict[str, list[str]]],
    dict[str, str],
]:
    """Collect and align the complete additive terms of one beta function."""

    ordered_signatures: list[str] = []
    model_terms: dict[int, dict[str, list[str]]] = {}
    representative: dict[str, str] = {}

    for row_index, (_record, payload) in enumerate(rows):
        latex_betas = payload.get("report_beta_latex", {}) or {}
        latex = latex_betas.get(coupling)

        if not latex or not str(latex).strip():
            model_terms[row_index] = {}
            continue

        cleaned = _normalise_rgbeta_latex(str(latex))
        terms = split_latex_terms(cleaned) or [cleaned]

        grouped: dict[str, list[str]] = {}

        for term in terms:
            signature = _comparison_term_signature(term)

            if signature not in representative:
                representative[signature] = term

            if signature not in ordered_signatures:
                ordered_signatures.append(signature)

            grouped.setdefault(signature, []).append(term)

        model_terms[row_index] = grouped

    return ordered_signatures, model_terms, representative


def _rge_term_cell(terms: list[str] | None) -> str:
    """Render the complete model-specific RGE term(s) in one comparison cell."""

    if not terms:
        return r"---"

    return matrix_cell(terms, empty_value=r"---")

def write_rge_comparison(records: list[RunRecord]) -> Path:
    """Create the across-model UV RGE comparison report.

    There is one section per running coupling.  Within that section, rows are
    model configurations and columns are aligned additive RGE terms.  Cells
    contain the complete term, including its model-dependent coefficient.
    """

    output_path = rge_report_path("UV")

    rows: list[tuple[RunRecord, dict[str, Any]]] = []
    coupling_names: set[str] = set()

    for record in records:
        if record.summary.get("UVRGEStatus") != "Success":
            continue

        payload = _load_uv_payload(record)
        if payload is None:
            continue

        raw_betas = payload.get("betas", {})
        latex_betas = payload.get("report_beta_latex", {})
        report_betas = payload.get("report_betas", {})

        if not isinstance(raw_betas, dict):
            raw_betas = {}
        if not isinstance(latex_betas, dict):
            latex_betas = {}
        if not isinstance(report_betas, dict):
            report_betas = {}

        coupling_names.update(str(name) for name in report_betas or raw_betas)
        coupling_names.update(str(name) for name in latex_betas)
        rows.append((record, payload))

    lines: list[str] = [
        r"\documentclass[8pt]{article}",
        r"\usepackage[margin=0.65cm]{geometry}",
        r"\usepackage{amsmath,amssymb,adjustbox,pdflscape,longtable,array,booktabs}",
        r"\usepackage[T1]{fontenc}",
        r"\setlength{\tabcolsep}{2pt}",
        r"\renewcommand{\arraystretch}{1.2}",
        r"\begin{document}",
        r"\begin{landscape}",
        r"\section*{T3 UV one-loop RGE term comparison}",
        (
            r"Each running coupling has its own table. Rows are model configurations "
            r"and columns are additive structures in the beta function. Every populated "
            r"cell contains the complete model-specific term, including its numerical "
            r"coefficient. A dash means that the term is absent for that configuration."
        ),
        r"\medskip",
        (
            r"All equations use "
            r"$16\pi^2\,\mu\,dX/d\mu=\beta_X^{(1)}$. "
            r"Gauge-coupling results use the same RGBeta $g^2\rightarrow g$ conversion "
            r"as in the previous report."
        ),
    ]

    if not rows:
        lines.append(r"No successful UV RGBeta results are available.")
    else:
        for coupling_number, coupling in enumerate(
            _ordered_couplings(coupling_names),
            start=1,
        ):
            if coupling_number > 1:
                lines.append(r"\clearpage")

            symbol = _coupling_symbol(coupling)

            signatures, model_terms, representative = _coupling_term_rows(
                rows,
                coupling,
            )

            lines.extend(
                [
                    rf"\section*{{$\beta_{{{symbol}}}^{{(1)}}$}}",
                    rf"\[16\pi^2\,\mu\frac{{d {symbol}}}{{d\mu}}="
                    rf"\beta_{{{symbol}}}^{{(1)}}\]",
                ]
            )

            if not signatures:
                lines.append(
                    r"\textit{No LaTeX beta-function terms were available for this coupling.}"
                )
                continue

            # Large RGEs can contain many additive terms. Split the term columns
            # across multiple tables while preserving the same model rows.
            term_chunk_size = 3
            signature_chunks = [
                signatures[index:index + term_chunk_size]
                for index in range(0, len(signatures), term_chunk_size)
            ]

            for chunk_number, signature_chunk in enumerate(
                signature_chunks,
                start=1,
            ):
                if chunk_number > 1:
                    lines.append(r"\clearpage")

                if len(signature_chunks) > 1:
                    lines.append(
                        rf"\subsection*{{Terms {chunk_number} of "
                        rf"{len(signature_chunks)}}}"
                    )

                # Show what each term column represents.  This is only a column
                # label/representative; the table cells below retain each model's
                # complete term.
                lines.append(r"\begin{center}")
                lines.append(
                    r"\begin{tabular}{@{}c >{\raggedright\arraybackslash}p{0.82\linewidth}@{}}"
                )
                lines.append(r"\toprule")
                lines.append(r"Column & Representative term structure \\")
                lines.append(r"\midrule")

                for local_index, signature in enumerate(signature_chunk, start=1):
                    representative_term = representative[signature]
                    lines.append(
                        rf"T{local_index} & "
                        + r"\(\displaystyle "
                        + representative_term
                        + r"\) \\"
                    )
                    lines.append(r"\midrule")

                lines.extend(
                    [
                        r"\bottomrule",
                        r"\end{tabular}",
                        r"\end{center}",
                        r"\smallskip",
                    ]
                )

                widths = " ".join(
                    r">{\raggedright\arraybackslash}p{0.21\linewidth}"
                    for _ in signature_chunk
                )

                column_spec = (
                    r"@{}llcccccc "
                    + widths
                    + r"@{}"
                )

                term_headers = [
                    rf"T{index}"
                    for index in range(1, len(signature_chunk) + 1)
                ]

                header = [
                    r"Model",
                    r"$\alpha$",
                    r"$d_{S_1}$",
                    r"$Y_{S_1}$",
                    r"$d_{S_2}$",
                    r"$Y_{S_2}$",
                    r"$d_F$",
                    r"$Y_F$",
                    *term_headers,
                ]

                lines.extend(
                    [
                        r"\tiny",
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

                for row_index, (record, payload) in enumerate(rows):
                    d_s1, y_s1, d_s2, y_s2, d_f, y_f = record_quantum_numbers(record)
                    grouped = model_terms.get(row_index, {})

                    cells = [
                        latex_escape_text(record.name),
                        rf"${record.alpha}$",
                        rf"${d_s1}$",
                        rf"${latex_fraction(y_s1)}$",
                        rf"${d_s2}$",
                        rf"${latex_fraction(y_s2)}$",
                        rf"${d_f}$",
                        rf"${latex_fraction(y_f)}$",
                        *[
                            _rge_term_cell(grouped.get(signature))
                            for signature in signature_chunk
                        ],
                    ]

                    lines.append(" & ".join(cells) + r" \\")
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

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"\nUV RGE term-comparison report:\n{output_path}")
    return output_path

def write_and_compile_rge_comparison(records: list[RunRecord]) -> Path:
    """Write and compile the across-model UV RGE comparison report."""

    report_tex = write_rge_comparison(records)
    compile_latex_document(report_tex)
    return report_tex



def _renormalisable_rows_from_loader(
    records: list[RunRecord],
    loader,
) -> tuple[list[tuple[RunRecord, dict[str, Any]]], set[str]]:
    """Collect successful RGBeta rows and coupling names for one EFT stage."""

    rows: list[tuple[RunRecord, dict[str, Any]]] = []
    coupling_names: set[str] = set()

    for record in records:
        payload = loader(record)
        if payload is None:
            continue

        raw_betas = payload.get("betas", {})
        latex_betas = payload.get("report_beta_latex", {})
        report_betas = payload.get("report_betas", {})

        if not isinstance(raw_betas, dict):
            raw_betas = {}
        if not isinstance(latex_betas, dict):
            latex_betas = {}
        if not isinstance(report_betas, dict):
            report_betas = {}

        coupling_names.update(str(name) for name in report_betas or raw_betas)
        coupling_names.update(str(name) for name in latex_betas)
        rows.append((record, payload))

    return rows, coupling_names


def _sympify_report_expression(raw: str) -> sp.Expr | None:
    """Parse one saved SymPy expression from the Wilson-RGE JSON."""

    try:
        return sp.expand(
            sp.sympify(
                str(raw),
                locals={
                    "conjugate": sp.conjugate,
                    "sqrt": sp.sqrt,
                    "I": sp.I,
                },
            )
        )
    except (TypeError, ValueError, SyntaxError, sp.SympifyError):
        return None


def _wilson_component_sort_key(component: str) -> tuple[int, ...]:
    try:
        return tuple(int(piece) for piece in component.split(","))
    except ValueError:
        return (10**9,)


def _wilson_component_symbol(component: str) -> str:
    pieces = component.split(",")
    if len(pieces) != 4:
        return r"C_{\mathrm{" + latex_escape_text(component) + r"}}"
    i, j, a, b = pieces
    return rf"C_{{{i}{j}{a}{b}}}"


def _eft1_wilson_component_rows(
    records: list[RunRecord],
    component: str,
) -> tuple[
    list[str],
    dict[int, dict[str, list[str]]],
    dict[str, str],
]:
    """Align additive terms of one Wilson beta component across models."""

    signature_order: list[str] = []
    model_terms: dict[int, dict[str, list[str]]] = {}
    representative: dict[str, str] = {}

    for row_index, record in enumerate(records):
        payload = _load_eft1_wilson_payload(record)
        if payload is None:
            model_terms[row_index] = {}
            continue

        beta_entry = (payload.get("betas", {}) or {}).get(component)
        if not isinstance(beta_entry, dict):
            model_terms[row_index] = {}
            continue

        expression = _sympify_report_expression(beta_entry.get("total", ""))
        if expression is None:
            model_terms[row_index] = {}
            continue

        grouped: dict[str, list[str]] = {}
        for term in sp.Add.make_args(sp.expand(expression)):
            signature = _sympy_term_signature(term)
            latex_term = sp.latex(term)
            grouped.setdefault(signature, []).append(latex_term)

            if signature not in representative:
                representative[signature] = sp.latex(
                    sp.expand(term.as_coeff_Mul()[1])
                )
                signature_order.append(signature)

        model_terms[row_index] = grouped

    return signature_order, model_terms, representative


def write_eft1_rge_comparison(records: list[RunRecord]) -> Path:
    """Create the complete intermediate-EFT RGE report after integrating out F.

    The report contains both:
      1. the renormalisable RGBeta running of SM + S1 + S2; and
      2. the one-loop running/mixing of the tree-generated dimension-five
         psi^2 phi^2 Wilson tensor.

    The Wilson section is component-level because the current master-RGE
    calculation is stored in the real scalar basis C_{ijab}.
    """

    stage_label = "EFT_1_after_F"
    output_path = rge_report_path(stage_label)

    ren_rows, coupling_names = _renormalisable_rows_from_loader(
        records,
        _load_eft1_renormalisable_payload,
    )

    wilson_records = [
        record
        for record in records
        if _load_eft1_wilson_payload(record) is not None
    ]

    wilson_components: set[str] = set()
    for record in wilson_records:
        payload = _load_eft1_wilson_payload(record)
        if payload is None:
            continue
        betas = payload.get("betas", {}) or {}
        if isinstance(betas, dict):
            wilson_components.update(str(name) for name in betas)

    lines: list[str] = [
        r"\documentclass[8pt]{article}",
        r"\usepackage[margin=0.65cm]{geometry}",
        r"\usepackage{amsmath,amssymb,adjustbox,pdflscape,longtable,array,booktabs}",
        r"\usepackage[T1]{fontenc}",
        r"\setlength{\tabcolsep}{2pt}",
        r"\renewcommand{\arraystretch}{1.2}",
        r"\begin{document}",
        r"\begin{landscape}",
        r"\section*{T3 intermediate-EFT one-loop RGE term comparison: EFT\_1 after $F$}",
        (
            r"This stage has active field content SM+$S_1+S_2$. "
            r"The report combines the renormalisable RGBeta running with the "
            r"dimension-five $\psi^2\phi^2$ Wilson-coefficient running generated "
            r"after integrating out $F$."
        ),
        r"\medskip",
        (
            r"Throughout, $16\pi^2\,\mu\,dX/d\mu=\beta_X^{(1)}$. "
            r"For the Wilson sector only the tree-generated stage-1 coefficient "
            r"$C^{(0)}$ is inserted into the one-loop anomalous dimension; the "
            r"one-loop threshold piece is retained as a boundary term and is not "
            r"run again at this order."
        ),
        r"\section*{Renormalisable EFT1 couplings}",
    ]

    if not ren_rows:
        lines.append(
            r"\textit{No successful renormalisable EFT1 RGBeta results are available.}"
        )
    else:
        for coupling_number, coupling in enumerate(
            _ordered_couplings(coupling_names),
            start=1,
        ):
            if coupling_number > 1:
                lines.append(r"\clearpage")

            symbol = _coupling_symbol(coupling)
            signatures, model_terms, representative = _coupling_term_rows(
                ren_rows,
                coupling,
            )

            lines.extend(
                [
                    rf"\section*{{$\beta_{{{symbol}}}^{{(1)}}$}}",
                    rf"\[16\pi^2\,\mu\frac{{d {symbol}}}{{d\mu}}="
                    rf"\beta_{{{symbol}}}^{{(1)}}\]",
                ]
            )

            if not signatures:
                lines.append(
                    r"\textit{No LaTeX beta-function terms were available for this coupling.}"
                )
                continue

            term_chunk_size = 3
            signature_chunks = [
                signatures[index:index + term_chunk_size]
                for index in range(0, len(signatures), term_chunk_size)
            ]

            for chunk_number, signature_chunk in enumerate(
                signature_chunks,
                start=1,
            ):
                if chunk_number > 1:
                    lines.append(r"\clearpage")

                if len(signature_chunks) > 1:
                    lines.append(
                        rf"\subsection*{{Terms {chunk_number} of "
                        rf"{len(signature_chunks)}}}"
                    )

                lines.extend(
                    [
                        r"\begin{center}",
                        r"\begin{tabular}{@{}c >{\raggedright\arraybackslash}p{0.82\linewidth}@{}}",
                        r"\toprule",
                        r"Column & Representative term structure \\",
                        r"\midrule",
                    ]
                )
                for local_index, signature in enumerate(signature_chunk, start=1):
                    lines.append(
                        rf"T{local_index} & "
                        + r"\(\displaystyle "
                        + representative[signature]
                        + r"\) \\"
                    )
                    lines.append(r"\midrule")

                lines.extend(
                    [
                        r"\bottomrule",
                        r"\end{tabular}",
                        r"\end{center}",
                        r"\smallskip",
                    ]
                )

                widths = " ".join(
                    r">{\raggedright\arraybackslash}p{0.21\linewidth}"
                    for _ in signature_chunk
                )
                column_spec = r"@{}llcccccc " + widths + r"@{}"
                headers = [
                    r"Model",
                    r"$\alpha$",
                    r"$d_{S_1}$",
                    r"$Y_{S_1}$",
                    r"$d_{S_2}$",
                    r"$Y_{S_2}$",
                    r"$d_F$",
                    r"$Y_F$",
                    *[
                        rf"T{index}"
                        for index in range(1, len(signature_chunk) + 1)
                    ],
                ]

                lines.extend(
                    [
                        r"\tiny",
                        rf"\begin{{longtable}}{{{column_spec}}}",
                        r"\toprule",
                        " & ".join(headers) + r" \\",
                        r"\midrule",
                        r"\endfirsthead",
                        r"\toprule",
                        " & ".join(headers) + r" \\",
                        r"\midrule",
                        r"\endhead",
                    ]
                )

                for row_index, (record, _payload) in enumerate(ren_rows):
                    d_s1, y_s1, d_s2, y_s2, d_f, y_f = record_quantum_numbers(
                        record
                    )
                    grouped = model_terms.get(row_index, {})
                    cells = [
                        latex_escape_text(record.name),
                        rf"${record.alpha}$",
                        rf"${d_s1}$",
                        rf"${latex_fraction(y_s1)}$",
                        rf"${d_s2}$",
                        rf"${latex_fraction(y_s2)}$",
                        rf"${d_f}$",
                        rf"${latex_fraction(y_f)}$",
                        *[
                            _rge_term_cell(grouped.get(signature))
                            for signature in signature_chunk
                        ],
                    ]
                    lines.append(" & ".join(cells) + r" \\")
                    lines.append(r"\midrule")

                lines.extend([r"\bottomrule", r"\end{longtable}"])

    lines.extend(
        [
            r"\clearpage",
            r"\section*{Dimension-five $\psi^2\phi^2$ Wilson coefficients}",
            (
                r"The Wilson tensor is reported in the real-scalar component basis "
                r"$C_{ijab}$. Each component with a nonzero one-loop beta function "
                r"has its own table. Components absent at the matching boundary but "
                r"generated by operator mixing are included automatically."
            ),
        ]
    )

    if not wilson_records:
        lines.append(
            r"\textit{No successful EFT1 Wilson-RGE results are available.}"
        )
    else:
        ordered_components = sorted(
            wilson_components,
            key=_wilson_component_sort_key,
        )

        for component_number, component in enumerate(
            ordered_components,
            start=1,
        ):
            if component_number > 1:
                lines.append(r"\clearpage")

            symbol = _wilson_component_symbol(component)
            signatures, model_terms, representative = _eft1_wilson_component_rows(
                wilson_records,
                component,
            )

            lines.extend(
                [
                    rf"\section*{{$\beta_{{{symbol}}}^{{(1)}}$}}",
                    rf"\[16\pi^2\,\mu\frac{{d {symbol}}}{{d\mu}}="
                    rf"\beta_{{{symbol}}}^{{(1)}}\]",
                ]
            )

            if not signatures:
                lines.append(
                    r"\textit{No nonzero beta-function terms for this component.}"
                )
                continue

            term_chunk_size = 3
            chunks = [
                signatures[index:index + term_chunk_size]
                for index in range(0, len(signatures), term_chunk_size)
            ]

            for chunk_number, signature_chunk in enumerate(chunks, start=1):
                if chunk_number > 1:
                    lines.append(r"\clearpage")

                if len(chunks) > 1:
                    lines.append(
                        rf"\subsection*{{Terms {chunk_number} of {len(chunks)}}}"
                    )

                lines.extend(
                    [
                        r"\begin{center}",
                        r"\begin{tabular}{@{}c >{\raggedright\arraybackslash}p{0.82\linewidth}@{}}",
                        r"\toprule",
                        r"Column & Representative term structure \\",
                        r"\midrule",
                    ]
                )
                for local_index, signature in enumerate(signature_chunk, start=1):
                    lines.append(
                        rf"T{local_index} & "
                        + r"\(\displaystyle "
                        + representative[signature]
                        + r"\) \\"
                    )
                    lines.append(r"\midrule")

                lines.extend(
                    [
                        r"\bottomrule",
                        r"\end{tabular}",
                        r"\end{center}",
                        r"\smallskip",
                    ]
                )

                widths = " ".join(
                    r">{\raggedright\arraybackslash}p{0.21\linewidth}"
                    for _ in signature_chunk
                )
                column_spec = r"@{}llcccccc " + widths + r"@{}"
                headers = [
                    r"Model",
                    r"$\alpha$",
                    r"$d_{S_1}$",
                    r"$Y_{S_1}$",
                    r"$d_{S_2}$",
                    r"$Y_{S_2}$",
                    r"$d_F$",
                    r"$Y_F$",
                    *[
                        rf"T{index}"
                        for index in range(1, len(signature_chunk) + 1)
                    ],
                ]

                lines.extend(
                    [
                        r"\tiny",
                        rf"\begin{{longtable}}{{{column_spec}}}",
                        r"\toprule",
                        " & ".join(headers) + r" \\",
                        r"\midrule",
                        r"\endfirsthead",
                        r"\toprule",
                        " & ".join(headers) + r" \\",
                        r"\midrule",
                        r"\endhead",
                    ]
                )

                for row_index, record in enumerate(wilson_records):
                    d_s1, y_s1, d_s2, y_s2, d_f, y_f = record_quantum_numbers(
                        record
                    )
                    grouped = model_terms.get(row_index, {})
                    cells = [
                        latex_escape_text(record.name),
                        rf"${record.alpha}$",
                        rf"${d_s1}$",
                        rf"${latex_fraction(y_s1)}$",
                        rf"${d_s2}$",
                        rf"${latex_fraction(y_s2)}$",
                        rf"${d_f}$",
                        rf"${latex_fraction(y_f)}$",
                        *[
                            _rge_term_cell(grouped.get(signature))
                            for signature in signature_chunk
                        ],
                    ]
                    lines.append(" & ".join(cells) + r" \\")
                    lines.append(r"\midrule")

                lines.extend([r"\bottomrule", r"\end{longtable}"])

    lines.extend(
        [
            r"\end{landscape}",
            r"\end{document}",
            "",
        ]
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"\nIntermediate-EFT RGE term-comparison report:\n{output_path}")
    return output_path


def write_and_compile_eft1_rge_comparison(
    records: list[RunRecord],
) -> Path:
    """Write and compile the complete EFT1 RGE comparison report."""

    report_tex = write_eft1_rge_comparison(records)
    compile_latex_document(report_tex)
    return report_tex


def _parse_beta_ratio(summary: dict[str, Any]) -> sp.Expr | None:
    """Parse the saved matched-EFT beta/C5 expression."""

    raw = summary.get("BetaOverC5")
    if not raw:
        return None

    try:
        return sp.expand(
            sp.sympify(
                str(raw),
                locals={
                    "conjugate": sp.conjugate,
                },
            )
        )
    except (TypeError, ValueError, SyntaxError, sp.SympifyError):
        return None


def _sympy_term_signature(term: sp.Expr) -> str:
    """Coefficient-insensitive structural key for one additive SymPy term."""

    _, structure = sp.sympify(term).as_coeff_Mul()
    return sp.srepr(structure)


def _final_eft_term_rows(
    records: list[RunRecord],
) -> tuple[
    list[str],
    dict[int, dict[str, list[str]]],
    dict[str, str],
]:
    """Align additive terms of the final SMEFT C5 beta across models."""

    signature_order: list[str] = []
    model_terms: dict[int, dict[str, list[str]]] = {}
    representative: dict[str, str] = {}

    c5 = sp.Symbol("C_5")

    for row_index, record in enumerate(records):
        ratio = _parse_beta_ratio(record.summary)
        if ratio is None:
            continue

        grouped: dict[str, list[str]] = {}

        for term in sp.Add.make_args(sp.expand(ratio)):
            signature = _sympy_term_signature(term)

            # The report is written as beta_C5, not beta_C5/C5, so every
            # populated cell contains the complete additive RGE term.
            latex_term = sp.latex(sp.expand(c5 * term))

            grouped.setdefault(signature, []).append(latex_term)

            if signature not in representative:
                representative[signature] = sp.latex(
                    sp.expand(c5 * term.as_coeff_Mul()[1])
                )
                signature_order.append(signature)

        model_terms[row_index] = grouped

    return signature_order, model_terms, representative


def write_final_eft_rge_comparison(records: list[RunRecord]) -> Path:
    """Create the across-model final-EFT Weinberg RGE comparison report."""

    stage_label = final_eft_stage_label(records)
    output_path = rge_report_path(stage_label)

    successful_records = [
        record
        for record in records
        if record.summary.get("RGEStatus") == "Success"
        and _parse_beta_ratio(record.summary) is not None
    ]

    signatures, model_terms, representative = _final_eft_term_rows(
        successful_records
    )

    lines: list[str] = [
        r"\documentclass[8pt]{article}",
        r"\usepackage[margin=0.65cm]{geometry}",
        r"\usepackage{amsmath,amssymb,adjustbox,pdflscape,longtable,array,booktabs}",
        r"\usepackage[T1]{fontenc}",
        r"\setlength{\tabcolsep}{2pt}",
        r"\renewcommand{\arraystretch}{1.2}",
        r"\begin{document}",
        r"\begin{landscape}",
        rf"\section*{{T3 final-EFT one-loop RGE term comparison: {latex_escape_text(stage_label)}}}",
        (
            r"All T3 particles have been integrated out at this stage. "
            r"The active theory is the SM plus the dimension-five Weinberg operator. "
            r"Rows are the UV model configurations that produced the matched coefficient, "
            r"while columns align the additive structures in its one-loop beta function. "
            r"Each populated cell contains the complete term, including its coefficient "
            r"and the factor $C_5$."
        ),
        r"\medskip",
        (
            r"We use "
            r"$16\pi^2\,\mu\,dC_5/d\mu=\beta_{C_5}^{(1)}$."
        ),
        r"\section*{$\beta_{C_5}^{(1)}$}",
        r"\[16\pi^2\,\mu\frac{dC_5}{d\mu}=\beta_{C_5}^{(1)}\]",
    ]

    if not successful_records:
        lines.append(
            r"\textit{No successful matched-EFT RGE results are available.}"
        )
    elif not signatures:
        lines.append(
            r"\textit{No additive beta-function terms could be parsed.}"
        )
    else:
        term_chunk_size = 3
        chunks = [
            signatures[index:index + term_chunk_size]
            for index in range(0, len(signatures), term_chunk_size)
        ]

        for chunk_number, signature_chunk in enumerate(chunks, start=1):
            if chunk_number > 1:
                lines.append(r"\clearpage")

            if len(chunks) > 1:
                lines.append(
                    rf"\subsection*{{Terms {chunk_number} of {len(chunks)}}}"
                )

            lines.extend(
                [
                    r"\begin{center}",
                    r"\begin{tabular}{@{}c >{\raggedright\arraybackslash}p{0.82\linewidth}@{}}",
                    r"\toprule",
                    r"Column & Representative term structure \\",
                    r"\midrule",
                ]
            )

            for local_index, signature in enumerate(signature_chunk, start=1):
                lines.append(
                    rf"T{local_index} & "
                    + r"\(\displaystyle "
                    + representative[signature]
                    + r"\) \\"
                )
                lines.append(r"\midrule")

            lines.extend(
                [
                    r"\bottomrule",
                    r"\end{tabular}",
                    r"\end{center}",
                    r"\smallskip",
                ]
            )

            widths = " ".join(
                r">{\raggedright\arraybackslash}p{0.21\linewidth}"
                for _ in signature_chunk
            )
            column_spec = r"@{}llcccccc " + widths + r"@{}"

            headers = [
                r"Model",
                r"$\alpha$",
                r"$d_{S_1}$",
                r"$Y_{S_1}$",
                r"$d_{S_2}$",
                r"$Y_{S_2}$",
                r"$d_F$",
                r"$Y_F$",
                *[
                    rf"T{index}"
                    for index in range(1, len(signature_chunk) + 1)
                ],
            ]

            lines.extend(
                [
                    r"\tiny",
                    rf"\begin{{longtable}}{{{column_spec}}}",
                    r"\toprule",
                    " & ".join(headers) + r" \\",
                    r"\midrule",
                    r"\endfirsthead",
                    r"\toprule",
                    " & ".join(headers) + r" \\",
                    r"\midrule",
                    r"\endhead",
                ]
            )

            for row_index, record in enumerate(successful_records):
                d_s1, y_s1, d_s2, y_s2, d_f, y_f = record_quantum_numbers(
                    record
                )
                grouped = model_terms.get(row_index, {})

                cells = [
                    latex_escape_text(record.name),
                    rf"${record.alpha}$",
                    rf"${d_s1}$",
                    rf"${latex_fraction(y_s1)}$",
                    rf"${d_s2}$",
                    rf"${latex_fraction(y_s2)}$",
                    rf"${d_f}$",
                    rf"${latex_fraction(y_f)}$",
                    *[
                        _rge_term_cell(grouped.get(signature))
                        for signature in signature_chunk
                    ],
                ]

                lines.append(" & ".join(cells) + r" \\")
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

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"\nFinal-EFT RGE term-comparison report:\n{output_path}")
    return output_path


def write_and_compile_final_eft_rge_comparison(
    records: list[RunRecord],
) -> Path:
    """Write and compile the final-EFT Weinberg RGE comparison report."""

    report_tex = write_final_eft_rge_comparison(records)
    compile_latex_document(report_tex)
    return report_tex
