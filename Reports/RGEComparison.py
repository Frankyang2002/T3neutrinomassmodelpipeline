from __future__ import annotations

"""Across-model comparison report for the UV one-loop T3 RGEs from RGBeta."""

import json
from pathlib import Path
from typing import Any

from common.Paths import REPORT_OUTPUT_DIR
from common.Records import RunRecord
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


def write_rge_comparison(records: list[RunRecord]) -> Path:
    """Create the across-model UV one-loop RGE comparison report."""

    output_path = REPORT_OUTPUT_DIR / "rge_comparison.tex"

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

        if not isinstance(raw_betas, dict):
            raw_betas = {}
        if not isinstance(latex_betas, dict):
            latex_betas = {}

        report_betas = payload.get("report_betas", {})
        if not isinstance(report_betas, dict):
            report_betas = {}

        coupling_names.update(str(name) for name in report_betas or raw_betas)
        coupling_names.update(str(name) for name in latex_betas)
        rows.append((record, payload))

    lines: list[str] = [
        r"\documentclass{article}",
        r"\usepackage[margin=0.8cm]{geometry}",
        r"\usepackage{amsmath,amssymb,adjustbox,pdflscape,longtable,array,booktabs}",
        r"\usepackage[T1]{fontenc}",
        r"\setlength{\tabcolsep}{2pt}",
        r"\renewcommand{\arraystretch}{1.2}",
        r"\begin{document}",
        r"\begin{landscape}",
        r"\section*{T3 UV one-loop RGE comparison}",
        (
            r"Each table is written in the conventional form "
            r"$16\pi^2\,\mu\,dX/d\mu=\beta_X^{(1)}$. "
            r"A dash means that the coupling is not registered for that representation assignment."
        ),
        r"\medskip",
        r"\subsection*{Conventions and the RGBeta gauge-coupling conversion}",
        (
            r"For ordinary couplings the report uses the quantity returned by RGBeta directly. "
            r"For gauge couplings, RGBeta uses $g^2$ as its running variable, so its one-loop output is proportional to "
            r"$d(g^2)/d\ln\mu$. By the chain rule,"
        ),
        r"\[\frac{d(g^2)}{d\ln\mu}=2g\frac{dg}{d\ln\mu}=2g\,\mu\frac{dg}{d\mu}.\]",
        (
            r"Therefore, to present the gauge RGE in the same convention as every other table, "
            r"$16\pi^2\,\mu\,dg/d\mu=\beta_g^{(1)}$, the RGBeta gauge result is divided by $2g$. "
            r"This is only a change of running variable; no physics or loop factor is being altered."
        ),
        r"\medskip",
        r"\subsection*{Notation}",
        r"$S_1$ and $S_2$ denote the two BSM scalar multiplets and $F$ the BSM fermion multiplet. "
        r"The columns $d_X$ give the $SU(2)_L$ representation dimensions and $Y_X$ the hypercharges. "
        r"A dagger denotes Hermitian conjugation, $T$ transpose, $*$ complex conjugation, and $\operatorname{Tr}$ a trace over flavour indices.",
        r"\medskip",
        r"\subsection*{Coupling glossary}",
        r"The interaction terms below are schematic: gauge-index contractions and conjugations are fixed by each model's quantum numbers. "
        r"For representation-dependent quartics, labels such as Adj and Cross distinguish independent $SU(2)_L$ invariant contractions.",
        r"\begin{longtable}{@{}p{0.11\linewidth}p{0.17\linewidth}p{0.37\linewidth}p{0.29\linewidth}@{}}",
        r"\toprule",
        r"Symbol & Meaning & Schematic Lagrangian term & Notes \\",
        r"\midrule",
        r"\endfirsthead",
        r"\toprule",
        r"Symbol & Meaning & Schematic Lagrangian term & Notes \\",
        r"\midrule",
        r"\endhead",
    ]

    for symbol, meaning, lagrangian_term, notes in GLOSSARY_ROWS:
        lines.append(f"{symbol} & {meaning} & {lagrangian_term} & {notes} " + r"\\")
        lines.append(r"\midrule")

    lines.extend([r"\bottomrule", r"\end{longtable}", r"\clearpage"])

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
            lines.extend(
                [
                    rf"\section*{{$\beta_{{{symbol}}}^{{(1)}}$}}",
                    rf"\[16\pi^2\,\mu\frac{{d {symbol}}}{{d\mu}}="
                    rf"\beta_{{{symbol}}}^{{(1)}}\]",
                    r"\tiny",
                    r"\begin{longtable}{@{}llcccccc >{\raggedright\arraybackslash}p{0.54\linewidth}@{}}",
                    r"\toprule",
                    (
                        r"Model & $\alpha$ & $d_{S_1}$ & $Y_{S_1}$ & "
                        r"$d_{S_2}$ & $Y_{S_2}$ & $d_F$ & $Y_F$ & "
                        rf"$\beta_{{{symbol}}}^{{(1)}}$ \\" 
                    ),
                    r"\midrule",
                    r"\endfirsthead",
                    r"\toprule",
                    (
                        r"Model & $\alpha$ & $d_{S_1}$ & $Y_{S_1}$ & "
                        r"$d_{S_2}$ & $Y_{S_2}$ & $d_F$ & $Y_F$ & "
                        rf"$\beta_{{{symbol}}}^{{(1)}}$ \\" 
                    ),
                    r"\midrule",
                    r"\endhead",
                ]
            )

            for record, payload in rows:
                raw_betas = payload.get("betas", {}) or {}
                latex_betas = payload.get("report_beta_latex", {}) or {}
                d_s1, y_s1, d_s2, y_s2, d_f, y_f = record_quantum_numbers(record)

                cell = _beta_cell(
                    latex_betas.get(coupling),
                    raw_betas.get(coupling),
                )

                cells = [
                    latex_escape_text(record.name),
                    rf"${record.alpha}$",
                    rf"${d_s1}$",
                    rf"${latex_fraction(y_s1)}$",
                    rf"${d_s2}$",
                    rf"${latex_fraction(y_s2)}$",
                    rf"${d_f}$",
                    rf"${latex_fraction(y_f)}$",
                    cell,
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
    print(f"\nUV RGE comparison report:\n{output_path}")
    return output_path


def write_and_compile_rge_comparison(records: list[RunRecord]) -> Path:
    """Write and compile the across-model UV RGE comparison report."""

    report_tex = write_rge_comparison(records)
    compile_latex_document(report_tex)
    return report_tex
