from __future__ import annotations

"""Stage-aware analytic group-factor comparison reports for the T3 pipeline.

The report layer is downstream of the physics implementations.  All numerical
group factors are read from the analytic modules in ``RGE.group_factors``.
The report uses a single paper-style notation throughout:

* numbered scalar couplings and masses are used consistently in every table;
* one beta function is compared across runs in one table;
* run/model labels are columns, while RGE structures are rows;
* EFT1 is named ``EFT_1_after_F`` because only F has been integrated out;
* the F-threshold Wilson tensor is reconstructed from the exported seed rather
  than hard-coded from A--E benchmark values.

For wide studies the run columns are split into chunks so that the generated
PDF remains readable.
"""

from dataclasses import asdict
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import sympy as sp

from common.Records import RunRecord
from Reports.StageReports import group_factor_report_path
from Reports.ReportGeneration import (
    compile_latex_document,
    latex_escape_text,
    paper_notation_key_lines,
    split_latex_terms,
)
from RGE.group_factors.DirectWeinbergGroupFactors import direct_weinberg_group_factor
from RGE.group_factors.FermionMassGroupFactors import fermion_mass_group_factors
from RGE.group_factors.FullYukawaBetaGroupFactors import complete_yukawa_group_factors
from RGE.group_factors.PortalQuarticBetaGroupFactors import (
    portal_quartic_beta_group_factors,
)
from RGE.group_factors.MixingQuarticGroupFactors import (
    mixing_quartic_group_factors,
)
from RGE.group_factors.ScalarMassGroupFactors import scalar_mass_group_factors
from RGE.running.EFT1WilsonAdapter import build_eft1_wilson_tensor
from Reports.RGEComparison import (
    load_eft1_renormalisable_rge_payload,
    load_uv_rge_payload,
    normalise_rgbeta_latex,
)


MAX_RUN_COLUMNS = 5


def _tex_expr(value: object) -> str:
    """Render exact rational/SymPy/string values as mathematical LaTeX."""
    try:
        return sp.latex(sp.sympify(str(value), locals={"sqrt": sp.sqrt}))
    except (sp.SympifyError, TypeError, ValueError):
        return latex_escape_text(value)


def _model_math_label(record: RunRecord) -> str:
    """Compact mathematical run label used as a table column heading."""
    name = latex_escape_text(record.name)
    return rf"\mathrm{{{name}}}\;(\alpha={record.alpha})"


def _run_key(record: RunRecord) -> str:
    """Return a key unique across hypercharge and dimension comparison scans.

    ``record.name`` alone is not unique in the hypercharge study: the five
    alpha values for T3-A all have the same model name.  Using the full
    representation/hypercharge assignment prevents later runs from overwriting
    earlier columns in the comparison dictionaries.
    """
    return (
        f"{record.name}|"
        f"dS1={record.d_s1}|dS2={record.d_s2}|dF={record.d_f}|"
        f"alpha={record.alpha}"
    )


def _hypercharges(record: RunRecord) -> tuple[sp.Rational, sp.Rational, sp.Rational]:
    alpha = int(record.alpha)
    return (
        sp.Rational(alpha, 2),
        sp.Rational(alpha + 2, 2),
        sp.Rational(alpha + 1, 2),
    )


def _valid_records(records: Iterable[RunRecord]) -> list[RunRecord]:
    return [record for record in records if record.summary.get("BuildStatus") == "Success"]


def _chunks(values: Sequence[RunRecord], size: int = MAX_RUN_COLUMNS):
    for start in range(0, len(values), size):
        yield list(values[start : start + size])


def _union_keys(per_run: Mapping[str, Mapping[str, object]]) -> list[str]:
    """Preserve first-seen physics ordering while forming a union of structures."""
    keys: list[str] = []
    seen: set[str] = set()
    for values in per_run.values():
        for key in values:
            if key not in seen:
                seen.add(key)
                keys.append(key)
    return keys


_TOKEN_TEX = {
    "lambdaH": r"\lambda_1",
    "lambdaS1": r"\lambda_2^{(1)}",
    "lambdaS2": r"\lambda_2^{(2)}",
    "lambdaH1": r"\lambda_3^{(1)}",
    "lambdaH2": r"\lambda_3^{(2)}",
    "lambdaH1Adj": r"\lambda_4^{(1)}",
    "lambdaH2Adj": r"\lambda_4^{(2)}",
    "lambdaT3": r"\lambda_5",
    "lambda12": r"\lambda_6",
    "lambda12Adj": r"\lambda_6^{(A)}",
    "lambda12Cross": r"\lambda_6^{(\times)}",
    "lambdaS1Adj": r"\lambda_2^{(1,A)}",
    "lambdaS2Adj": r"\lambda_2^{(2,A)}",
    "mS1Sq": r"m_1^{2}",
    "mS2Sq": r"m_2^{2}",
    "g2_sq": r"g_2^2",
    "gY_sq": r"g_1^2",
    "g2_four": r"g_2^4",
    "gY_four": r"g_1^4",
    "Tr_y1": r"T_{\nu}^{(1)}",
    "Tr_y2": r"T_{\nu}^{(2)}",
    "Tr_Yd": r"\operatorname{Tr}(Y_d^\dagger Y_d)",
    "Tr_Ye": r"\operatorname{Tr}(Y_e^\dagger Y_e)",
    "Tr_Yu": r"\operatorname{Tr}(Y_u^\dagger Y_u)",
    "T": r"T",
    "Tr_y1Ye": r"\operatorname{Tr}(y_1^\dagger y_1Y_e^\dagger Y_e)",
    "Tr_y2Ye": r"\operatorname{Tr}(y_2^\dagger y_2Y_e^\dagger Y_e)",
    "Tr_y1y2": r"\operatorname{Tr}(y_1^\dagger y_1y_2^\dagger y_2)",
    "Tr_MF_y1": r"M_F^\dagger M_F\,\operatorname{Tr}(y_1^\dagger y_1)",
    "Tr_MF_y2": r"M_F^\dagger M_F\,\operatorname{Tr}(y_2^\dagger y_2)",
    "abs_lambdaT3_sq": r"|\lambda_5|^2",
}


_SPECIAL_STRUCTURE_TEX = {
    "trace": r"\mathcal{G}_{\mathrm{tr}}",
    "self_matrix": r"\mathcal{G}_{\mathrm{self}}",
    "cross_matrix": r"\mathcal{G}_{\mathrm{cross}}",
    "charged_lepton_matrix": r"\mathcal{G}_{Y_e}",
    "su2_gauge": r"\mathcal{G}_{SU(2)_L}",
    "u1_gauge": r"\mathcal{G}_{U(1)_Y}",
}


def _structure_tex(key: str) -> str:
    """Convert implementation structure names into mathematical LaTeX."""
    if key in _SPECIAL_STRUCTURE_TEX:
        return _SPECIAL_STRUCTURE_TEX[key]

    if key in _TOKEN_TEX:
        return _TOKEN_TEX[key]

    # Generated non-singlet sectors use ``beta_x:term``.
    if ":" in key:
        beta_name, term = key.split(":", 1)
        beta_tex = _beta_name_tex(beta_name.removeprefix("beta_"))
        return rf"{beta_tex}\;:\;{_structure_tex(term)}"

    pieces = key.split("*")
    rendered: list[str] = []
    for piece in pieces:
        if piece.endswith("_sq") and piece not in _TOKEN_TEX:
            base = piece[:-3]
            rendered.append(rf"\left({_structure_tex(base)}\right)^2")
        else:
            rendered.append(_TOKEN_TEX.get(piece, rf"\mathrm{{{latex_escape_text(piece)}}}"))
    return r"\,".join(rendered)


def _beta_name_tex(name: str) -> str:
    mapping = {
        "y1": r"\beta_{y_1}",
        "y2": r"\beta_{y_2}",
        "MF": r"\beta_{M_F}",
        "mS1Sq": r"\beta_{m_1^{2}}",
        "mS2Sq": r"\beta_{m_2^{2}}",
        "lambdaH": r"\beta_{\lambda_1}",
        "lambdaS1": r"\beta_{\lambda_2^{(1)}}",
        "lambdaS2": r"\beta_{\lambda_2^{(2)}}",
        "lambdaH1": r"\beta_{\lambda_3^{(1)}}",
        "lambdaH2": r"\beta_{\lambda_3^{(2)}}",
        "lambdaH1Adj": r"\beta_{\lambda_4^{(1)}}",
        "lambdaH2Adj": r"\beta_{\lambda_4^{(2)}}",
        "lambdaT3": r"\beta_{\lambda_5}",
        "lambda12": r"\beta_{\lambda_6}",
        "lambda12Adj": r"\beta_{\lambda_6^{(A)}}",
        "lambda12Cross": r"\beta_{\lambda_6^{(\times)}}",
        "lambdaS1Adj": r"\beta_{\lambda_2^{(1,A)}}",
        "lambdaS2Adj": r"\beta_{\lambda_2^{(2,A)}}",
    }
    return mapping.get(name, rf"\beta_{{\mathrm{{{latex_escape_text(name)}}}}}")
def _as_rational(value: object) -> sp.Rational | None:
    """Return an exact rational when a group-factor value is numeric."""
    try:
        expr = sp.sympify(str(value))
    except (sp.SympifyError, TypeError, ValueError):
        return None
    return sp.Rational(expr) if expr.is_Rational else None


def _compact_trace_group_factors(values: Mapping[str, object]) -> dict[str, object]:
    """Replace the paper combination 3 Tr(Yd)+Tr(Ye)+3 Tr(Yu) by T.

    This is applied only when all three coefficients have exactly the ratio
    required by the definition of T, so no physics information is discarded.
    """
    result = dict(values)
    keys = list(result)

    for ye_key in keys:
        if "Tr_Ye" not in ye_key:
            continue
        yd_key = ye_key.replace("Tr_Ye", "Tr_Yd")
        yu_key = ye_key.replace("Tr_Ye", "Tr_Yu")
        if yd_key not in result or yu_key not in result:
            continue

        c_e = _as_rational(result[ye_key])
        c_d = _as_rational(result[yd_key])
        c_u = _as_rational(result[yu_key])
        if c_e is None or c_d is None or c_u is None:
            continue
        if c_d != 3 * c_e or c_u != 3 * c_e:
            continue

        t_key = ye_key.replace("Tr_Ye", "T")
        result[t_key] = c_e
        del result[ye_key]
        del result[yd_key]
        del result[yu_key]

    return result


def _comparison_table(
    beta_tex: str,
    records: Sequence[RunRecord],
    per_record_values: Mapping[str, Mapping[str, object]],
) -> list[str]:
    """One beta function, structures in rows and run coefficients in columns."""
    if not records:
        return []

    lines = [rf"\subsection*{{$ {beta_tex} $}}"]
    compact_values = {
        run_key: _compact_trace_group_factors(values)
        for run_key, values in per_record_values.items()
    }
    all_keys = _union_keys(compact_values)

    for chunk_index, chunk in enumerate(_chunks(records), start=1):
        if len(records) > MAX_RUN_COLUMNS:
            lines.append(
                rf"\textit{{Run block {chunk_index}: "
                rf"{latex_escape_text(chunk[0].name)}--{latex_escape_text(chunk[-1].name)}}}"
            )

        colspec = "@{}p{0.31\\linewidth}" + "c" * len(chunk) + "@{}"
        lines.append(rf"\begin{{longtable}}{{{colspec}}}")
        lines.append(r"\toprule")
        header = [r"$\mathcal{S}$"] + [
            rf"$ {_model_math_label(record)} $" for record in chunk
        ]
        lines.append(" & ".join(header) + r" \\")
        lines.append(r"\midrule")

        for key in all_keys:
            row = [rf"$ {_structure_tex(key)} $"]
            for record in chunk:
                values = compact_values.get(_run_key(record), {})
                if key in values:
                    row.append(rf"$ {_tex_expr(values[key])} $")
                else:
                    row.append(r"$0$")
            lines.append(" & ".join(row) + r" \\")

        lines.append(r"\bottomrule")
        lines.append(r"\end{longtable}")

    return lines


def _mixing_quartic_terms(mix, *, include_heavy_yukawa: bool) -> dict[str, object]:
    terms: dict[str, object] = {
        "lambdaT3*Tr_Yd": mix.yd_trace,
        "lambdaT3*Tr_Ye": mix.ye_trace,
        "lambdaT3*Tr_Yu": mix.yu_trace,
        "lambdaT3*lambdaH": mix.lambdaH,
        "lambdaT3*lambdaH1": mix.lambdaH1,
        "lambdaT3*lambdaH2": mix.lambdaH2,
        "lambdaT3*lambda12": mix.lambda12,
        "g2_sq*lambdaT3": mix.su2_gauge,
        "gY_sq*lambdaT3": mix.u1_gauge,
    }
    if include_heavy_yukawa:
        terms["lambdaT3*Tr_y1"] = mix.y1_trace
        terms["lambdaT3*Tr_y2"] = mix.y2_trace
    if mix.has_lambdaH1Adj:
        terms["lambdaT3*lambdaH1Adj"] = mix.lambdaH1Adj
    if mix.has_lambdaH2Adj:
        terms["lambdaT3*lambdaH2Adj"] = mix.lambdaH2Adj
    if mix.has_lambda12Adj:
        terms["lambdaT3*lambda12Adj"] = mix.lambda12Adj
    if mix.has_lambda12Cross:
        terms["lambdaT3*lambda12Cross"] = mix.lambda12Cross
    return terms


def _uv_group_factor_data(record: RunRecord) -> dict[str, dict[str, object]]:
    y = complete_yukawa_group_factors(
        dS1=record.d_s1,
        dS2=record.d_s2,
        dF=record.d_f,
        alpha=record.alpha,
    )
    mf = fermion_mass_group_factors(
        record.d_s1, record.d_s2, record.d_f, record.alpha
    )
    ms = scalar_mass_group_factors(
        record.d_s1, record.d_s2, record.d_f, record.alpha
    )
    portal = portal_quartic_beta_group_factors(
        record.d_s1, record.d_s2, record.d_f, record.alpha
    )
    mix = mixing_quartic_group_factors(
        record.d_s1, record.d_s2, record.d_f, record.alpha
    )

    return {
        "y1": asdict(y.y1),
        "y2": asdict(y.y2),
        "MF": {
            "y1_left": mf.y1_left,
            "y2_right": mf.y2_right,
            "su2_gauge": mf.su2_gauge,
            "u1_gauge": mf.u1_gauge,
        },
        "mS1Sq": dict(ms.beta_mS1Sq),
        "mS2Sq": dict(ms.beta_mS2Sq),
        "lambdaH1": dict(portal.beta_lambdaH1),
        "lambdaH2": dict(portal.beta_lambdaH2),
        "lambda12": dict(portal.beta_lambda12),
        "lambdaT3": _mixing_quartic_terms(mix, include_heavy_yukawa=True),
        "generated_non_singlet": dict(portal.generated_non_singlet),
    }


def _mf_structure_tex(key: str) -> str:
    mapping = {
        "y1_left": r"(y_1^T y_1^*)M_F",
        "y2_right": r"M_F(y_2^\dagger y_2)",
        "su2_gauge": r"g_2^2M_F",
        "u1_gauge": r"g_1^2M_F",
    }
    return mapping.get(key, _structure_tex(key))


def _comparison_table_mf(
    records: Sequence[RunRecord],
    per_record_values: Mapping[str, Mapping[str, object]],
) -> list[str]:
    if not records:
        return []

    lines = [r"\subsection*{$ \beta_{M_F} $}"]
    all_keys = _union_keys(per_record_values)

    for chunk_index, chunk in enumerate(_chunks(records), start=1):
        if len(records) > MAX_RUN_COLUMNS:
            lines.append(rf"\textit{{Run block {chunk_index}}}")

        colspec = "@{}p{0.31\\linewidth}" + "c" * len(chunk) + "@{}"
        lines.append(rf"\begin{{longtable}}{{{colspec}}}")
        lines.append(r"\toprule")
        header = [r"$\mathcal{S}$"] + [
            rf"$ {_model_math_label(record)} $" for record in chunk
        ]
        lines.append(" & ".join(header) + r" \\")
        lines.append(r"\midrule")

        for key in all_keys:
            row = [rf"$ {_mf_structure_tex(key)} $"]
            for record in chunk:
                values = per_record_values.get(_run_key(record), {})
                row.append(
                    rf"$ {_tex_expr(values[key])} $" if key in values else r"$0$"
                )
            lines.append(" & ".join(row) + r" \\")

        lines.append(r"\bottomrule")
        lines.append(r"\end{longtable}")

    return lines


def _eft1_active_scalar_terms(values: Mapping[str, object]) -> dict[str, object]:
    """Keep only structures built from fields/couplings active after F decouples."""
    forbidden_tokens = ("y1", "y2", "MF")
    return {
        key: value
        for key, value in values.items()
        if not any(token in key for token in forbidden_tokens)
    }


def _eft1_group_factor_data(record: RunRecord) -> dict[str, dict[str, object]]:
    ms = scalar_mass_group_factors(
        record.d_s1, record.d_s2, record.d_f, record.alpha
    )
    portal = portal_quartic_beta_group_factors(
        record.d_s1, record.d_s2, record.d_f, record.alpha
    )
    mix = mixing_quartic_group_factors(
        record.d_s1, record.d_s2, record.d_f, record.alpha
    )

    return {
        "mS1Sq": _eft1_active_scalar_terms(ms.beta_mS1Sq),
        "mS2Sq": _eft1_active_scalar_terms(ms.beta_mS2Sq),
        "lambdaH1": _eft1_active_scalar_terms(portal.beta_lambdaH1),
        "lambdaH2": _eft1_active_scalar_terms(portal.beta_lambdaH2),
        "lambda12": _eft1_active_scalar_terms(portal.beta_lambda12),
        "lambdaT3": _mixing_quartic_terms(mix, include_heavy_yukawa=False),
        "generated_non_singlet": _eft1_active_scalar_terms(
            portal.generated_non_singlet
        ),
    }


def _representation_table(records: Sequence[RunRecord], *, include_f: bool) -> list[str]:
    lines = [
        r"\subsection*{Field representations}",
        r"\begin{longtable}{@{}lccc" + ("c" if include_f else "") + r"@{}}",
        r"\toprule",
    ]
    header = [r"run", r"$S_1$", r"$S_2$", r"$\alpha$"]
    if include_f:
        header.insert(3, r"$F$")
    lines.append(" & ".join(header) + r" \\")
    lines.append(r"\midrule")

    for record in records:
        ys1, ys2, yf = _hypercharges(record)
        row = [
            rf"$\mathrm{{{latex_escape_text(record.name)}}}$",
            rf"$({record.d_s1},{sp.latex(ys1)})$",
            rf"$({record.d_s2},{sp.latex(ys2)})$",
        ]
        if include_f:
            row.append(rf"$({record.d_f},{sp.latex(yf)})$")
        row.append(rf"${record.alpha}$")
        lines.append(" & ".join(row) + r" \\")

    lines.extend([r"\bottomrule", r"\end{longtable}"])
    return lines


def _matching_tensor_rows(
    record: RunRecord,
    max_entries: int = 10,
) -> tuple[list[str], str]:
    seed_path = record.output_dir / "data" / "eft1_after_F_wilson_seed.json"
    if not seed_path.is_file():
        return [], "No exported EFT1 Wilson seed was available for this run."

    import json

    seed = json.loads(seed_path.read_text(encoding="utf-8"))
    tensor = build_eft1_wilson_tensor(
        seed,
        d_s1=record.d_s1,
        d_s2=record.d_s2,
        project_operator_symmetry=True,
    )
    items = sorted(tensor.nonzero_items())
    rows: list[str] = []
    for (i, j, a, b), value in items[:max_entries]:
        rows.append(
            rf"$C_{{{i}{j}{a}{b}}}$ & $ {sp.latex(sp.simplify(value))} $ \\"
        )
    note = (
        f"{len(items)} non-zero real-basis tensor components reconstructed; "
        f"showing the first {min(max_entries, len(items))}."
    )
    return rows, note


def _direct_weinberg_comparison(records: Sequence[RunRecord]) -> list[str]:
    values: dict[str, dict[str, object]] = {}
    for record in records:
        rw = direct_weinberg_group_factor(record.d_s1, record.d_s2, record.d_f)
        values[_run_key(record)] = {"RW": rw.reduced_factor}

    lines = [
        r"\subsection*{Direct Weinberg mixing}",
        r"The active EFT1 scalar interaction produces",
        r"\["
        r"16\pi^2\,\beta_{C_5}\supset "
        r"R_W\,\lambda_5\,C_{LLS_1S_2}."
        r"\]",
    ]

    for chunk in _chunks(records):
        colspec = "@{}l" + "c" * len(chunk) + "@{}"
        lines.append(rf"\begin{{longtable}}{{{colspec}}}")
        lines.append(r"\toprule")
        header = [r"$\mathcal{G}$"] + [
            rf"$ {_model_math_label(record)} $" for record in chunk
        ]
        lines.append(" & ".join(header) + r" \\")
        lines.append(r"\midrule")
        row = [r"$R_W$"] + [
            rf"$ {_tex_expr(values[_run_key(record)]['RW'])} $" for record in chunk
        ]
        lines.append(" & ".join(row) + r" \\")
        lines.extend([r"\bottomrule", r"\end{longtable}"])

    lines.extend(
        [
            r"The displayed $R_W$ values use the validated "
            r"Wigner--$6j$ reduction and inherited tensor convention.  The "
            r"overall tensor-level normalization/phase has not yet been "
            r"independently derived from first principles.",
        ]
    )
    return lines


def _extract_braced(text: str, start: int) -> tuple[str, int] | None:
    """Extract one balanced {...} group starting at ``start``."""
    if start >= len(text) or text[start] != "{":
        return None
    depth = 0
    for index in range(start, len(text)):
        if text[index] == "{":
            depth += 1
        elif text[index] == "}":
            depth -= 1
            if depth == 0:
                return text[start + 1:index], index
    return None


def _numeric_prefactor(term: str) -> tuple[sp.Rational, str]:
    """Split a rendered additive term into exact numeric coefficient and structure."""
    value = term.strip()
    sign = sp.Integer(1)
    if value.startswith("+"):
        value = value[1:].lstrip()
    elif value.startswith("-"):
        sign = -1
        value = value[1:].lstrip()

    # \frac{n}{d} X or \frac{n X}{d}.
    if value.startswith(r"\frac"):
        pos = len(r"\frac")
        first = _extract_braced(value, pos)
        if first is not None:
            numerator, first_end = first
            second_start = first_end + 1
            second = _extract_braced(value, second_start)
            if second is not None:
                denominator, second_end = second
                try:
                    den = sp.Integer(denominator.strip())
                except (TypeError, ValueError):
                    den = None
                if den is not None:
                    num_match = re.match(r"^\s*(\d+)(?:\s+|(?=\\|\{|\())(.*)$", numerator, re.S)
                    if num_match:
                        num = sp.Integer(num_match.group(1))
                        inside_structure = num_match.group(2).strip()
                        remainder = value[second_end + 1:].strip()
                        structure = " ".join(
                            piece for piece in (inside_structure, remainder) if piece
                        )
                        return sp.Rational(sign * num, den), structure
                    if numerator.strip().isdigit():
                        num = sp.Integer(numerator.strip())
                        structure = value[second_end + 1:].strip()
                        return sp.Rational(sign * num, den), structure
                    remainder = value[second_end + 1:].strip()
                    structure = " ".join(
                        piece for piece in (numerator.strip(), remainder) if piece
                    )
                    return sp.Rational(sign, den), structure

    integer = re.match(r"^(\d+)(?:\s+)(.*)$", value, re.S)
    if integer:
        return sign * sp.Integer(integer.group(1)), integer.group(2).strip()

    return sp.Rational(sign), value


def _strip_outer_parentheses(text: str) -> str | None:
    value = text.strip()
    if not (value.startswith("(") and value.endswith(")")):
        return None
    depth = 0
    for index, char in enumerate(value):
        if char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
            if depth == 0 and index != len(value) - 1:
                return None
    return value[1:-1].strip() if depth == 0 else None


def _expanded_rge_terms(latex: str) -> list[tuple[sp.Rational, str]]:
    """Expand simple parenthesised sums after a common numerical prefactor."""
    output: list[tuple[sp.Rational, str]] = []
    for top_term in split_latex_terms(latex) or [latex]:
        outer_coeff, structure = _numeric_prefactor(top_term)
        inner = _strip_outer_parentheses(structure)
        if inner is None:
            output.append((outer_coeff, structure))
            continue
        inner_terms = split_latex_terms(inner)
        if len(inner_terms) <= 1:
            output.append((outer_coeff, structure))
            continue
        for inner_term in inner_terms:
            inner_coeff, inner_structure = _numeric_prefactor(inner_term)
            output.append((sp.Rational(outer_coeff * inner_coeff), inner_structure))
    return output


def _clean_rge_structure_latex(structure: str) -> str:
    """Canonicalize harmless RGBeta TeX formatting for cross-run alignment."""
    value = structure.strip()
    value = value.replace(r"\left", "").replace(r"\right", "")
    value = value.replace(r"\,", " ")
    value = re.sub(r"\{(Y_[ude]|y_[12]|g_[123])\}", r"\1", value)
    value = re.sub(r"\s+", " ", value).strip()
    return value


def _trace_factor(structure: str) -> tuple[str | None, str]:
    """Identify a leading single-species trace and return its remaining structure."""
    value = _clean_rge_structure_latex(structure)
    if not value.startswith(r"\operatorname{Tr}"):
        return None, value

    # RGBeta's normalized form is Tr( Y_x (Y_x)^dagger ) times the rest.
    open_index = value.find("(")
    if open_index < 0:
        return None, value
    depth = 0
    close_index = None
    for index in range(open_index, len(value)):
        if value[index] == "(":
            depth += 1
        elif value[index] == ")":
            depth -= 1
            if depth == 0:
                close_index = index
                break
    if close_index is None:
        return None, value

    trace_body = value[open_index + 1:close_index]
    remainder = value[close_index + 1:].strip()
    species = None
    for candidate in ("Y_d", "Y_e", "Y_u", "y_1", "y_2"):
        if candidate in trace_body:
            others = [
                item for item in ("Y_d", "Y_e", "Y_u", "y_1", "y_2")
                if item != candidate and item in trace_body
            ]
            if not others:
                species = candidate
                break
    return species, remainder


def _compact_saved_trace_terms(
    terms: list[tuple[sp.Rational, str]],
) -> list[tuple[sp.Rational, str]]:
    """Use T and T_nu^(i) whenever the saved RGE terms permit it."""
    ordinary: list[tuple[sp.Rational, str]] = []
    traces: dict[tuple[str, str], sp.Rational] = {}

    for coeff, structure in terms:
        species, remainder = _trace_factor(structure)
        if species is None:
            ordinary.append((coeff, _clean_rge_structure_latex(structure)))
        else:
            key = (species, remainder)
            traces[key] = traces.get(key, sp.Rational(0)) + coeff

    # Combine the exact paper T = Tr(Ye†Ye+3Yu†Yu+3Yd†Yd).
    remainders = {remainder for _, remainder in traces}
    used: set[tuple[str, str]] = set()
    for remainder in sorted(remainders):
        kd, ke, ku = ("Y_d", remainder), ("Y_e", remainder), ("Y_u", remainder)
        if kd in traces and ke in traces and ku in traces:
            ce = traces[ke]
            if traces[kd] == 3 * ce and traces[ku] == 3 * ce:
                structure = r"T" + (r"\," + remainder if remainder else "")
                ordinary.append((ce, structure))
                used.update((kd, ke, ku))

    # The generalized T3 theory has two neutrino-sector Yukawas, so use
    # T_nu^(1,2) rather than inventing one aggregate T_nu that is not universal.
    for (species, remainder), coeff in traces.items():
        if (species, remainder) in used:
            continue
        if species == "y_1":
            structure = r"T_{\nu}^{(1)}" + (r"\," + remainder if remainder else "")
        elif species == "y_2":
            structure = r"T_{\nu}^{(2)}" + (r"\," + remainder if remainder else "")
        else:
            trace_symbol = {
                "Y_d": r"\operatorname{Tr}(Y_d^\dagger Y_d)",
                "Y_e": r"\operatorname{Tr}(Y_e^\dagger Y_e)",
                "Y_u": r"\operatorname{Tr}(Y_u^\dagger Y_u)",
            }[species]
            structure = trace_symbol + (r"\," + remainder if remainder else "")
        ordinary.append((coeff, structure))

    return ordinary


def _saved_rge_term_table(
    beta_tex: str,
    records: Sequence[RunRecord],
    per_run: Mapping[str, Mapping[str, object]],
) -> list[str]:
    """Render saved RGBeta beta functions as structure/coefficients tables."""
    if not records:
        return []

    all_structures = _union_keys(per_run)
    lines = [rf"\subsection*{{$ {beta_tex} $}}"]

    for chunk in _chunks(records):
        colspec = "@{}p{0.38\\linewidth}" + "c" * len(chunk) + "@{}"
        lines.extend([
            rf"\begin{{longtable}}{{{colspec}}}",
            r"\toprule",
            " & ".join(
                [r"$\mathcal{S}$"]
                + [rf"$ {_model_math_label(record)} $" for record in chunk]
            ) + r" \\",
            r"\midrule",
        ])

        for structure in all_structures:
            row = [rf"$ {structure} $"]
            for record in chunk:
                value = per_run.get(_run_key(record), {}).get(structure)
                row.append(r"$0$" if value is None else rf"$ {_tex_expr(value)} $")
            lines.append(" & ".join(row) + r" \\")

        lines.extend([r"\bottomrule", r"\end{longtable}"])

    return lines


def _saved_rge_beta_sections(
    records: Sequence[RunRecord],
    *,
    loader,
    couplings: Sequence[str],
    title: str,
) -> list[str]:
    """Render gauge/Yukawa RGBeta results in the same comparison-table format."""
    symbol_map = {
        "gY": r"g_1",
        "g2": r"g_2",
        "g3": r"g_3",
        "yu": r"Y_u",
        "yd": r"Y_d",
        "ye": r"Y_e",
        "y1": r"y_1",
        "y2": r"y_2",
    }

    lines = [rf"\section*{{{title}}}"]

    for coupling in couplings:
        per_run: dict[str, dict[str, object]] = {}
        have_any = False

        for record in records:
            run_key = _run_key(record)
            payload = loader(record)
            per_run[run_key] = {}
            if payload is None:
                continue
            latex = (payload.get("report_beta_latex", {}) or {}).get(coupling)
            if not latex or not str(latex).strip():
                continue

            normalized = normalise_rgbeta_latex(str(latex))
            terms = _compact_saved_trace_terms(_expanded_rge_terms(normalized))
            for coefficient, structure in terms:
                structure = _clean_rge_structure_latex(structure)
                per_run[run_key][structure] = (
                    sp.Rational(per_run[run_key].get(structure, 0))
                    + sp.Rational(coefficient)
                )
                have_any = True

        if not have_any:
            continue

        symbol = symbol_map[coupling]
        lines.extend(
            _saved_rge_term_table(
                rf"\beta_{{{symbol}}}^{{(1)}}",
                records,
                per_run,
            )
        )

    return lines


def _weinberg_index_conventions() -> list[str]:
    """Explain flavour versus scalar indices used by the EFT1 Wilson tensor."""
    return [
        r"\section*{Weinberg and scalar-index conventions}",
        r"The index sets in the EFT1 Wilson coefficient have different physical meanings:",
        r"\begin{itemize}",
        r"\item $i,j$ are lepton-flavour indices.  These are the indices that survive "
        r"into the low-energy Weinberg coefficient $C_5^{ij}$.",
        r"\item In a complex irreducible-representation basis, $A$ labels a component "
        r"of $S_1$ and $B$ labels a component of $S_2$.  The operator coefficient may "
        r"therefore be written schematically as $C_{LLS_1S_2}^{ij;AB}$.",
        r"\item The EFT1 Wilson adapter converts the scalar fields to a global real-scalar "
        r"basis.  The stored tensor is then $C_{ijab}$.  The indices $a,b$ label real "
        r"scalar components in that basis; they are neither flavour indices nor generator indices.",
        r"\item At the scalar threshold the scalar indices are contracted through the "
        r"$\lambda_5$ interaction and the $SU(2)_L$ recoupling.  The result is the "
        r"scalar-free SMEFT coefficient $C_5^{ij}$, symmetric in $i,j$.",
        r"\end{itemize}",
    ]


def _matching_tensor_appendix(records: Sequence[RunRecord]) -> list[str]:
    """Keep basis-dependent tensor samples separate from comparison beta tables."""
    lines = [
        r"\section*{F-threshold matching tensor samples}",
        r"The component bases have representation-dependent dimensions, so the "
        r"raw tensor entries are not placed in the cross-run beta-function tables.",
    ]
    for record in records:
        rows, note = _matching_tensor_rows(record)
        lines.append(rf"\subsection*{{$ {_model_math_label(record)} $}}")
        lines.append(rf"\textit{{{latex_escape_text(note)}}}")
        if rows:
            lines.extend(
                [
                    r"\begin{longtable}{@{}ll@{}}",
                    r"\toprule",
                    r"$C_{ijab}$ & matched coefficient \\",
                    r"\midrule",
                    *rows,
                    r"\bottomrule",
                    r"\end{longtable}",
                ]
            )
    return lines


def _generated_non_singlet_section(
    records: Sequence[RunRecord],
    stage_data: Mapping[str, Mapping[str, Mapping[str, object]]],
) -> list[str]:
    """Render one comparison table for each generated extra quartic coupling.

    The analytic portal module stores these entries as
    ``beta_<coupling>:<structure>``.  Here they are regrouped so every beta
    function has structures in rows and models/runs in columns, matching the
    rest of the report.
    """
    coupling_order = (
        "lambdaH1Adj",
        "lambdaH2Adj",
        "lambda12Adj",
        "lambda12Cross",
        "lambdaS1Adj",
        "lambdaS2Adj",
    )

    grouped: dict[str, dict[str, dict[str, object]]] = {
        coupling: {_run_key(record): {} for record in records}
        for coupling in coupling_order
    }

    discovered: list[str] = []

    for record in records:
        run_key = _run_key(record)
        flat = stage_data[run_key].get("generated_non_singlet", {})
        for full_key, value in flat.items():
            if ":" not in full_key:
                continue
            beta_key, structure = full_key.split(":", 1)
            if not beta_key.startswith("beta_"):
                continue
            coupling = beta_key.removeprefix("beta_")
            if coupling not in grouped:
                grouped[coupling] = {
                    _run_key(item): {}
                    for item in records
                }
            grouped[coupling][run_key][structure] = value
            if coupling not in discovered:
                discovered.append(coupling)

    ordered = [
        coupling
        for coupling in coupling_order
        if coupling in discovered
    ] + [
        coupling
        for coupling in discovered
        if coupling not in coupling_order
    ]

    if not ordered:
        return []

    lines = [r"\section*{Additional scalar-coupling beta functions}"]
    for coupling in ordered:
        per_run = grouped[coupling]
        if not any(per_run.values()):
            continue
        lines.extend(
            _comparison_table(
                _beta_name_tex(coupling),
                records,
                per_run,
            )
        )
    return lines


def _notation_key() -> list[str]:
    """Shared field-contraction notation and index key."""
    return paper_notation_key_lines()


def _document_header(title: str, stage_text: str) -> list[str]:
    return [
        r"\documentclass[9pt]{article}",
        r"\usepackage[a4paper,margin=1.15cm]{geometry}",
        r"\usepackage{amsmath,amssymb,booktabs,longtable,array}",
        r"\usepackage[T1]{fontenc}",
        r"\setlength{\parindent}{0pt}",
        r"\setlength{\parskip}{0.45em}",
        r"\setlength{\tabcolsep}{3.5pt}",
        r"\renewcommand{\arraystretch}{1.16}",
        r"\begin{document}",
        rf"\section*{{{title}}}",
        stage_text,
        r"\["
        r"j=\frac{d-1}{2},\qquad "
        r"C_2(d)=\frac{d^2-1}{4},\qquad "
        r"T(d)=\frac{d(d^2-1)}{12}."
        r"\]",
        r"In every comparison table, $\mathcal{S}$ denotes the RGE structure "
        r"multiplying the displayed group factor in $16\pi^2\beta$.",
    ]


def write_gf_uv(records: list[RunRecord], *, report_root: Path) -> Path:
    valid = _valid_records(records)
    path = group_factor_report_path("UV", report_root=report_root)
    path.parent.mkdir(parents=True, exist_ok=True)

    lines = _document_header(
        "T3 group factors: UV",
        (
            r"Active theory: $\mathrm{SM}+S_1+S_2+F$.  "
            r"The same beta function is shown across runs in a single comparison "
            r"table, with exact analytic group factors in the run columns."
        ),
    )
    lines.extend(_notation_key())
    lines.extend(_representation_table(valid, include_f=True))
    lines.extend(
        _saved_rge_beta_sections(
            valid,
            loader=load_uv_rge_payload,
            couplings=("gY", "g2", "g3", "yu", "yd", "ye", "y1", "y2"),
            title="Gauge and Yukawa beta functions",
        )
    )

    uv_data = {_run_key(record): _uv_group_factor_data(record) for record in valid}

    lines.extend(
        _comparison_table(
            r"\beta_{y_1}",
            valid,
            {_run_key(record): uv_data[_run_key(record)]["y1"] for record in valid},
        )
    )
    lines.extend(
        _comparison_table(
            r"\beta_{y_2}",
            valid,
            {_run_key(record): uv_data[_run_key(record)]["y2"] for record in valid},
        )
    )
    lines.extend(
        _comparison_table_mf(
            valid,
            {_run_key(record): uv_data[_run_key(record)]["MF"] for record in valid},
        )
    )

    for key in ("mS1Sq", "mS2Sq", "lambdaH1", "lambdaH2", "lambda12", "lambdaT3"):
        lines.extend(
            _comparison_table(
                _beta_name_tex(key),
                valid,
                {_run_key(record): uv_data[_run_key(record)][key] for record in valid},
            )
        )

    lines.extend(_generated_non_singlet_section(valid, uv_data))

    lines.extend(
        [
            r"\section*{Threshold interpretation}",
            r"At $\mu=M_F$, the fermion $F$ is integrated out.  Therefore "
            r"$y_1$, $y_2$ and $M_F$ cease to be active EFT parameters and "
            r"their information is transferred to the matched "
            r"$C_{LLS_1S_2}$ boundary tensor.",
            r"\end{document}",
            "",
        ]
    )

    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def write_gf_eft1(records: list[RunRecord], *, report_root: Path) -> Path | None:
    relevant: list[RunRecord] = []
    for record in _valid_records(records):
        stages = record.summary.get("EFTStages", [])
        if (
            stages
            and stages[0].get("IntegratedFields") == ["F"]
            and set(stages[0].get("ActiveHeavyFields", [])) == {"S1", "S2"}
        ):
            relevant.append(record)

    if not relevant:
        return None

    # Only F has been integrated out at this stage.
    path = group_factor_report_path("EFT_1_after_F", report_root=report_root)
    path.parent.mkdir(parents=True, exist_ok=True)

    lines = _document_header(
        "T3 group factors: EFT1 after integrating out $F$",
        (
            r"Active theory: $\mathrm{SM}+S_1+S_2+C_{LLS_1S_2}$.  "
            r"The fermion $F$ is absent, so $y_1$, $y_2$ and $M_F$ are not "
            r"running EFT couplings."
        ),
    )
    lines.extend(_notation_key())
    lines.extend(_representation_table(relevant, include_f=False))
    lines.extend(
        _saved_rge_beta_sections(
            relevant,
            loader=load_eft1_renormalisable_rge_payload,
            couplings=("gY", "g2", "g3", "yu", "yd", "ye"),
            title="Gauge and Yukawa beta functions in EFT1",
        )
    )

    eft1_data = {
        _run_key(record): _eft1_group_factor_data(record)
        for record in relevant
    }

    for key in ("mS1Sq", "mS2Sq", "lambdaH1", "lambdaH2", "lambda12", "lambdaT3"):
        lines.extend(
            _comparison_table(
                _beta_name_tex(key),
                relevant,
                {_run_key(record): eft1_data[_run_key(record)][key] for record in relevant},
            )
        )

    lines.extend(_generated_non_singlet_section(relevant, eft1_data))
    lines.extend(_weinberg_index_conventions())
    lines.extend(_direct_weinberg_comparison(relevant))
    lines.extend(_matching_tensor_appendix(relevant))

    lines.extend(
        [
            r"\section*{Next threshold}",
            r"At the next threshold, $S_1$ and $S_2$ are integrated out.  "
            r"Their representation dependence is then absorbed into the "
            r"matching value of $C_5$; below that threshold the SMEFT "
            r"Weinberg-operator anomalous dimension is universal.",
            r"\end{document}",
            "",
        ]
    )

    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def write_and_compile_stage_group_factor_reports(
    records: list[RunRecord],
    *,
    report_root: Path,
) -> list[Path]:
    """Write and compile the implemented stage-aware GroupFactors reports."""
    outputs: list[Path] = []

    uv = write_gf_uv(records, report_root=report_root)
    compile_latex_document(uv)
    outputs.append(uv)

    eft1 = write_gf_eft1(records, report_root=report_root)
    if eft1 is not None:
        compile_latex_document(eft1)
        outputs.append(eft1)

    return outputs
