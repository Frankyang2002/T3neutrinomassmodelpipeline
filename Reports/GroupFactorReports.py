from __future__ import annotations

import re
"""Stage-aware analytic group-factor comparison reports for the T3 pipeline.

The report layer is downstream of the physics implementations.  All numerical
group factors are read from the analytic modules in ``RGE.group_factors``.
The report uses a single paper-style notation throughout:

* numbered scalar couplings and masses are used consistently in every table;
* one beta function is compared across runs in one table;
* run/model labels are columns, while RGE structures are rows;
* every EFT stage recorded by the threshold pipeline receives a GroupFactors
  report, independent of the threshold ordering;
* the special F-first stage retains the full analytic EFT1 treatment;
* fully decoupled stages display the universal SMEFT Weinberg-operator group
  factors;
* other threshold orderings receive a stage-content/group-theory report without
  reusing beta-function coefficients from a different EFT;
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
        f"shared={int(record.shared_scalar)}|"
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
    "mSSq": r"m_S^{2}",
    "mS1Sq": r"m_1^{2}",
    "mS2Sq": r"m_2^{2}",
    "g2_sq": r"g_2^2",
    "gY_sq": r"g_1^2",
    "g2_four": r"g_2^4",
    "gY_four": r"g_1^4",
    "g2sq_gYsq": r"g_1^2g_2^2",
    "C5": r"C_5",
    "Tr_y1": r"T_{\nu}^{(1)}",
    "Tr_y2": r"T_{\nu}^{(2)}",
    "Tnu": r"T_\nu",
    "Tr_Yd": r"\operatorname{Tr}(Y_d^\dagger Y_d)",
    "Tr_Ye": r"\operatorname{Tr}(Y_e^\dagger Y_e)",
    "Tr_Yu": r"\operatorname{Tr}(Y_u^\dagger Y_u)",
    "T": r"T",
    "Tr_y1Ye": r"\operatorname{Tr}(y_1^\dagger y_1Y_e^\dagger Y_e)",
    "Tr_y2Ye": r"\operatorname{Tr}(y_2^\dagger y_2Y_e^\dagger Y_e)",
    "Tr_y1y2": r"\operatorname{Tr}(y_1^\dagger y_1y_2^\dagger y_2)",
    "Tr_MF_y1": r"M_F^\dagger M_F\,T_{\nu}^{(1)}",
    "Tr_MF_y2": r"M_F^\dagger M_F\,T_{\nu}^{(2)}",
    "abs_lambdaT3_sq": r"|\lambda_5|^2",
}


_SPECIAL_STRUCTURE_TEX = {
    # Context-independent special structures only.
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
        "h": r"\beta_h",
        "mSSq": r"\beta_{m_S^{2}}",
        "lambdaS": r"\beta_{\lambda_2}",
        "lambda3": r"\beta_{\lambda_3}",
        "lambda4": r"\beta_{\lambda_4}",
        "lambda5": r"\beta_{\lambda_5}",
        "mS1Sq": r"\beta_{m_1^{2}}",
        "mS2Sq": r"\beta_{m_2^{2}}",
        "lambdaH": r"\beta_{\lambda_1}",
        "lambdaS1": r"\beta_{\lambda_{S_1}^{(1)}}",
        "lambdaS2": r"\beta_{\lambda_{S_2}^{(1)}}",
        "lambdaH1": r"\beta_{\lambda_{HS_1}^{(1)}}",
        "lambdaH2": r"\beta_{\lambda_{HS_2}^{(1)}}",
        "lambdaH1Adj": r"\beta_{\lambda_{HS_1}^{(A)}}",
        "lambdaH2Adj": r"\beta_{\lambda_{HS_2}^{(A)}}",
        "lambdaT3": r"\beta_{\lambda_5}",
        "lambda12": r"\beta_{\lambda_{12}^{(1)}}",
        "lambda12Adj": r"\beta_{\lambda_{12}^{(A)}}",
        "lambda12Cross": r"\beta_{\lambda_{12}^{(\times)}}",
        "lambdaS1Adj": r"\beta_{\lambda_{S_1}^{(A)}}",
        "lambdaS2Adj": r"\beta_{\lambda_{S_2}^{(A)}}",
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

    keys = list(result)
    for y1_key in keys:
        if "Tr_y1" not in y1_key or y1_key not in result:
            continue
        y2_key = y1_key.replace("Tr_y1", "Tr_y2")
        if y2_key not in result:
            continue
        c1 = _as_rational(result[y1_key])
        c2 = _as_rational(result[y2_key])
        if c1 is None or c2 is None or c1 != c2:
            continue
        result[y1_key.replace("Tr_y1", "Tnu")] = c1
        del result[y1_key]
        del result[y2_key]

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


def _comparison_table_yukawa(
    beta_name: str,
    records: Sequence[RunRecord],
    per_record_values: Mapping[str, Mapping[str, object]],
) -> list[str]:
    """Render analytic y1/y2 group factors against the actual RGE structures."""
    if beta_name == "y1":
        structures = {
            "trace": r"T_{\nu}^{(1)}\,y_1",
            "self_matrix": r"y_1y_1^\dagger y_1",
            "cross_matrix": r"y_2y_2^\dagger y_1",
            "charged_lepton_matrix": r"Y_eY_e^\dagger y_1",
            "su2_gauge": r"g_2^2\,y_1",
            "u1_gauge": r"g_1^2\,y_1",
        }
    elif beta_name == "y2":
        structures = {
            "trace": r"T_{\nu}^{(2)}\,y_2",
            "self_matrix": r"y_2y_2^\dagger y_2",
            "cross_matrix": r"y_1y_1^\dagger y_2",
            "charged_lepton_matrix": r"Y_eY_e^\dagger y_2",
            "su2_gauge": r"g_2^2\,y_2",
            "u1_gauge": r"g_1^2\,y_2",
        }
    else:
        raise ValueError(f"Unsupported Yukawa beta name: {beta_name}")

    lines = [rf"\subsection*{{$ \beta_{{{beta_name}}} $}}"]
    for chunk in _chunks(records):
        colspec = "@{}p{0.31\\linewidth}" + "c" * len(chunk) + "@{}"
        lines.extend([
            rf"\begin{{longtable}}{{{colspec}}}",
            r"\toprule",
            " & ".join(
                [r"$\mathcal{S}$"]
                + [rf"$ {_model_math_label(record)} $" for record in chunk]
            ) + r" \\",
            r"\midrule",
        ])
        for key, structure in structures.items():
            row = [rf"$ {structure} $"]
            for record in chunk:
                value = per_record_values.get(_run_key(record), {}).get(key)
                row.append(r"$0$" if value is None else rf"$ {_tex_expr(value)} $")
            lines.append(" & ".join(row) + r" \\")
        lines.extend([r"\bottomrule", r"\end{longtable}"])
    return lines


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
        rendered_value = sp.latex(sp.simplify(value)).replace("MF", r"M_F")
        rows.append(
            rf"$C_{{{i}{j}{a}{b}}}$ & $ {rendered_value} $ \\"
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
        header = [r"$\mathcal{S}$"] + [
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
    for candidate in ("Y_d", "Y_e", "Y_u", "y_1", "y_2", "h"):
        if candidate in trace_body:
            others = [
                item for item in ("Y_d", "Y_e", "Y_u", "y_1", "y_2", "h")
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

    for remainder in sorted(remainders):
        k1, k2 = ("y_1", remainder), ("y_2", remainder)
        if k1 in traces and k2 in traces and k1 not in used and k2 not in used:
            if traces[k1] == traces[k2]:
                structure = r"T_\nu" + (r"\," + remainder if remainder else "")
                ordinary.append((traces[k1], structure))
                used.update((k1, k2))

    for (species, remainder), coeff in traces.items():
        if (species, remainder) in used:
            continue
        if species == "y_1":
            structure = r"T_{\nu}^{(1)}" + (r"\," + remainder if remainder else "")
        elif species == "y_2":
            structure = r"T_{\nu}^{(2)}" + (r"\," + remainder if remainder else "")
        elif species == "h":
            structure = r"T_\nu" + (r"\," + remainder if remainder else "")
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
        "h": r"h",
        "MF": r"M_F",
        "mSSq": r"m_S^2",
        "lambdaH": r"\lambda_1",
        "lambdaS": r"\lambda_2",
        "lambda3": r"\lambda_3",
        "lambda4": r"\lambda_4",
        "lambda5": r"\lambda_5",
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


def _notation_table(rows: Sequence[tuple[str, str]]) -> list[str]:
    lines = [
        r"\section*{Notation key}",
        r"\begin{longtable}{@{}p{0.23\linewidth}p{0.69\linewidth}@{}}",
        r"\toprule",
        r"symbol & interaction / definition \\",
        r"\midrule",
    ]
    lines.extend(rf"{symbol} & {term} \\" for symbol, term in rows)
    lines.extend([r"\bottomrule", r"\end{longtable}"])
    return lines


def _uv_notation_key() -> list[str]:
    return _notation_table([
        (r"$g_1,g_2,g_3$", r"SM gauge interactions"),
        (r"$Y_u,Y_d,Y_e$", r"SM Yukawa interactions"),
        (r"$y_1$", r"$LFS_1+\mathrm{h.c.}$"),
        (r"$y_2$", r"$LFS_2+\mathrm{h.c.}$"),
        (r"$M_F$", r"$-\bar FM_FF$"),
        (r"$T$", r"$\operatorname{Tr}(Y_e^\dagger Y_e+3Y_u^\dagger Y_u+3Y_d^\dagger Y_d)$"),
        (r"$T_\nu^{(1)}$", r"$\operatorname{Tr}(y_1^\dagger y_1)$"),
        (r"$T_\nu^{(2)}$", r"$\operatorname{Tr}(y_2^\dagger y_2)$"),
        (r"$T_\nu$", r"$T_\nu^{(1)}+T_\nu^{(2)}$"),
        (r"$m_1^2,m_2^2$", r"$m_1^2S_1^\dagger S_1,\ m_2^2S_2^\dagger S_2$"),
        (r"$\lambda_1$", r"$\frac12\lambda_1(H^\dagger H)^2$"),
        (r"$\lambda_{S_i}^{(1)}$", r"$\frac12\lambda_{S_i}^{(1)}(S_i^\dagger S_i)^2$"),
        (r"$\lambda_{HS_i}^{(1)}$", r"$(H^\dagger H)(S_i^\dagger S_i)$"),
        (r"$\lambda_{HS_i}^{(A)}$", r"$(H^\dagger T^AH)(S_i^\dagger T^AS_i)$"),
        (r"$\lambda_5$", r"$HHS_1S_2^\dagger+\mathrm{h.c.}$"),
        (r"$\lambda_{12}^{(1)}$", r"$(S_1^\dagger S_1)(S_2^\dagger S_2)$"),
        (r"$\lambda_{12}^{(A)}$", r"$(S_1^\dagger T^AS_1)(S_2^\dagger T^AS_2)$"),
        (r"$\lambda_{12}^{(\times)}$", r"crossed independent $S_1$--$S_2$ contraction"),
    ])


def _eft1_notation_key() -> list[str]:
    return _notation_table([
        (r"$g_1,g_2,g_3$", r"SM gauge interactions"),
        (r"$Y_u,Y_d,Y_e$", r"SM Yukawa interactions"),
        (r"$T$", r"$\operatorname{Tr}(Y_e^\dagger Y_e+3Y_u^\dagger Y_u+3Y_d^\dagger Y_d)$"),
        (r"$m_1^2,m_2^2$", r"$m_1^2S_1^\dagger S_1,\ m_2^2S_2^\dagger S_2$"),
        (r"$\lambda_1$", r"$\frac12\lambda_1(H^\dagger H)^2$"),
        (r"$\lambda_{S_i}^{(1)}$", r"$\frac12\lambda_{S_i}^{(1)}(S_i^\dagger S_i)^2$"),
        (r"$\lambda_{HS_i}^{(1)}$", r"$(H^\dagger H)(S_i^\dagger S_i)$"),
        (r"$\lambda_{HS_i}^{(A)}$", r"$(H^\dagger T^AH)(S_i^\dagger T^AS_i)$"),
        (r"$\lambda_5$", r"$HHS_1S_2^\dagger+\mathrm{h.c.}$"),
        (r"$\lambda_{12}^{(1)}$", r"$(S_1^\dagger S_1)(S_2^\dagger S_2)$"),
        (r"$\lambda_{12}^{(A)}$", r"$(S_1^\dagger T^AS_1)(S_2^\dagger T^AS_2)$"),
        (r"$\lambda_{12}^{(\times)}$", r"crossed independent $S_1$--$S_2$ contraction"),
        (r"$C_{LLS_1S_2}^{ij;AB}$", r"$L_iL_jS_1^AS_2^B$"),
        (r"$C_{ijab}$", r"$L_iL_j\phi_a\phi_b$"),
        (r"$C_5^{ij}$", r"$(L_i^TCL_j)HH$"),
        (r"$R_W$", r"$16\pi^2\beta_{C_5}\supset R_W\lambda_5C_{LLS_1S_2}$"),
    ])


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



def _stage_payload(record: RunRecord, stage_label: str) -> Mapping[str, object] | None:
    """Return one serialized EFT-stage payload by label."""
    for stage in record.summary.get("EFTStages", []) or []:
        if str(stage.get("Label", "")).strip() == stage_label:
            return stage
    return None


def _stage_labels(records: Sequence[RunRecord]) -> list[str]:
    """Return all EFT-stage labels in first-seen threshold order."""
    labels: list[str] = []
    for record in records:
        for stage in record.summary.get("EFTStages", []) or []:
            label = str(stage.get("Label", "")).strip()
            if label and label not in labels:
                labels.append(label)
    return labels


def _records_for_stage(
    records: Sequence[RunRecord],
    stage_label: str,
) -> list[RunRecord]:
    """Return successful records that actually contain ``stage_label``."""
    return [
        record
        for record in _valid_records(records)
        if _stage_payload(record, stage_label) is not None
    ]


def _stage_active_fields(
    record: RunRecord,
    stage_label: str,
) -> tuple[str, ...]:
    payload = _stage_payload(record, stage_label) or {}
    return tuple(str(field) for field in payload.get("ActiveHeavyFields", []) or [])


def _stage_integrated_fields(
    record: RunRecord,
    stage_label: str,
) -> tuple[str, ...]:
    payload = _stage_payload(record, stage_label) or {}
    return tuple(str(field) for field in payload.get("IntegratedFields", []) or [])


def _stage_is_f_first_eft1(
    record: RunRecord,
    stage_label: str,
) -> bool:
    """Return True for the analytically supported SM+S1+S2 EFT after F."""
    return (
        set(_stage_integrated_fields(record, stage_label)) == {"F"}
        and set(_stage_active_fields(record, stage_label))
        == ({"S"} if record.shared_scalar else {"S1", "S2"})
    )


def _stage_is_fully_decoupled(
    record: RunRecord,
    stage_label: str,
) -> bool:
    """Return True when no T3 heavy field remains dynamical."""
    return len(_stage_active_fields(record, stage_label)) == 0


def _stage_content_table(
    records: Sequence[RunRecord],
    stage_label: str,
) -> list[str]:
    """Show the threshold action and surviving heavy-field content."""
    if not records:
        return []

    shared = all(record.shared_scalar for record in records)
    if shared:
        lines = [
            r"\section*{Stage field content}",
            r"\begin{longtable}{@{}lcc p{0.24\linewidth}p{0.24\linewidth}@{}}",
            r"\toprule",
            r"run & $S$ & $F$ & integrated at this stage & active after matching \\",
            r"\midrule",
        ]
        for record in records:
            y_s = sp.Rational(1, 2)
            y_f = sp.Rational(0)
            integrated = _stage_integrated_fields(record, stage_label)
            active = _stage_active_fields(record, stage_label)
            lines.append(
                " & ".join(
                    [
                        latex_escape_text(record.name),
                        rf"$({record.d_s1},{_tex_expr(y_s)})$",
                        rf"$({record.d_f},{_tex_expr(y_f)})$",
                        latex_escape_text(", ".join(integrated) if integrated else "none"),
                        latex_escape_text(", ".join(active) if active else "none"),
                    ]
                )
                + r" \\"
            )
        lines.extend([r"\bottomrule", r"\end{longtable}"])
        return lines

    lines = [
        r"\section*{Stage field content}",
        r"\begin{longtable}{@{}lccc p{0.20\linewidth}p{0.20\linewidth}@{}}",
        r"\toprule",
        r"run & $S_1$ & $S_2$ & $F$ & integrated at this stage & active after matching \\",
        r"\midrule",
    ]

    for record in records:
        y_s1, y_s2, y_f = _hypercharges(record)
        representations = {
            "S1": rf"$({record.d_s1},{_tex_expr(y_s1)})$",
            "S2": rf"$({record.d_s2},{_tex_expr(y_s2)})$",
            "F": rf"$({record.d_f},{_tex_expr(y_f)})$",
        }
        integrated = _stage_integrated_fields(record, stage_label)
        active = _stage_active_fields(record, stage_label)
        lines.append(
            " & ".join(
                [
                    latex_escape_text(record.name),
                    representations["S1"],
                    representations["S2"],
                    representations["F"],
                    latex_escape_text(", ".join(integrated) if integrated else "none"),
                    latex_escape_text(", ".join(active) if active else "none"),
                ]
            )
            + r" \\"
        )

    lines.extend([r"\bottomrule", r"\end{longtable}"])
    return lines


def _active_quantity_table(
    records: Sequence[RunRecord],
    stage_label: str,
) -> list[str]:
    """List which T3 renormalisable quantities remain active in a generic stage.

    This table deliberately does not recycle UV beta-function coefficients after
    a field has been removed.  The coefficient set changes with the EFT field
    content and must come from a dedicated stage-RGE calculation.
    """
    if not records:
        return []

    quantity_rules = (
        (r"$y_1$", lambda active: {"F", "S1"}.issubset(active)),
        (r"$y_2$", lambda active: {"F", "S2"}.issubset(active)),
        (r"$M_F$", lambda active: "F" in active),
        (r"$m_1^2$", lambda active: "S1" in active),
        (r"$m_2^2$", lambda active: "S2" in active),
        (r"$\lambda_{S_1}^{(1)}$", lambda active: "S1" in active),
        (r"$\lambda_{S_2}^{(1)}$", lambda active: "S2" in active),
        (r"$\lambda_{HS_1}^{(1)},\,\lambda_{HS_1}^{(A)}$", lambda active: "S1" in active),
        (r"$\lambda_{HS_2}^{(1)},\,\lambda_{HS_2}^{(A)}$", lambda active: "S2" in active),
        (
            r"$\lambda_5,\,\lambda_{12}^{(1)},\,\lambda_{12}^{(A)},\,\lambda_{12}^{(\times)}$",
            lambda active: {"S1", "S2"}.issubset(active),
        ),
    )

    lines = [
        r"\section*{Active renormalisable T3 quantities}",
        r"$g_1$, $g_2$, $g_3$, $Y_u$, $Y_d$, $Y_e$ and the SM Higgs "
        r"couplings remain active at every stage.",
        r"\begin{longtable}{@{}p{0.45\linewidth}"
        + "c" * len(records)
        + r"@{}}",
        r"\toprule",
        " & ".join(
            [r"quantity"]
            + [rf"$ {_model_math_label(record)} $" for record in records]
        )
        + r" \\",
        r"\midrule",
    ]

    for symbol, predicate in quantity_rules:
        cells = [symbol]
        for record in records:
            active = set(_stage_active_fields(record, stage_label))
            cells.append(r"$\checkmark$" if predicate(active) else r"$-$")
        lines.append(" & ".join(cells) + r" \\")

    lines.extend([r"\bottomrule", r"\end{longtable}"])
    return lines


def _smeft_notation_key() -> list[str]:
    """Minimal notation key for a stage with only the SMEFT Weinberg operator."""
    return [
        r"\section*{Notation key}",
        r"\begin{longtable}{@{}p{0.19\linewidth}p{0.73\linewidth}@{}}",
        r"\toprule",
        r"symbol & interaction / definition \\",
        r"\midrule",
        r"$g_1$ & $D_\mu\supset i g_1 Y B_\mu$ \\",
        r"$g_2$ & $D_\mu\supset i g_2 T^A W_\mu^A$ \\",
        r"$g_3$ & $D_\mu\supset i g_3 t^A G_\mu^A$ \\",
        r"$Y_u$ & $-\bar Q\,Y_u\,\widetilde H\,u+\mathrm{h.c.}$ \\",
        r"$Y_d$ & $-\bar Q\,Y_d\,H\,d+\mathrm{h.c.}$ \\",
        r"$Y_e$ & $-\bar L\,Y_e\,H\,e+\mathrm{h.c.}$ \\",
        r"$T$ & $\operatorname{Tr}\!\left(Y_e^\dagger Y_e"
        r"+3Y_u^\dagger Y_u+3Y_d^\dagger Y_d\right)$ \\",
        r"$\lambda_1$ & $\frac12\lambda_1(H^\dagger H)^2$ \\",
        r"$C_5^{ij}$ & $(L_i^T C L_j)\,H H$ \\",
        r"\bottomrule",
        r"\end{longtable}",
    ]


def _smeft_c5_group_factor_table(records: Sequence[RunRecord]) -> list[str]:
    """Universal one-loop SMEFT Weinberg group factors after full decoupling."""
    per_run = {
        _run_key(record): {
            "lambdaH*C5": sp.Integer(2),
            "g2_sq*C5": sp.Integer(-3),
            "T*C5": sp.Integer(2),
            "YeLeft": sp.Rational(-3, 2),
            "YeRight": sp.Rational(-3, 2),
        }
        for record in records
    }

    previous = dict(_SPECIAL_STRUCTURE_TEX)
    _SPECIAL_STRUCTURE_TEX.update(
        {
            "YeLeft": r"(Y_e^\dagger Y_e)^T C_5",
            "YeRight": r"C_5(Y_e^\dagger Y_e)",
        }
    )
    try:
        return _comparison_table(
            r"\beta_{C_5}",
            records,
            per_run,
        )
    finally:
        _SPECIAL_STRUCTURE_TEX.clear()
        _SPECIAL_STRUCTURE_TEX.update(previous)


def _write_gf_f_first_stage_shared(
    records: list[RunRecord],
    *,
    stage_label: str,
    report_root: Path,
) -> Path:
    path = group_factor_report_path(stage_label, report_root=report_root)
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = _document_header(
        "Scotogenic/shared-scalar group factors: " + stage_label.replace("_", r"\_"),
        r"Active theory: $\mathrm{SM}+S+C_{LLSS}$. The heavy fermion is absent.",
    )
    lines.extend(_shared_scalar_notation_key())
    lines.extend(_stage_content_table(records, stage_label))
    lines.extend(_saved_rge_beta_sections(
        records,
        loader=load_eft1_renormalisable_rge_payload,
        couplings=("gY", "g2", "g3", "yu", "yd", "ye", "mSSq",
                   "lambdaH", "lambdaS", "lambda3", "lambda4", "lambda5"),
        title="Renormalisable one-loop beta functions in this EFT",
    ))
    lines.extend(_direct_weinberg_comparison(records))
    lines.extend([
        r"\section*{Matching-basis note}",
        r"The Matchete threshold calculation still uses the two formal T3 scalar "
        r"legs internally.  In this report they are identified with one physical "
        r"field through $S_1=i\sigma_2S^*$ and $S_2=S$; raw formal-basis "
        r"$y_1/y_2$ component samples are therefore not displayed.",
    ])
    lines.extend([r"\end{document}", ""])
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def _write_gf_f_first_stage(
    records: list[RunRecord],
    *,
    stage_label: str,
    report_root: Path,
) -> Path:
    """Write the full analytically supported F-first intermediate-EFT report."""
    if records and all(record.shared_scalar for record in records):
        return _write_gf_f_first_stage_shared(
            records, stage_label=stage_label, report_root=report_root
        )
    path = group_factor_report_path(stage_label, report_root=report_root)
    path.parent.mkdir(parents=True, exist_ok=True)

    lines = _document_header(
        "T3 group factors: " + stage_label.replace("_", r"\_"),
        (
            r"Active theory: $\mathrm{SM}+S_1+S_2+C_{LLS_1S_2}$.  "
            r"The fermion $F$ is absent, so $y_1$, $y_2$ and $M_F$ are not "
            r"running EFT couplings."
        ),
    )
    lines.extend(_uv_notation_key())
    lines.extend(_stage_content_table(records, stage_label))
    lines.extend(_representation_table(records, include_f=False))
    lines.extend(
        _saved_rge_beta_sections(
            records,
            loader=load_eft1_renormalisable_rge_payload,
            couplings=("gY", "g2", "g3", "yu", "yd", "ye"),
            title="Gauge and Yukawa beta functions in this EFT",
        )
    )

    eft_data = {
        _run_key(record): _eft1_group_factor_data(record)
        for record in records
    }

    for key in ("mS1Sq", "mS2Sq", "lambdaH1", "lambdaH2", "lambda12", "lambdaT3"):
        lines.extend(
            _comparison_table(
                _beta_name_tex(key),
                records,
                {
                    _run_key(record): eft_data[_run_key(record)][key]
                    for record in records
                },
            )
        )

    lines.extend(_generated_non_singlet_section(records, eft_data))
    lines.extend(_direct_weinberg_comparison(records))
    lines.extend(_matching_tensor_appendix(records))

    lines.extend(
        [
            r"\section*{Next threshold}",
            r"Further matching depends on the next field or fields selected by "
            r"the threshold plan.  Once all T3 fields are removed, their "
            r"representation dependence is carried by the matched $C_5$ "
            r"boundary condition and the subsequent SMEFT running is universal.",
            r"\end{document}",
            "",
        ]
    )

    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def _write_gf_fully_decoupled_stage(
    records: list[RunRecord],
    *,
    stage_label: str,
    report_root: Path,
) -> Path:
    """Write the universal SMEFT GroupFactors report after all T3 fields decouple."""
    path = group_factor_report_path(stage_label, report_root=report_root)
    path.parent.mkdir(parents=True, exist_ok=True)

    lines = _document_header(
        "T3 group factors: " + stage_label.replace("_", r"\_"),
        (
            r"Active theory: SMEFT with the Weinberg operator $C_5$.  "
            r"No T3 heavy field remains dynamical."
        ),
    )
    lines.extend(_smeft_notation_key())
    lines.extend(_stage_content_table(records, stage_label))
    if records and all(record.shared_scalar for record in records):
        lines.extend(
            [
                r"\section*{Decoupling statement}",
                r"All dependence on the heavy scalar $S$ and fermion $F$ is "
                r"contained in the matched boundary value "
                r"$C_5(\mu_{\mathrm{th}})$.  Below the final threshold the "
                r"one-loop Weinberg-operator running is the universal SMEFT result.",
            ]
        )
    else:
        lines.extend(
            [
                r"\section*{Decoupling statement}",
                r"All dependence on the original $S_1$, $S_2$ and $F$ "
                r"representations is contained in the matched boundary value "
                r"$C_5(\mu_{\mathrm{th}})$.  The one-loop running below the final "
                r"threshold is common to all T3 classes.",
            ]
        )
    lines.extend(_smeft_c5_group_factor_table(records))
    lines.extend(
        [
            r"\section*{Flavor form}",
            r"\begin{align}",
            r"16\pi^2\frac{dC_5}{d\ln\mu}"
            r"&=(2\lambda_1-3g_2^2+2T)C_5\\",
            r"&\quad-\frac32\left["
            r"(Y_e^\dagger Y_e)^T C_5+C_5(Y_e^\dagger Y_e)"
            r"\right].",
            r"\end{align}",
            r"\end{document}",
            "",
        ]
    )

    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def _write_gf_generic_intermediate_stage(
    records: list[RunRecord],
    *,
    stage_label: str,
    report_root: Path,
) -> Path:
    """Write a safe GroupFactors report for any other threshold ordering.

    A dedicated RGE coefficient table is intentionally not fabricated here:
    removing S1 or S2 before F changes the dynamical field content, so the UV
    or F-first EFT beta coefficients cannot simply be recycled.
    """
    path = group_factor_report_path(stage_label, report_root=report_root)
    path.parent.mkdir(parents=True, exist_ok=True)

    active_examples = sorted(
        {
            tuple(_stage_active_fields(record, stage_label))
            for record in records
        }
    )
    active_text = "; ".join(
        ", ".join(active) if active else "none"
        for active in active_examples
    )

    lines = _document_header(
        "T3 group factors: " + stage_label.replace("_", r"\_"),
        (
            r"Threshold-order-aware EFT stage.  Active heavy-field content: "
            + latex_escape_text(active_text)
            + r"."
        ),
    )
    lines.extend(_eft1_notation_key())
    lines.extend(_stage_content_table(records, stage_label))
    lines.extend(_active_quantity_table(records, stage_label))
    lines.extend(
        [
            r"\section*{Stage-specific group-factor status}",
            r"The threshold ordering at this stage is supported by the report "
            r"layer, so a PDF is produced for it.  A dedicated analytic "
            r"renormalisable-RGE group-factor module is not currently available "
            r"for this reduced heavy-field content.  Coefficients from the UV "
            r"theory or from the special $F$-first EFT are therefore not reused "
            r"here, because doing so would assign beta-function coefficients "
            r"from the wrong dynamical theory.",
            r"\section*{Representation invariants}",
        ]
    )
    lines.extend(_representation_table(records, include_f=True))
    lines.extend(
        [
            r"\section*{Next threshold}",
            r"The next GroupFactors PDF is generated from the next serialized "
            r"EFT-stage label in the threshold plan, independent of which field "
            r"is removed next.",
            r"\end{document}",
            "",
        ]
    )

    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def write_gf_stage(
    records: list[RunRecord],
    *,
    stage_label: str,
    report_root: Path,
) -> Path | None:
    """Write the GroupFactors report for one actual serialized EFT stage."""
    relevant = _records_for_stage(records, stage_label)
    if not relevant:
        return None

    if all(_stage_is_f_first_eft1(record, stage_label) for record in relevant):
        return _write_gf_f_first_stage(
            relevant,
            stage_label=stage_label,
            report_root=report_root,
        )

    if all(_stage_is_fully_decoupled(record, stage_label) for record in relevant):
        return _write_gf_fully_decoupled_stage(
            relevant,
            stage_label=stage_label,
            report_root=report_root,
        )

    return _write_gf_generic_intermediate_stage(
        relevant,
        stage_label=stage_label,
        report_root=report_root,
    )


def _shared_scalar_notation_key() -> list[str]:
    return [
        r"\section*{Notation key}",
        r"\begin{longtable}{@{}p{0.21\linewidth}p{0.71\linewidth}@{}}",
        r"\toprule",
        r"symbol & interaction / definition \\",
        r"\midrule",
        r"$S$ & physical inert scalar doublet, $S\sim(2,+1/2)$ \\",
        r"$\widetilde S$ & $i\sigma_2S^*$; formal matching leg $S_1=\widetilde S$, while $S_2=S$ \\",
        r"$F$ & neutral singlet or triplet Majorana fermion \\",
        r"$h$ & $LFS+\mathrm{h.c.}$ \\",
        r"$M_F$ & heavy-fermion Majorana mass \\",
        r"$T_\nu$ & $\operatorname{Tr}(h^\dagger h)$ \\",
        r"$m_S^2$ & $m_S^2 S^\dagger S$ \\",
        r"$\lambda_1$ & $\frac12\lambda_1(H^\dagger H)^2$ \\",
        r"$\lambda_2$ & $\frac12\lambda_2(S^\dagger S)^2$ \\",
        r"$\lambda_3$ & $\lambda_3(H^\dagger H)(S^\dagger S)$ \\",
        r"$\lambda_4$ & $\lambda_4(H^\dagger S)(S^\dagger H)$ \\",
        r"$\lambda_5$ & $\frac12\lambda_5[(H^\dagger S)^2+\mathrm{h.c.}]$ \\",
        r"\bottomrule",
        r"\end{longtable}",
    ]


def _write_gf_uv_shared(records: list[RunRecord], *, report_root: Path) -> Path:
    path = group_factor_report_path("UV", report_root=report_root)
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = _document_header(
        "Scotogenic/shared-scalar group factors: UV",
        (
            r"Active theory: $\mathrm{SM}+S+F$.  The formal T3 legs are "
            r"identified as $S_1=i\sigma_2S^*$ and $S_2=S$, so the physical "
            r"scalar is counted once in loops."
        ),
    )
    lines.extend(_shared_scalar_notation_key())
    lines.extend(_saved_rge_beta_sections(
        records,
        loader=load_uv_rge_payload,
        couplings=("gY", "g2", "g3", "yu", "yd", "ye", "h", "MF", "mSSq",
                   "lambdaH", "lambdaS", "lambda3", "lambda4", "lambda5"),
        title="One-loop beta functions in the physical one-scalar theory",
    ))
    lines.extend([
        r"\section*{Formal matching bridge}",
        r"Matchete keeps the validated two-leg T3 topology, but those legs are two descriptions of the same physical scalar in this branch.",
        r"\end{document}", "",
    ])
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def write_gf_uv(records: list[RunRecord], *, report_root: Path) -> Path:
    valid = _valid_records(records)
    if valid and all(record.shared_scalar for record in valid):
        return _write_gf_uv_shared(valid, report_root=report_root)
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
    lines.extend(paper_notation_key_lines())
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
            r"\section*{First threshold}",
            r"The first EFT GroupFactors report below is selected from the "
            r"actual threshold plan stored in the run.  No assumption is made "
            r"here about whether $F$, $S_1$, $S_2$, or a degenerate set is "
            r"integrated out first.",
            r"\end{document}",
            "",
        ]
    )

    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def write_gf_eft1(
    records: list[RunRecord],
    *,
    report_root: Path,
) -> Path | None:
    """Compatibility wrapper for the special F-first intermediate EFT."""
    for stage_label in _stage_labels(_valid_records(records)):
        relevant = _records_for_stage(records, stage_label)
        if relevant and all(
            _stage_is_f_first_eft1(record, stage_label)
            for record in relevant
        ):
            return _write_gf_f_first_stage(
                relevant,
                stage_label=stage_label,
                report_root=report_root,
            )
    return None


def write_and_compile_stage_group_factor_reports(
    records: list[RunRecord],
    *,
    report_root: Path,
) -> list[Path]:
    """Write and compile UV plus every EFT-stage GroupFactors report.

    Stage labels are taken directly from ``summary["EFTStages"]``.  Therefore
    the report set follows the threshold plan actually used in the run:
    F-first, S1-first, S2-first, simultaneous thresholds, or longer sequential
    chains all receive one PDF per recorded EFT stage.
    """
    valid = _valid_records(records)
    outputs: list[Path] = []

    uv = write_gf_uv(valid, report_root=report_root)
    compile_latex_document(uv)
    outputs.append(uv)

    for stage_label in _stage_labels(valid):
        stage_report = write_gf_stage(
            valid,
            stage_label=stage_label,
            report_root=report_root,
        )
        if stage_report is None:
            continue
        compile_latex_document(stage_report)
        outputs.append(stage_report)

    return outputs

