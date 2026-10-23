"""Insert representation-specific C12 factors and C12/C5 mixing into GroupFactors PDF.

Uses the actual Matchete-CG tensor and saved one-loop Wilson tensor beta.
Does not regenerate EFT physics or create a separate report family.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable, Sequence
from collections import defaultdict

import sympy as sp

from RGE.group_factors.validation.C12RepresentationGroupFactors import derive_from_saved_data

BEGIN = "% BEGIN REPRESENTATION-SPECIFIC C12 GROUP FACTORS (generated)"
END = "% END REPRESENTATION-SPECIFIC C12 GROUP FACTORS (generated)"
INPUT_FILES = (
    "eft1_after_F_wilson_seed.json",
    "eft1_rgbeta_rge.json",
    "eft1_wilson_rge.json",
)
SECTORS = (
    ("yukawa_wavefunction", r"Yukawa wavefunction"),
    ("scalar_pair", r"Scalar-pair gauge/quartic"),
    ("mixed_yukawa_gauge", r"Mixed Yukawa and gauge"),
    ("crossed_yukawa", r"Crossed Yukawa"),
    ("conjugate_coefficient_yukawa", r"Conjugate Wilson/Yukawa"),
    ("scalar_anomalous_dimension", r"Scalar anomalous dimension"),
    ("fermion_anomalous_dimension", r"Fermion anomalous dimension"),
    ("total", r"Total beta on $C_{12}$ seed support"),
)


def _latex_escape(value: str) -> str:
    """Escape text (not the explicitly authored math in this module)."""
    table = {"\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "$": r"\$",
             "#": r"\#", "_": r"\_", "{": r"\{", "}": r"\}",
             "~": r"\textasciitilde{}", "^": r"\textasciicircum{}"}
    return "".join(table.get(c, c) for c in str(value))


def _inputs_exist(data_dir: Path) -> bool:
    return all((Path(data_dir) / name).is_file() for name in INPUT_FILES)


def _coefficient_latex(value: str) -> str:
    """Render the exact saved SymPy expression in readable RGBeta notation."""
    expression = sp.sympify(value, locals={"sqrt": sp.sqrt, "conjugate": sp.conjugate})
    names = {
        "g1": r"g_1", "g2": r"g_2", "g3": r"g_3",
        "ye": r"y_e", "yu": r"y_u", "yd": r"y_d",
        "lambdaH": r"\lambda_H", "lambdaT3": r"\lambda_{T3}",
        "lambda12": r"\lambda_{12}",
        "lambda12Adj": r"\lambda_{12,\mathrm{Adj}}",
        "lambda12Inv1": r"\lambda_{12,\mathrm{Inv1}}",
        "lambda12Inv2": r"\lambda_{12,\mathrm{Inv2}}",
        "lambdaS1Inv1": r"\lambda_{S_1,\mathrm{Inv1}}",
        "lambdaS2Inv1": r"\lambda_{S_2,\mathrm{Inv1}}",
        "MF": r"M_F",
    }
    return sp.latex(expression, symbol_names={
        symbol: names.get(symbol.name, symbol.name)
        for symbol in expression.free_symbols
    })


def _model_class(label: str, rep: dict) -> str:
    """Use the named T3 class where available; do not infer physics from it."""
    import re
    match = re.search(r"T3[-_]([A-E])(?:-|_|\b)", label, re.I)
    if match:
        return match.group(1).upper()
    # The representation is the authoritative metadata, class is presentation.
    ds1, ds2, df = (int(rep[k]) for k in ("dS1", "dS2", "dF"))
    return {(1, 3, 2): "A", (2, 2, 1): "B", (2, 2, 3): "C",
            (3, 1, 2): "D", (3, 3, 2): "E"}.get((ds1, ds2, df), "OTHER")


def _model_header(label: str, rep: dict) -> str:
    import re
    cls = _model_class(label, rep)
    if cls not in "ABCDE" or len(cls) != 1:
        return _latex_escape(label)
    a = int(rep["alpha"])
    suffix = f"m{-a}" if a < 0 else f"p{a}"
    # Keep alpha visible, as in existing C5 block tables.
    return rf"T3-{cls}-{suffix} ($\alpha={a}$)"


def _factor_terms(expression: str) -> dict[str, sp.Expr]:
    """Split an exact polynomial into symbolic structures and numeric factors.

    Uses the full projected beta, NOT the isolated master-tensor terms: gauge
    pieces can be non-proportional individually yet sum to a collinear total.
    """
    expr = sp.expand(sp.sympify(expression, locals={"sqrt": sp.sqrt, "conjugate": sp.conjugate}))
    factors: dict[str, sp.Expr] = {}
    if expr == 0:
        return factors
    for term in sp.Add.make_args(expr):
        coupling_symbols = sorted(term.free_symbols, key=lambda x: str(x))
        coefficient = sp.simplify(term.subs({symbol: 1 for symbol in coupling_symbols}))
        if coefficient == 0 or coefficient.free_symbols:
            raise ValueError(f"Non-numeric group factor in projected coefficient: {term}")
        structure = sp.simplify(term / coefficient)
        key = sp.sstr(structure)
        factors[key] = sp.simplify(factors.get(key, 0) + coefficient)
    return {name: value for name, value in factors.items() if value != 0}


def _group_structure_latex(raw: str) -> str:
    expr = sp.sympify(raw, locals={"sqrt": sp.sqrt, "conjugate": sp.conjugate})
    latex = _coefficient_latex(sp.sstr(expr))
    if expr == 1:
        return r"C_{12}"
    return latex + r"\,C_{12}"


def _sort_structure(name: str) -> tuple[int, str]:
    # Follow gauge, Yukawa, and quartic ordering used by the existing reports.
    if name.startswith("g1"):
        return 0, name
    if name.startswith("g2"):
        return 1, name
    if name.startswith("g3"):
        return 2, name
    if "ye" in name or "yu" in name or "yd" in name:
        return 3, name
    if "lambda" in name:
        return 4, name
    return 5, name


def _run_blocks(results: list[dict]) -> list[list[dict]]:
    buckets: dict[str, list[dict]] = defaultdict(list)
    for item in results:
        cls = _model_class(item["model"], item["representation"])
        buckets[cls].append(item)
    blocks = []
    for cls in ("A", "B", "C", "D", "E", "OTHER"):
        group = sorted(buckets.get(cls, []), key=lambda r: int(r["representation"]["alpha"]))
        for i in range(0, len(group), 3):
            blocks.append(group[i:i + 3])
    return blocks


def _comparison_table(block: list[dict], idx: int) -> list[str]:
    heads = [_model_header(item["model"], item["representation"]) for item in block]
    n = len(block)
    title = _latex_escape(block[0]["model"] + " -- " + block[-1]["model"])
    colspec = "p{0.26\\linewidth}" + "".join([rf"p{{{0.68/n:.4f}\linewidth}}" for _ in block])
    rows = sorted({key for item in block for key in item["terms"]}, key=_sort_structure)
    lines = [
        rf"\paragraph*{{Run block {idx}: {title}}}",
        rf"\begin{{longtable}}{{@{{}}{colspec}@{{}}}}",
        r"\toprule",
        r"$S$ & " + " & ".join(heads) + r" \\",
        r"\midrule",
    ]
    for key in rows:
        cells = []
        for item in block:
            if item["status"] == "MissingSavedTensorRGE":
                cells.append(r"\textit{no data}")
            else:
                cells.append("$" + sp.latex(item["terms"].get(key, sp.S.Zero)) + "$")
        lines.append("$" + _group_structure_latex(key) + "$ & " + " & ".join(cells) + r" \\")
    # Keep one structure-by-model table per run block (as in C5), with a
    # visibly separate subsection for the off-diagonal anomalous dimensions.
    lines.append(r"\midrule")
    lines.append(
        rf"\multicolumn{{{n+1}}}{{@{{}}l}}{{\textit{{Off-diagonal "
        r"$C_{12}\leftrightarrow C_5$ mixing}} \\"
    )
    for structure, factor_fn in (
        (r"$\lambda_5 C_{12}\to\beta_{C_5}$", _forward_mixing_factor),
        (r"$\lambda_5^* C_5\to\beta_{C_{12}}$", _reverse_mixing_factor),
    ):
        cells = []
        for item in block:
            factor = factor_fn(item)
            cells.append(
                rf"${sp.latex(factor)}$" if factor is not None
                else r"\textit{not verified}"
            )
        lines.append(structure + " & " + " & ".join(cells) + r" \\")
    lines.extend([r"\bottomrule", r"\end{longtable}"])
    # Non-collinearity is an operator-space diagnostic, not a group factor.
    # Explain it outside the numerical table rather than as beta_perp != 0.
    not_closed = [x["model"] for x in block
                  if x["status"] == "Projected" and not x["total_collinear"]]
    if not_closed:
        names = ", ".join(_latex_escape(x) for x in not_closed)
        lines.extend([
            r"\noindent\textit{Projection status:} " + names + r" are not collinear",
            r"with their matched $C_{12}$ tensor on seed support. The entries",
            r"above then show only individually aligned source terms, not a",
            r"complete closed $C_{12}$ RGE.",
            r"\par\medskip",
        ])
    return lines


def _forward_mixing_factor(item: dict) -> sp.Expr | None:
    """Existing analytic R_W for C12 -> C5; does not rederive that factor."""
    if item["status"] != "Projected":
        return None
    from RGE.group_factors.validation.ValidateDirectWeinbergGroupFactors import (
        direct_weinberg_group_factor,
    )
    rep = item["representation"]
    factor = direct_weinberg_group_factor(
        int(rep["dS1"]), int(rep["dS2"]), int(rep["dF"]),
    )
    return sp.sympify(factor.reduced_factor, locals={"sqrt": sp.sqrt})


def _reverse_mixing_factor(item: dict) -> sp.Expr | None:
    """T3-B alpha=-1 reverse coefficient previously checked in tensor RGE.

    Other representations have not passed a comparable C5-seeded projection.
    Missing results are *unknown*, not zero.
    """
    if item["status"] != "Projected":
        return None
    rep = item["representation"]
    dimensions = tuple(int(rep[key]) for key in ("dS1", "dS2", "dF", "alpha"))
    if dimensions == (2, 2, 1, -1) and not rep.get("shared_scalar", False):
        return 8 * sp.sqrt(3) / 3
    return None


def render_c12_pdf_section(model_data: Sequence[tuple[str, Path]]) -> tuple[str, list[dict]]:
    """Create C5-style structure-by-model tables from saved exact C12 projections."""
    if not model_data:
        raise ValueError("No model data supplied for representation-specific C12 report")
    lines = [
        BEGIN,
        r"\clearpage",
        r"\section*{Representation-specific $C_{12}$ one-loop group factors}",
        r"The table layout matches the $C_5$ group-factor comparison: rows are",
        r"independent symbolic RGE structures, while columns are model representations.",
        r"Entries are exact coefficients in $16\pi^2\beta_{C_{12}}$,",
        r"derived from each model's saved Matchete-CG Wilson seed and one-loop",
        r"tensor RGE. All heavy fermions are integrated out at one common scale.",
        r"A table entry is a coefficient of the aligned $C_{12}$ tensor component.",
        r"Non-collinear operator mixing is recorded separately as a closure warning.",
    ]
    items: list[dict] = []
    for label, raw_dir in model_data:
        data_dir = Path(raw_dir)
        missing = [name for name in INPUT_FILES if not (data_dir / name).is_file()]
        if missing:
            import re
            match = re.search(r"T3[-_]([A-E])[-_]([mp])(\d+)", label)
            model_alpha = (-1 if match and match.group(2) == 'm' else 1) * int(match.group(3)) if match else 0
            dims = {"A": (1, 3, 2), "B": (2, 2, 1), "C": (2, 2, 3), "D": (3, 1, 2), "E": (3, 3, 2)}
            ds = dims.get(match.group(1), (0, 0, 0)) if match else (0, 0, 0)
            rep = {"dS1": ds[0], "dS2": ds[1], "dF": ds[2], "alpha": model_alpha}
            items.append({"model": label, "status": "MissingSavedTensorRGE",
                          "representation": rep, "terms": {}, "total_collinear": False,
                          "missing": missing})
            continue
        result = derive_from_saved_data(data_dir)
        rep = result["representation"]
        total = result["sectors"]["total"]
        collinear = bool(total["collinear_on_c12_seed_support"])
        if collinear:
            terms = _factor_terms(str(total["coefficient"]))
        else:
            # Display only those individual sources which ARE proportional to
            # the matched tensor, as a partial result, with beta_perp marked.
            # Sum before splitting, to handle cancellations between sectors.
            aligned = sp.S.Zero
            for name, source in result["sectors"].items():
                if name == "total" or not source["collinear_on_c12_seed_support"]:
                    continue
                aligned += sp.sympify(str(source["coefficient"]), locals={
                    "sqrt": sp.sqrt, "conjugate": sp.conjugate})
            terms = _factor_terms(sp.sstr(sp.expand(aligned)))
        items.append({"model": label, "status": "Projected",
                      "representation": rep, "terms": terms,
                      "total_collinear": collinear,
                      "seed_support_component_count": result["seed_support_component_count"],
                      "generated_outside_support": len(result["generated_mixed_components_outside_c12_seed_support"])})
    usable = [item for item in items if item["status"] == "Projected"]
    if not usable:
        raise ValueError("No C12 group factors could be derived: saved tensor data missing for every model")
    for idx, block in enumerate(_run_blocks(items), 1):
        lines.extend(_comparison_table(block, idx))
    lines.extend([
        r"\paragraph{Mixing conventions.} All factors are defined in",
        r"$16\pi^2\,dC/d\ln\mu$. The $C_{12}\to C_5$ row displays",
        r"the existing analytic $R_W$ factor (repeated from the Weinberg section);",
        r"its representation-dependent phase/normalisation validation remains",
        r"a separate issue. The $C_5\to C_{12}$ row uses the independently",
        r"checked T3-B, $\alpha=-1$ coefficient in the retained",
        r"ordered-tensor two-operator convention. This is not an additional",
        r"contribution to the earlier Weinberg table.",
        r"\textit{Not verified} means no validated reverse projection is",
        r"available for that representation; it does not imply vanishing mixing.",
        r"In the numerical neutrino-mass construction a convention-dependent",
        r"factor of two relates the ordered $C_5$ tensor to the physical",
        r"Majorana coefficient.",
    ])
    for item in items:
        if item["status"] == "MissingSavedTensorRGE":
            lines.append(r"\noindent " + _latex_escape(item["model"])
                         + r": saved tensor inputs unavailable (" + _latex_escape(", ".join(item["missing"])) + r").")
        elif not item["total_collinear"]:
            pass  # Already explained below the relevant comparison table.
    lines.extend([
        r"\paragraph{Scope.} The coefficients are obtained by extracting the",
        r"$y_1^*y_2^*$ channel of the saved one-generation scalar-only EFT tensor",
        r"RGE. They do not prove full three-generation operator closure.",
        r"The SM-like $C_5$ tables elsewhere in this report are unchanged.",
        END,
    ])
    details = [{k: v for k, v in item.items() if k != "terms"} for item in items]
    return "\n".join(lines) + "\n", details

def replace_pdf_tex_section(text: str, generated: str) -> str:
    """Idempotently replace only our marked section immediately before end document."""
    if text.count(BEGIN) != text.count(END):
        raise ValueError("Mismatched C12 PDF report markers")
    if text.count(BEGIN) > 1:
        raise ValueError("Duplicate C12 PDF report markers")
    if BEGIN in text:
        start = text.index(BEGIN)
        end = text.index(END) + len(END)
        if end <= start:
            raise ValueError("Malformed C12 PDF report markers")
        return text[:start] + generated.rstrip("\n") + text[end:]
    trailer = r"\end{document}"
    if text.count(trailer) != 1:
        raise ValueError("GroupFactors TeX requires exactly one end-of-document marker")
    return text.replace(trailer, generated + trailer, 1)


def update_existing_group_factor_pdf(*, tex_path: Path,
                                     model_data: Sequence[tuple[str, Path]],
                                     compile_pdf: bool = True) -> dict:
    """Update existing GroupFactors TeX in place, then rebuild SAME-STEM PDF."""
    path = Path(tex_path)
    if not path.is_file():
        raise FileNotFoundError(f"Existing GroupFactors TeX not found: {path}")
    generated, details = render_c12_pdf_section(model_data)
    original = path.read_text(encoding="utf-8-sig")
    updated = replace_pdf_tex_section(original, generated)
    changed = updated != original
    if changed:
        path.write_text(updated, encoding="utf-8")
    if compile_pdf and (changed or not path.with_suffix(".pdf").is_file()):
        from Reports.ReportGeneration import compile_latex_document
        compile_latex_document(path)
        if not path.with_suffix(".pdf").is_file():
            raise RuntimeError(f"PDF compilation did not produce {path.with_suffix('.pdf')}")
    return {"tex": str(path), "pdf": str(path.with_suffix(".pdf")),
            "changed": changed, "models": details}


def refresh_stage_group_factor_pdfs(records: Iterable[object],
                                   tex_files: Iterable[Path]) -> list[dict]:
    """Called after the normal stage PDF generator, using actual F-first records."""
    from Reports.GroupFactorReports import _records_at_stage, _stage_is_f_first_eft1

    rows = list(records)
    result = []
    for tex in tex_files:
        path = Path(tex)
        if not path.name.startswith("GF_") or path.suffix.lower() != ".tex":
            continue
        label = path.stem[3:]
        relevant = _records_at_stage(rows, label)
        relevant = [r for r in relevant if _stage_is_f_first_eft1(r, label) and not r.shared_scalar]
        if not relevant:
            continue
        model_data = [(r.name, Path(r.output_dir) / "data") for r in relevant]
        result.append(update_existing_group_factor_pdf(tex_path=path, model_data=model_data))
    return result


def _discover_models(source: Path) -> list[tuple[str, Path]]:
    """Recover all saved representations from the usual analytical study data."""
    from scripts.BuildCombinedAnalyticalReports import reconstruct_records
    records = reconstruct_records(Path(source))
    return [(record.name, Path(record.output_dir) / "data") for record in records]


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tex", type=Path, default=Path("Reports/output/analytical/GroupFactors/GF_EFT_1_after_F.tex"))
    parser.add_argument("--data-dir", action="append", type=Path,
                        help="Saved model data directory; repeat for each model")
    parser.add_argument("--label", action="append", help="Model name corresponding to each --data-dir")
    parser.add_argument("--source", type=Path,
                        default=Path("output/full/comparison/smallY_smallL"),
                        help="Saved full-study scenario; used automatically if no --data-dir supplied")
    parser.add_argument("--compile", action="store_true", help="Compile PDF; default also compiles")
    parser.add_argument("--no-compile", action="store_true")
    args = parser.parse_args(argv)
    if args.compile and args.no_compile:
        parser.error("--compile and --no-compile are mutually exclusive")
    if args.data_dir:
        labels = args.label or []
        if labels and len(labels) != len(args.data_dir):
            parser.error("One --label required per --data-dir")
        models = [(labels[i] if labels else data.parent.name, data)
                  for i, data in enumerate(args.data_dir)]
    else:
        if args.label:
            parser.error("--label requires --data-dir")
        models = _discover_models(args.source)
    result = update_existing_group_factor_pdf(
        tex_path=args.tex, model_data=models, compile_pdf=not args.no_compile,
    )
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
