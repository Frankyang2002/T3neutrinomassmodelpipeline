from __future__ import annotations

"""Dependency-free beta_y2 nongauge diagnostic.

Reads saved RGBeta UV JSON files and checks whether the nongauge part of the
reported beta_y2 expression is alpha-independent for each T3 representation.

No project-internal report/parser module is imported, so this script can be run
directly from the repository root regardless of where RGEComparison.py lives.
"""

import argparse
import json
from pathlib import Path
import re
from fractions import Fraction


MODEL_DIMS = {
    "A": (1, 3, 2),
    "B": (2, 2, 1),
    "C": (2, 2, 3),
    "D": (3, 1, 2),
    "E": (3, 3, 2),
}


def ftxt(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def leg_factors(dS: int, dF: int):
    dgt = max(dS, dF)
    return {
        "GS": Fraction(dgt, dS),
        "GL": Fraction(dgt, 2),
        "GF": Fraction(dgt, dF),
    }


def predicted(d1: int, d2: int, dF: int):
    f1 = leg_factors(d1, dF)
    f2 = leg_factors(d2, dF)
    return {
        "Tr(y2 y2dag) y2": f2["GS"],
        "y2 y2dag y2": Fraction(1, 2) * (f2["GL"] + f2["GF"]),
        "y1 y1dag y2": Fraction(1, 2) * f1["GL"],
        "Ye Yedag y2": Fraction(1, 2),
    }


def model_alpha(path: Path):
    name = path.parents[1].name
    m = re.fullmatch(r"T3_([A-E])_alpha_([mp])(\d+)", name)
    if not m:
        raise ValueError(name)
    model, sign, mag = m.groups()
    alpha = int(mag)
    return model, -alpha if sign == "m" else alpha


def split_top_level_additive(text: str) -> list[str]:
    """Split TeXForm expression at top-level + or - while preserving signs."""
    parts = []
    start = 0
    par = brk = brc = 0

    for i, ch in enumerate(text):
        if ch == "(":
            par += 1
        elif ch == ")":
            par -= 1
        elif ch == "[":
            brk += 1
        elif ch == "]":
            brk -= 1
        elif ch == "{":
            brc += 1
        elif ch == "}":
            brc -= 1
        elif ch in "+-" and par == brk == brc == 0 and i > start:
            parts.append(text[start:i].strip())
            start = i

    tail = text[start:].strip()
    if tail:
        parts.append(tail)
    return parts


def canonical_term(term: str) -> str:
    # Formatting-only canonicalisation; no physics parsing.
    t = re.sub(r"\s+", "", term)
    t = t.replace(r"\left", "").replace(r"\right", "")
    return t


def is_gauge_term(term: str) -> bool:
    t = canonical_term(term)
    # Covers raw TeXForm variants such as \text{g2}, g2, \text{gY}, gY.
    return (
        "g2" in t
        or "gY" in t
        or r"g_2" in t
        or r"g_Y" in t
        or r"g_{2}" in t
        or r"g_{Y}" in t
    )


def nongauge_from_latex(raw_latex: str) -> tuple[list[str], str]:
    terms = split_top_level_additive(str(raw_latex))
    nongauge_terms = [term for term in terms if not is_gauge_term(term)]
    canonical = "|".join(canonical_term(t) for t in nongauge_terms)
    return nongauge_terms, canonical


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--root", type=Path, default=Path("output/full/hypercharge")
    )
    parser.add_argument(
        "--json-output",
        type=Path,
        default=Path("output/group_factors/y2_nongauge_diagnostic.json"),
    )
    args = parser.parse_args()

    paths = sorted(args.root.glob("T3_*_alpha_*/data/uv_rgbeta_rge.json"))
    if not paths:
        raise SystemExit(f"No UV RGBeta outputs found below {args.root}")

    rows = []
    by_model = {}

    print("beta_y2 nongauge diagnostic")
    print()

    for path in paths:
        payload = json.loads(path.read_text(encoding="utf-8-sig"))
        if payload.get("status") != "Success":
            continue

        model, alpha = model_alpha(path)
        d1, d2, dF = MODEL_DIMS[model]

        latex_map = payload.get("report_beta_latex", {}) or {}
        raw_latex = latex_map.get("y2")
        if not raw_latex:
            raise RuntimeError(f"No report_beta_latex['y2'] in {path}")

        terms, canonical = nongauge_from_latex(raw_latex)
        pred = predicted(d1, d2, dF)

        print(f"T3-{model} alpha={alpha:+d}")
        print("  RGBeta nongauge terms:")
        for term in terms:
            print(f"    {term}")
        print(
            "  predicted coeffs: "
            + ", ".join(f"{name}={ftxt(value)}" for name, value in pred.items())
        )

        rows.append({
            "model": model,
            "alpha": alpha,
            "dS1": d1,
            "dS2": d2,
            "dF": dF,
            "nongauge_terms": terms,
            "canonical_nongauge": canonical,
            "predicted_coefficients": {
                name: ftxt(value) for name, value in pred.items()
            },
            "source": str(path),
        })
        by_model.setdefault(model, set()).add(canonical)

    print()
    print("alpha-independence check")
    overall = True
    for model in sorted(by_model):
        n = len(by_model[model])
        ok = n == 1
        overall &= ok
        print(f"  T3-{model}: {'PASS' if ok else 'FAIL'} ({n} distinct nongauge expressions)")

    result = {
        "status": "Success" if rows and overall else "NeedsInspection",
        "alpha_independent": overall,
        "distinct_nongauge_expression_count": {
            model: len(values) for model, values in sorted(by_model.items())
        },
        "rows": rows,
    }

    args.json_output.parent.mkdir(parents=True, exist_ok=True)
    args.json_output.write_text(json.dumps(result, indent=2), encoding="utf-8")

    print()
    print(f"JSON summary: {args.json_output}")
    print(f"OVERALL:      {result['status']}")
    return 0 if rows else 1


if __name__ == "__main__":
    raise SystemExit(main())
