from __future__ import annotations

"""Validate the generic T3 heavy-fermion mass beta function against RGBeta.

This version parses RGBeta's saved ``report_betas["MF"]`` structurally by
replacing each complete ``Matrix[...][heavy[i],heavy[j]]`` object with a
symbolic token and letting SymPy extract exact coefficients.  This avoids the
sign/parenthesis failures of the earlier regex-only validator and also handles
the neutral-F branch, where RGBeta writes both MF and Trans[MF].
"""

import argparse
from fractions import Fraction
import json
from pathlib import Path
import re
import sys
from typing import Any

import sympy as sp

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from RGE.group_factors.FermionMassGroupFactors import fermion_mass_group_factors


MODEL_DIMS = {
    "A": (1, 3, 2),
    "B": (2, 2, 1),
    "C": (2, 2, 3),
    "D": (3, 1, 2),
    "E": (3, 3, 2),
}


def _f(text: str) -> Fraction:
    return Fraction(str(text))


def _txt(value: Fraction | None) -> str | None:
    if value is None:
        return None
    return (
        str(value.numerator)
        if value.denominator == 1
        else f"{value.numerator}/{value.denominator}"
    )


def _model_alpha(path: Path) -> tuple[str, int]:
    name = path.parents[1].name
    match = re.fullmatch(r"T3_([A-E])_alpha_([mp])(\d+)", name)
    if not match:
        raise ValueError(f"Unrecognised run directory: {name}")
    model, sign, magnitude = match.groups()
    alpha = int(magnitude)
    if sign == "m":
        alpha = -alpha
    return model, alpha


def _load(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8-sig"))
    if payload.get("status") != "Success":
        raise ValueError(f"RGBeta status is not Success: {path}")
    return payload


def _mf_report(payload: dict[str, Any]) -> str:
    raw = (payload.get("report_betas", {}) or {}).get("MF")
    if not raw:
        raise KeyError("Missing report_betas['MF']")
    return str(raw)


def _matching_bracket(text: str, open_index: int) -> int:
    if text[open_index] != "[":
        raise ValueError("Expected '['")
    depth = 0
    for i in range(open_index, len(text)):
        if text[i] == "[":
            depth += 1
        elif text[i] == "]":
            depth -= 1
            if depth == 0:
                return i
    raise ValueError("Unbalanced square brackets")


def _matrix_objects(expr: str) -> list[tuple[int, int, str]]:
    """Return (start,end,inner) for complete Matrix[inner][indices] objects."""
    found: list[tuple[int, int, str]] = []
    pos = 0
    while True:
        start = expr.find("Matrix[", pos)
        if start < 0:
            break

        first_open = start + len("Matrix")
        first_close = _matching_bracket(expr, first_open)
        inner = expr[first_open + 1:first_close]

        end = first_close + 1
        if end < len(expr) and expr[end] == "[":
            second_close = _matching_bracket(expr, end)
            end = second_close + 1

        found.append((start, end, inner))
        pos = end
    return found


def _replace_matrix_objects(
    expr: str,
    classifier,
) -> tuple[str, dict[str, str]]:
    """Replace every complete Matrix object by one simple SymPy symbol."""
    objects = _matrix_objects(expr)
    replacements: list[tuple[int, int, str]] = []
    labels: dict[str, str] = {}

    for n, (start, end, inner) in enumerate(objects):
        token = classifier(inner, n)
        if token is None:
            token = f"Z{n}"
        replacements.append((start, end, token))
        labels[token] = inner

    result = expr
    for start, end, token in reversed(replacements):
        result = result[:start] + token + result[end:]

    return result, labels


def _sympy_scalar_expr(text: str, extra_symbols: tuple[str, ...]) -> sp.Expr:
    cleaned = text.replace("^", "**")
    locals_ = {name: sp.Symbol(name) for name in extra_symbols}
    locals_.update({
        "g2": sp.Symbol("g2"),
        "gY": sp.Symbol("gY"),
    })
    # Any Z<n> dummy matrix tokens are independent symbols.
    for name in set(re.findall(r"\bZ\d+\b", cleaned)):
        locals_[name] = sp.Symbol(name)

    return sp.expand(sp.sympify(cleaned, locals=locals_))


def _gauge_coeff(expr: str, gauge: str) -> Fraction:
    """Coefficient of g^2 times the physical MF direction.

    The neutral-F branch contains Matrix[MF] + Matrix[Trans[MF]].
    Both are mapped onto the same physical mass token M, i.e. MF^T=MF.
    """
    def classify(inner: str, _n: int) -> str | None:
        compact = inner.replace(" ", "")
        if compact in {"MF", "Trans[MF]"}:
            return "M"
        return None

    reduced, _ = _replace_matrix_objects(expr, classify)
    parsed = _sympy_scalar_expr(reduced, ("M",))
    g = sp.Symbol(gauge)
    M = sp.Symbol("M")

    coeff = sp.expand(parsed).coeff(M, 1).coeff(g, 2)

    # Eliminate all non-gauge dummy symbols; the desired term must be scalar.
    for symbol in coeff.free_symbols:
        if symbol not in {g}:
            coeff = coeff.subs(symbol, 0)
    coeff = sp.simplify(coeff)

    return Fraction(int(sp.numer(coeff)), int(sp.denom(coeff)))


def _yukawa_coeff(expr: str, target_inner: str) -> Fraction:
    """Coefficient of one exact matrix structure in beta_MF."""
    def classify(inner: str, _n: int) -> str | None:
        if inner.replace(" ", "") == target_inner:
            return "X"
        return None

    reduced, _ = _replace_matrix_objects(expr, classify)
    parsed = _sympy_scalar_expr(reduced, ("X",))
    X = sp.Symbol("X")

    coeff = sp.expand(parsed).coeff(X, 1)

    # Set every unrelated symbol/coupling to zero.
    for symbol in list(coeff.free_symbols):
        if symbol != X:
            coeff = coeff.subs(symbol, 0)
    coeff = sp.simplify(coeff)

    return Fraction(int(sp.numer(coeff)), int(sp.denom(coeff)))


def _vectorlike_y1_coeff(expr: str) -> Fraction:
    return _yukawa_coeff(
        expr,
        "Trans[y1],Bar[y1],MF",
    )


def _vectorlike_y2_coeff(expr: str) -> Fraction:
    return _yukawa_coeff(
        expr,
        "MF,Trans[Bar[y2]],y2",
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Validate T3 beta_MF group factors against saved RGBeta results."
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path("output/full/hypercharge"),
    )
    parser.add_argument(
        "--json-output",
        type=Path,
        default=Path("output/group_factors/fermion_mass_validation.json"),
    )
    args = parser.parse_args()

    paths = sorted(args.root.glob("T3_*_alpha_*/data/uv_rgbeta_rge.json"))
    if not paths:
        raise FileNotFoundError(f"No RGBeta files found below {args.root}")

    checks: list[dict[str, Any]] = []
    all_ok = True

    print("Gauge validation across hypercharge scan")
    print("model alpha  gY(exp/rgb/res)   g2(exp/rgb/res)  status")
    print("-" * 70)

    for path in paths:
        model, alpha = _model_alpha(path)
        dS1, dS2, dF = MODEL_DIMS[model]
        pred = fermion_mass_group_factors(dS1, dS2, dF, alpha)
        expr = _mf_report(_load(path))

        got_u1 = _gauge_coeff(expr, "gY")
        got_su2 = _gauge_coeff(expr, "g2")
        exp_u1 = _f(pred.u1_gauge)
        exp_su2 = _f(pred.su2_gauge)

        du1 = got_u1 - exp_u1
        dsu2 = got_su2 - exp_su2
        ok = du1 == 0 and dsu2 == 0
        all_ok &= ok

        print(
            f"T3-{model:1} {alpha:>5}  "
            f"{_txt(exp_u1):>5}/{_txt(got_u1):>5}/{_txt(du1):>3}   "
            f"{_txt(exp_su2):>5}/{_txt(got_su2):>5}/{_txt(dsu2):>3}  "
            f"{'PASS' if ok else 'FAIL'}"
        )

        checks.append({
            "kind": "gauge",
            "model": model,
            "alpha": alpha,
            "expected_u1": _txt(exp_u1),
            "rgbeta_u1": _txt(got_u1),
            "u1_residual": _txt(du1),
            "expected_su2": _txt(exp_su2),
            "rgbeta_su2": _txt(got_su2),
            "su2_residual": _txt(dsu2),
            "match": ok,
            "source": str(path),
        })

    print()
    print("Vector-like Yukawa-leg validation at alpha=0")
    print("model  y1L(exp/rgb/res)  y2R(exp/rgb/res)  status")
    print("-" * 62)

    for model, (dS1, dS2, dF) in MODEL_DIMS.items():
        path = args.root / f"T3_{model}_alpha_p0" / "data" / "uv_rgbeta_rge.json"
        expr = _mf_report(_load(path))
        pred = fermion_mass_group_factors(dS1, dS2, dF, 0)

        got1 = _vectorlike_y1_coeff(expr)
        got2 = _vectorlike_y2_coeff(expr)
        exp1 = _f(pred.y1_left)
        exp2 = _f(pred.y2_right)

        d1 = got1 - exp1
        d2 = got2 - exp2
        ok = d1 == 0 and d2 == 0
        all_ok &= ok

        print(
            f"T3-{model:1}    "
            f"{_txt(exp1):>4}/{_txt(got1):>4}/{_txt(d1):>3}       "
            f"{_txt(exp2):>4}/{_txt(got2):>4}/{_txt(d2):>3}      "
            f"{'PASS' if ok else 'FAIL'}"
        )

        checks.append({
            "kind": "vectorlike_yukawa",
            "model": model,
            "alpha": 0,
            "expected_y1_left": _txt(exp1),
            "rgbeta_y1_left": _txt(got1),
            "y1_residual": _txt(d1),
            "expected_y2_right": _txt(exp2),
            "rgbeta_y2_right": _txt(got2),
            "y2_residual": _txt(d2),
            "match": ok,
            "source": str(path),
        })

    payload = {
        "status": "Success" if all_ok else "Failed",
        "formula_vectorlike_branch": (
            "16*pi^2 beta_MF = 1/2 GF1 (y1^T y1*) MF "
            "+ 1/2 GF2 MF (y2^dagger y2) "
            "- 6 g2^2 C2(F) MF - 6 gY^2 YF^2 MF"
        ),
        "neutral_branch_note": (
            "For YF=0 RGBeta defines the mass as {F,F}; its saved beta contains "
            "MF and Trans[MF]. Gauge validation projects onto MF^T=MF. "
            "The neutral-branch Yukawa tensor structure is not identified with "
            "the vector-like left/right formula in this validator."
        ),
        "checks": checks,
    }

    args.json_output.parent.mkdir(parents=True, exist_ok=True)
    args.json_output.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print()
    print(f"RGBeta files checked: {len(paths)}")
    print(f"JSON summary:         {args.json_output}")
    print(f"OVERALL:              {payload['status']}")
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
