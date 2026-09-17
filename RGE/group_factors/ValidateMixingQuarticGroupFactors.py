from __future__ import annotations

"""Validate analytic beta_lambdaT3 group factors.

This validator has two logically separate layers.

1. Physical analytic validation
   The newly implemented non-singlet recouplings are checked against the
   independently derived gauge-covariant A--E values:
       A: H2Adj = 4
       B/C: H1Adj = -2, H2Adj = 2, 12Adj = -1/2
       D: H1Adj = -4
       E: H1Adj = -2, H2Adj = 2, 12Adj = 2, Cross = -2

   These are derived checks, not the source of the implementation.

2. RGBeta diagnostic comparison
   If saved RGBeta outputs are available, gauge terms and the already-resolved
   singlet/Yukawa terms are compared exactly.  Non-singlet RGBeta values are
   reported separately because the raw RGBeta formal-index convention differs
   from the physical component-basis result in the pseudoreal/real cases.

A mismatch in the RGBeta non-singlet diagnostic table therefore does not make
the analytic physical validation fail.
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

from RGE.group_factors.MixingQuarticGroupFactors import (
    mixing_quartic_group_factors,
)


MODEL_DIMS = {
    "A": (1, 3, 2),
    "B": (2, 2, 1),
    "C": (2, 2, 3),
    "D": (3, 1, 2),
    "E": (3, 3, 2),
}

PHYSICAL_RECOUPLING_TARGETS = {
    "A": {
        "lambdaH2Adj": Fraction(4),
    },
    "B": {
        "lambdaH1Adj": Fraction(-2),
        "lambdaH2Adj": Fraction(2),
        "lambda12Adj": Fraction(-1, 2),
    },
    "C": {
        "lambdaH1Adj": Fraction(-2),
        "lambdaH2Adj": Fraction(2),
        "lambda12Adj": Fraction(-1, 2),
    },
    "D": {
        "lambdaH1Adj": Fraction(-4),
    },
    "E": {
        "lambdaH1Adj": Fraction(-2),
        "lambdaH2Adj": Fraction(2),
        "lambda12Adj": Fraction(2),
        "lambda12Cross": Fraction(-2),
    },
}


def _txt(value: Fraction) -> str:
    return (
        str(value.numerator)
        if value.denominator == 1
        else f"{value.numerator}/{value.denominator}"
    )


def _f(value: str) -> Fraction:
    return Fraction(str(value))


def _model_alpha(path: Path) -> tuple[str, int]:
    name = path.parents[1].name
    match = re.fullmatch(r"T3_([A-E])_alpha_([mp])(\d+)", name)
    if not match:
        raise ValueError(name)
    model, sign, mag = match.groups()
    alpha = int(mag)
    if sign == "m":
        alpha = -alpha
    return model, alpha


def _load(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8-sig"))
    if payload.get("status") != "Success":
        raise ValueError(f"RGBeta status is not Success: {path}")
    return payload


def _beta(payload: dict[str, Any]) -> str:
    raw = (payload.get("report_betas", {}) or {}).get("lambdaT3")
    if not raw:
        raise KeyError("Missing report_betas['lambdaT3']")
    return str(raw)


def _matching_bracket(text: str, open_index: int) -> int:
    depth = 0
    for i in range(open_index, len(text)):
        if text[i] == "[":
            depth += 1
        elif text[i] == "]":
            depth -= 1
            if depth == 0:
                return i
    raise ValueError("Unbalanced []")


def _replace_traces(expr: str) -> str:
    result = expr
    replacements: list[tuple[int, int, str]] = []
    pos = 0
    count = 0

    while True:
        start = result.find("Tr[", pos)
        if start < 0:
            break
        open_index = start + 2
        close = _matching_bracket(result, open_index)
        inner = result[open_index + 1:close].replace(" ", "")

        mapping = {
            "y1.Trans[Bar[y1]]": "Ty1",
            "y2.Trans[Bar[y2]]": "Ty2",
            "yd.Trans[Bar[yd]]": "Tyd",
            "ye.Trans[Bar[ye]]": "Tye",
            "yu.Trans[Bar[yu]]": "Tyu",
        }
        token = mapping.get(inner, f"Tother{count}")
        count += 1
        replacements.append((start, close + 1, token))
        pos = close + 1

    for start, end, token in reversed(replacements):
        result = result[:start] + token + result[end:]

    return result


def _replace_bar_calls(expr: str) -> str:
    result = expr
    replacements: list[tuple[int, int, str]] = []
    pos = 0
    count = 0

    while True:
        start = result.find("Bar[", pos)
        if start < 0:
            break
        open_index = start + len("Bar")
        close = _matching_bracket(result, open_index)
        inner = result[open_index + 1:close].strip()
        safe = re.sub(r"[^A-Za-z0-9_]", "_", inner)
        token = f"Bar_{safe}_{count}"
        count += 1
        replacements.append((start, close + 1, token))
        pos = close + 1

    for start, end, token in reversed(replacements):
        result = result[:start] + token + result[end:]

    return result


def _parse(expr: str) -> sp.Expr:
    cleaned = _replace_traces(expr)
    cleaned = _replace_bar_calls(cleaned)
    cleaned = cleaned.replace("^", "**")

    names = set(re.findall(r"\b[A-Za-z][A-Za-z0-9_]*\b", cleaned))
    locals_ = {name: sp.Symbol(name) for name in names}
    return sp.expand(sp.sympify(cleaned, locals=locals_))


def _coeff(expr: sp.Expr, monomial: dict[str, int]) -> Fraction:
    result = expr

    target_symbols = {
        sp.Symbol(name): power
        for name, power in monomial.items()
    }

    for sym, power in target_symbols.items():
        result = sp.expand(result).coeff(sym, power)

    for sym in list(result.free_symbols):
        if sym not in target_symbols:
            result = result.subs(sym, 0)

    result = sp.simplify(result)
    return Fraction(int(sp.numer(result)), int(sp.denom(result)))


def _single_coeff(parsed: sp.Expr, factor: str) -> Fraction:
    return _coeff(parsed, {"lambdaT3": 1, factor: 1})


def _physical_recoupling_validation() -> tuple[list[dict[str, Any]], bool]:
    checks: list[dict[str, Any]] = []
    all_ok = True

    print("Physical analytic non-singlet recouplings")
    print("model  H1Adj   H2Adj   12Adj   Cross   status")
    print("-" * 54)

    for model, (dS1, dS2, dF) in MODEL_DIMS.items():
        pred = mixing_quartic_group_factors(dS1, dS2, dF, 0)
        expected = PHYSICAL_RECOUPLING_TARGETS[model]

        actual = {
            "lambdaH1Adj": _f(pred.lambdaH1Adj),
            "lambdaH2Adj": _f(pred.lambdaH2Adj),
            "lambda12Adj": _f(pred.lambda12Adj),
            "lambda12Cross": _f(pred.lambda12Cross),
        }

        expected_full = {
            "lambdaH1Adj": expected.get("lambdaH1Adj", Fraction(0)),
            "lambdaH2Adj": expected.get("lambdaH2Adj", Fraction(0)),
            "lambda12Adj": expected.get("lambda12Adj", Fraction(0)),
            "lambda12Cross": expected.get("lambda12Cross", Fraction(0)),
        }

        residuals = {
            key: actual[key] - expected_full[key]
            for key in actual
        }
        ok = all(value == 0 for value in residuals.values())
        all_ok &= ok

        print(
            f"T3-{model}  "
            f"{_txt(actual['lambdaH1Adj']):>5}   "
            f"{_txt(actual['lambdaH2Adj']):>5}   "
            f"{_txt(actual['lambda12Adj']):>5}   "
            f"{_txt(actual['lambda12Cross']):>5}   "
            f"{'PASS' if ok else 'FAIL'}"
        )

        checks.append({
            "kind": "physical_recoupling",
            "model": model,
            "predicted": {k: _txt(v) for k, v in actual.items()},
            "independent_target": {
                k: _txt(v) for k, v in expected_full.items()
            },
            "residuals": {k: _txt(v) for k, v in residuals.items()},
            "match": ok,
        })

    return checks, all_ok


def _rgbeta_diagnostics(
    root: Path,
) -> tuple[list[dict[str, Any]], dict[str, dict[str, str]], bool]:
    paths = sorted(root.glob("T3_*_alpha_*/data/uv_rgbeta_rge.json"))
    if not paths:
        print()
        print(f"RGBeta diagnostics skipped: no files under {root}")
        return [], {}, True

    checks: list[dict[str, Any]] = []
    all_resolved_ok = True

    print()
    print("Resolved RGBeta comparison across hypercharge scan")
    print("model alpha  gY(exp/rgb/res)       g2(exp/rgb/res)       status")
    print("-" * 78)

    for path in paths:
        model, alpha = _model_alpha(path)
        dS1, dS2, dF = MODEL_DIMS[model]
        pred = mixing_quartic_group_factors(dS1, dS2, dF, alpha)
        parsed = _parse(_beta(_load(path)))

        got_y = _coeff(parsed, {"gY": 2, "lambdaT3": 1})
        got_2 = _coeff(parsed, {"g2": 2, "lambdaT3": 1})
        exp_y = _f(pred.u1_gauge)
        exp_2 = _f(pred.su2_gauge)
        dy = got_y - exp_y
        d2 = got_2 - exp_2
        ok = dy == 0 and d2 == 0
        all_resolved_ok &= ok

        print(
            f"T3-{model} {alpha:>5}  "
            f"{_txt(exp_y):>6}/{_txt(got_y):>6}/{_txt(dy):>3}       "
            f"{_txt(exp_2):>6}/{_txt(got_2):>6}/{_txt(d2):>3}       "
            f"{'PASS' if ok else 'FAIL'}"
        )

        checks.append({
            "kind": "rgbeta_gauge",
            "model": model,
            "alpha": alpha,
            "expected_u1": _txt(exp_y),
            "rgbeta_u1": _txt(got_y),
            "u1_residual": _txt(dy),
            "expected_su2": _txt(exp_2),
            "rgbeta_su2": _txt(got_2),
            "su2_residual": _txt(d2),
            "match": ok,
            "source": str(path),
        })

    print()
    print("Resolved non-gauge RGBeta comparison at alpha=0")
    print("model  y1  y2  yd  ye  yu  lH  lH1 lH2 l12  status")
    print("-" * 66)

    raw_recouplings: dict[str, dict[str, str]] = {}

    for model, (dS1, dS2, dF) in MODEL_DIMS.items():
        path = root / f"T3_{model}_alpha_p0" / "data" / "uv_rgbeta_rge.json"
        if not path.exists():
            continue

        parsed = _parse(_beta(_load(path)))
        pred = mixing_quartic_group_factors(dS1, dS2, dF, 0)

        names = (
            ("Ty1", pred.y1_trace),
            ("Ty2", pred.y2_trace),
            ("Tyd", pred.yd_trace),
            ("Tye", pred.ye_trace),
            ("Tyu", pred.yu_trace),
            ("lambdaH", pred.lambdaH),
            ("lambdaH1", pred.lambdaH1),
            ("lambdaH2", pred.lambdaH2),
            ("lambda12", pred.lambda12),
        )

        row_ok = True
        got_values = []
        residuals = {}
        for factor, expected_text in names:
            got = _single_coeff(parsed, factor)
            expected = _f(expected_text)
            residual = got - expected
            row_ok &= residual == 0
            got_values.append(_txt(got))
            residuals[factor] = _txt(residual)

        all_resolved_ok &= row_ok

        print(
            f"T3-{model}  "
            + " ".join(f"{v:>3}" for v in got_values)
            + f"  {'PASS' if row_ok else 'FAIL'}"
        )

        checks.append({
            "kind": "rgbeta_resolved_non_gauge",
            "model": model,
            "alpha": 0,
            "residuals": residuals,
            "match": row_ok,
            "source": str(path),
        })

        raw = {}
        for factor in (
            "lambdaH1Adj",
            "lambdaH2Adj",
            "lambda12Adj",
            "lambda12Cross",
        ):
            got = _single_coeff(parsed, factor)
            if got != 0:
                raw[factor] = _txt(got)
        raw_recouplings[model] = raw

    return checks, raw_recouplings, all_resolved_ok


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--root",
        type=Path,
        default=Path("output/full/hypercharge"),
    )
    parser.add_argument(
        "--json-output",
        type=Path,
        default=Path(
            "output/group_factors/lambdaT3_analytic_validation.json"
        ),
    )
    args = parser.parse_args()

    physical_checks, physical_ok = _physical_recoupling_validation()
    rgbeta_checks, rgbeta_raw, rgbeta_resolved_ok = _rgbeta_diagnostics(args.root)

    print()
    print("Raw RGBeta non-singlet coefficients (diagnostic only)")
    if rgbeta_raw:
        for model in MODEL_DIMS:
            print(f"T3-{model}: {rgbeta_raw.get(model, {})}")
    else:
        print("(not available)")

    overall_ok = physical_ok and rgbeta_resolved_ok

    payload = {
        "status": "Success" if overall_ok else "Failed",
        "analytic_formula": {
            "R_H1Adj": "C2(S2)-C2(S1)-2",
            "R_H2Adj": "2-C2(S1)+C2(S2)",
            "R_12Adj": "C2(S1)+C2(S2)-2",
            "R_12Cross": (
                "-2 for dS1=dS2=3 in the current T3 operator basis; "
                "0 otherwise"
            ),
            "y1_trace": "GS1",
            "y2_trace": "GS2",
            "SM_Higgs_Yukawa": (
                "6 Tr(YdYd†)+2 Tr(YeYe†)+6 Tr(YuYu†)"
            ),
            "singlet_quartics": (
                "2 lambdaH + 4 lambdaH1 + 4 lambdaH2 + 2 lambda12"
            ),
            "SU2": "-3*(2*C2(H)+C2(S1)+C2(S2))",
            "U1": "-3*(2*YH^2+YS1^2+YS2^2)",
        },
        "physical_recoupling_validation": physical_checks,
        "rgbeta_resolved_checks": rgbeta_checks,
        "raw_rgbeta_non_singlet_diagnostic_only": rgbeta_raw,
        "notes": [
            (
                "The physical non-singlet recoupling formulas are derived in "
                "the gauge-covariant component basis."
            ),
            (
                "Raw RGBeta non-singlet coefficients are retained only as a "
                "diagnostic because its formal-index convention differs in "
                "the pseudoreal/real representation cases."
            ),
        ],
    }

    args.json_output.parent.mkdir(parents=True, exist_ok=True)
    args.json_output.write_text(
        json.dumps(payload, indent=2),
        encoding="utf-8",
    )

    print()
    print(f"JSON summary: {args.json_output}")
    print(f"OVERALL:      {payload['status']}")
    return 0 if overall_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
