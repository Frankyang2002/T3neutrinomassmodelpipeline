from __future__ import annotations

"""Validate analytic beta_lambdaT3 group factors.

This validator now owns the compact closed-form lambdaT3 group-factor formula
that previously lived in ``MixingQuarticGroupFactors.py``.

There are still two logically independent validation layers:

1. Physical analytic validation
   The non-singlet recouplings are compared with the independently derived
   gauge-covariant A--E values:
       A: H2Adj = 4
       B/C: H1Adj = -2, H2Adj = 2, 12Adj = -1/2
       D: H1Adj = -4
       E: H1Adj = -2, H2Adj = 2, 12Adj = 2, Cross = -2

2. RGBeta diagnostic comparison
   If saved RGBeta outputs are available, gauge terms and the resolved
   singlet/Yukawa terms are compared exactly.  Non-singlet RGBeta values are
   reported separately because the raw RGBeta formal-index convention differs
   from the physical component-basis result in pseudoreal/real cases.

Independent derivations remain separate:
    RGBetaMixingQuarticRecoupling.py
    MixingQuarticRecoupling.py

A mismatch in the RGBeta non-singlet diagnostic table does not make the
analytic physical validation fail.
"""

import argparse
from dataclasses import dataclass
from fractions import Fraction
import json
from pathlib import Path
import re
from typing import Any

import sympy as sp

from RGE.group_factors.RepresentationFactors import (
    su2_quadratic_casimir_from_dimension,
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


@dataclass(frozen=True)
class MixingQuarticGroupFactors:
    dS1: int
    dS2: int
    dF: int
    alpha: int

    GS1: str
    GS2: str

    C2H: str
    C2S1: str
    C2S2: str

    YH: str
    YS1: str
    YS2: str

    y1_trace: str
    y2_trace: str
    yd_trace: str
    ye_trace: str
    yu_trace: str

    lambdaH: str
    lambdaH1: str
    lambdaH2: str
    lambda12: str

    lambdaH1Adj: str
    lambdaH2Adj: str
    lambda12Adj: str
    lambda12Cross: str

    has_lambdaH1Adj: bool
    has_lambdaH2Adj: bool
    has_lambda12Adj: bool
    has_lambda12Cross: bool

    su2_gauge: str
    u1_gauge: str


def _txt(value: Fraction) -> str:
    return (
        str(value.numerator)
        if value.denominator == 1
        else f"{value.numerator}/{value.denominator}"
    )


def _f(value: str) -> Fraction:
    return Fraction(str(value))


def _validate_t3_dimensions(dS1: int, dS2: int, dF: int) -> None:
    if min(dS1, dS2, dF) < 1:
        raise ValueError("SU(2) representation dimensions must be positive.")

    if abs(dS1 - dF) != 1 or abs(dS2 - dF) != 1:
        raise ValueError("T3 Yukawa invariance requires dSi=dF+/-1.")


def mixing_quartic_group_factors(
    dS1: int,
    dS2: int,
    dF: int,
    alpha: int,
) -> MixingQuarticGroupFactors:
    """Return the closed-form one-loop group factors in beta_lambdaT3."""

    dS1, dS2, dF, alpha = map(int, (dS1, dS2, dF, alpha))
    _validate_t3_dimensions(dS1, dS2, dF)

    GS1 = Fraction(max(dF, dS1), dS1)
    GS2 = Fraction(max(dF, dS2), dS2)

    C2H = Fraction(3, 4)
    C2S1 = su2_quadratic_casimir_from_dimension(dS1)
    C2S2 = su2_quadratic_casimir_from_dimension(dS2)

    YH = Fraction(1, 2)
    YS1 = Fraction(alpha, 2)
    YS2 = Fraction(alpha + 2, 2)

    # Gauge wave-function contribution:
    #   -3 g^2 sum_external C2
    # and analogously for U(1) charges.
    su2 = -3 * (2 * C2H + C2S1 + C2S2)
    u1 = -3 * (2 * YH * YH + YS1 * YS1 + YS2 * YS2)

    has_h1_adj = dS1 > 1
    has_h2_adj = dS2 > 1
    has_12_adj = dS1 > 1 and dS2 > 1
    has_cross = dS1 == 3 and dS2 == 3

    # Gauge-covariant scalar recouplings in the physical component basis.
    r_h1 = C2S2 - C2S1 - 2 if has_h1_adj else Fraction(0)
    r_h2 = 2 - C2S1 + C2S2 if has_h2_adj else Fraction(0)
    r_12 = C2S1 + C2S2 - 2 if has_12_adj else Fraction(0)
    r_cross = Fraction(-2) if has_cross else Fraction(0)

    return MixingQuarticGroupFactors(
        dS1=dS1,
        dS2=dS2,
        dF=dF,
        alpha=alpha,
        GS1=_txt(GS1),
        GS2=_txt(GS2),
        C2H=_txt(C2H),
        C2S1=_txt(C2S1),
        C2S2=_txt(C2S2),
        YH=_txt(YH),
        YS1=_txt(YS1),
        YS2=_txt(YS2),
        y1_trace=_txt(GS1),
        y2_trace=_txt(GS2),
        yd_trace="6",
        ye_trace="2",
        yu_trace="6",
        lambdaH="2",
        lambdaH1="4",
        lambdaH2="4",
        lambda12="2",
        lambdaH1Adj=_txt(r_h1),
        lambdaH2Adj=_txt(r_h2),
        lambda12Adj=_txt(r_12),
        lambda12Cross=_txt(r_cross),
        has_lambdaH1Adj=has_h1_adj,
        has_lambdaH2Adj=has_h2_adj,
        has_lambda12Adj=has_12_adj,
        has_lambda12Cross=has_cross,
        su2_gauge=_txt(su2),
        u1_gauge=_txt(u1),
    )


def _model_alpha(path: Path) -> tuple[str, int]:
    name = path.parents[1].name
    match = re.fullmatch(r"T3_([A-E])_alpha_([mp])(\d+)", name)
    if not match:
        raise ValueError(name)

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


def _beta(payload: dict[str, Any]) -> str:
    raw = (payload.get("report_betas", {}) or {}).get("lambdaT3")
    if not raw:
        raise KeyError("Missing report_betas['lambdaT3']")
    return str(raw)


def _matching_bracket(text: str, open_index: int) -> int:
    depth = 0
    for index in range(open_index, len(text)):
        if text[index] == "[":
            depth += 1
        elif text[index] == "]":
            depth -= 1
            if depth == 0:
                return index
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

    for symbol, power in target_symbols.items():
        result = sp.expand(result).coeff(symbol, power)

    for symbol in list(result.free_symbols):
        if symbol not in target_symbols:
            result = result.subs(symbol, 0)

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
        prediction = mixing_quartic_group_factors(
            dS1,
            dS2,
            dF,
            0,
        )
        expected = PHYSICAL_RECOUPLING_TARGETS[model]

        actual = {
            "lambdaH1Adj": _f(prediction.lambdaH1Adj),
            "lambdaH2Adj": _f(prediction.lambdaH2Adj),
            "lambda12Adj": _f(prediction.lambda12Adj),
            "lambda12Cross": _f(prediction.lambda12Cross),
        }

        expected_full = {
            "lambdaH1Adj": expected.get(
                "lambdaH1Adj",
                Fraction(0),
            ),
            "lambdaH2Adj": expected.get(
                "lambdaH2Adj",
                Fraction(0),
            ),
            "lambda12Adj": expected.get(
                "lambda12Adj",
                Fraction(0),
            ),
            "lambda12Cross": expected.get(
                "lambda12Cross",
                Fraction(0),
            ),
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

        checks.append(
            {
                "kind": "physical_recoupling",
                "model": model,
                "predicted": {
                    key: _txt(value)
                    for key, value in actual.items()
                },
                "independent_target": {
                    key: _txt(value)
                    for key, value in expected_full.items()
                },
                "residuals": {
                    key: _txt(value)
                    for key, value in residuals.items()
                },
                "match": ok,
            }
        )

    return checks, all_ok


def _rgbeta_diagnostics(
    root: Path,
) -> tuple[
    list[dict[str, Any]],
    dict[str, dict[str, str]],
    bool,
]:
    paths = sorted(
        root.glob("T3_*_alpha_*/data/uv_rgbeta_rge.json")
    )

    if not paths:
        print()
        print(f"RGBeta diagnostics skipped: no files under {root}")
        return [], {}, True

    checks: list[dict[str, Any]] = []
    all_resolved_ok = True

    print()
    print("Resolved RGBeta comparison across hypercharge scan")
    print(
        "model alpha  gY(exp/rgb/res)       "
        "g2(exp/rgb/res)       status"
    )
    print("-" * 78)

    for path in paths:
        model, alpha = _model_alpha(path)
        dS1, dS2, dF = MODEL_DIMS[model]

        prediction = mixing_quartic_group_factors(
            dS1,
            dS2,
            dF,
            alpha,
        )
        parsed = _parse(_beta(_load(path)))

        got_y = _coeff(parsed, {"gY": 2, "lambdaT3": 1})
        got_2 = _coeff(parsed, {"g2": 2, "lambdaT3": 1})
        exp_y = _f(prediction.u1_gauge)
        exp_2 = _f(prediction.su2_gauge)

        residual_y = got_y - exp_y
        residual_2 = got_2 - exp_2
        ok = residual_y == 0 and residual_2 == 0
        all_resolved_ok &= ok

        print(
            f"T3-{model} {alpha:>5}  "
            f"{_txt(exp_y):>6}/{_txt(got_y):>6}/"
            f"{_txt(residual_y):>3}       "
            f"{_txt(exp_2):>6}/{_txt(got_2):>6}/"
            f"{_txt(residual_2):>3}       "
            f"{'PASS' if ok else 'FAIL'}"
        )

        checks.append(
            {
                "kind": "rgbeta_gauge",
                "model": model,
                "alpha": alpha,
                "expected_u1": _txt(exp_y),
                "rgbeta_u1": _txt(got_y),
                "u1_residual": _txt(residual_y),
                "expected_su2": _txt(exp_2),
                "rgbeta_su2": _txt(got_2),
                "su2_residual": _txt(residual_2),
                "match": ok,
                "source": str(path),
            }
        )

    print()
    print("Resolved non-gauge RGBeta comparison at alpha=0")
    print("model  y1  y2  yd  ye  yu  lH  lH1 lH2 l12  status")
    print("-" * 66)

    raw_recouplings: dict[str, dict[str, str]] = {}

    for model, (dS1, dS2, dF) in MODEL_DIMS.items():
        path = (
            root
            / f"T3_{model}_alpha_p0"
            / "data"
            / "uv_rgbeta_rge.json"
        )
        if not path.exists():
            continue

        parsed = _parse(_beta(_load(path)))
        prediction = mixing_quartic_group_factors(
            dS1,
            dS2,
            dF,
            0,
        )

        names = (
            ("Ty1", prediction.y1_trace),
            ("Ty2", prediction.y2_trace),
            ("Tyd", prediction.yd_trace),
            ("Tye", prediction.ye_trace),
            ("Tyu", prediction.yu_trace),
            ("lambdaH", prediction.lambdaH),
            ("lambdaH1", prediction.lambdaH1),
            ("lambdaH2", prediction.lambdaH2),
            ("lambda12", prediction.lambda12),
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
            + " ".join(f"{value:>3}" for value in got_values)
            + f"  {'PASS' if row_ok else 'FAIL'}"
        )

        checks.append(
            {
                "kind": "rgbeta_resolved_non_gauge",
                "model": model,
                "alpha": 0,
                "residuals": residuals,
                "match": row_ok,
                "source": str(path),
            }
        )

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

    physical_checks, physical_ok = (
        _physical_recoupling_validation()
    )
    (
        rgbeta_checks,
        rgbeta_raw,
        rgbeta_resolved_ok,
    ) = _rgbeta_diagnostics(args.root)

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
        "independent_derivations": [
            "RGE/group_factors/RGBetaMixingQuarticRecoupling.py",
            "RGE/group_factors/MixingQuarticRecoupling.py",
        ],
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
