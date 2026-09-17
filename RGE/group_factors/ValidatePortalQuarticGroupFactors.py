from __future__ import annotations

"""Validate the resolved analytic backbone of beta_lambdaH1/H2/lambda12.

The validator compares the generic formulas against saved RGBeta UV outputs.

Validated here:
  * Yukawa wave-function trace factors
  * mixed Yukawa box factors
  * |lambdaT3|^2 factors
  * gauge-linear g^2 lambda factors
  * singlet scalar-quartic products

Not yet claimed here:
  * pure gauge g^4 structures
  * adjoint-quartic recouplings
  * T3-E Cross mixing structures

Those are the next tensor-recoupling stage.
"""

import argparse
from fractions import Fraction
import json
from pathlib import Path
import re
import sys

import sympy as sp

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from RGE.group_factors.PortalQuarticGroupFactors import portal_quartic_backbone


MODEL_DIMS = {
    "A": (1, 3, 2),
    "B": (2, 2, 1),
    "C": (2, 2, 3),
    "D": (3, 1, 2),
    "E": (3, 3, 2),
}


def F(x: str) -> Fraction:
    return Fraction(str(x))


def txt(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def matching_bracket(text: str, open_index: int) -> int:
    depth = 0
    for i in range(open_index, len(text)):
        if text[i] == "[":
            depth += 1
        elif text[i] == "]":
            depth -= 1
            if depth == 0:
                return i
    raise ValueError("unbalanced brackets")


def replace_traces(expr: str) -> str:
    mapping = {
        "y1.Trans[Bar[y1]]": "Ty1",
        "y2.Trans[Bar[y2]]": "Ty2",
        "yd.Trans[Bar[yd]]": "Tyd",
        "ye.Trans[Bar[ye]]": "Tye",
        "yu.Trans[Bar[yu]]": "Tyu",
        "y1.Trans[Bar[y1]].ye.Trans[Bar[ye]]": "Ty1Ye",
        "y2.Trans[Bar[y2]].ye.Trans[Bar[ye]]": "Ty2Ye",
        "y1.Trans[Bar[y1]].y2.Trans[Bar[y2]]": "Ty1y2",
    }

    result = expr
    replacements = []
    pos = 0
    unknown = 0
    while True:
        start = result.find("Tr[", pos)
        if start < 0:
            break
        open_index = start + 2
        close = matching_bracket(result, open_index)
        inner = result[open_index + 1:close].replace(" ", "")
        token = mapping.get(inner, f"Tother{unknown}")
        unknown += 1
        replacements.append((start, close + 1, token))
        pos = close + 1

    for start, end, token in reversed(replacements):
        result = result[:start] + token + result[end:]
    return result


def parse_beta(raw: str) -> sp.Expr:
    text = replace_traces(str(raw))
    text = text.replace("Bar[lambdaT3]", "lambdaT3Bar")
    text = text.replace("^", "**")

    # Remove any remaining Bar[...] objects without losing algebraic separation.
    count = 0
    while "Bar[" in text:
        start = text.find("Bar[")
        open_index = start + 3
        close = matching_bracket(text, open_index)
        inner = text[open_index + 1:close]
        safe = re.sub(r"[^A-Za-z0-9_]", "_", inner)
        token = f"Bar_{safe}_{count}"
        count += 1
        text = text[:start] + token + text[close + 1:]

    names = set(re.findall(r"\b[A-Za-z][A-Za-z0-9_]*\b", text))
    locals_ = {name: sp.Symbol(name) for name in names}
    return sp.expand(sp.sympify(text, locals=locals_))


def coeff(expr: sp.Expr, monomial: dict[str, int]) -> Fraction:
    out = expr
    target = {sp.Symbol(k): v for k, v in monomial.items()}
    for sym, power in target.items():
        out = sp.expand(out).coeff(sym, power)
    for sym in list(out.free_symbols):
        if sym not in target:
            out = out.subs(sym, 0)
    out = sp.simplify(out)
    return Fraction(int(sp.numer(out)), int(sp.denom(out)))


def model_alpha(path: Path) -> tuple[str, int]:
    name = path.parents[1].name
    m = re.fullmatch(r"T3_([A-E])_alpha_([mp])(\d+)", name)
    if not m:
        raise ValueError(name)
    model, sign, mag = m.groups()
    alpha = int(mag)
    return model, (-alpha if sign == "m" else alpha)


def load_betas(path: Path) -> dict[str, sp.Expr]:
    payload = json.loads(path.read_text(encoding="utf-8-sig"))
    if payload.get("status") != "Success":
        raise ValueError(f"RGBeta failed: {path}")
    report = payload.get("report_betas", {}) or {}
    return {
        name: parse_beta(report[name])
        for name in ("lambdaH1", "lambdaH2", "lambda12")
        if name in report
    }


def expected_checks(p):
    return {
        "lambdaH1": {
            "Ty1*lambdaH1": ({"Ty1": 1, "lambdaH1": 1}, p.H1_y1_trace),
            "Tyd*lambdaH1": ({"Tyd": 1, "lambdaH1": 1}, p.H1_Yd_trace),
            "Tye*lambdaH1": ({"Tye": 1, "lambdaH1": 1}, p.H1_Ye_trace),
            "Tyu*lambdaH1": ({"Tyu": 1, "lambdaH1": 1}, p.H1_Yu_trace),
            "Ty1Ye": ({"Ty1Ye": 1}, p.H1_y1Ye_mixed),
            "|lambdaT3|^2": ({"lambdaT3": 1, "lambdaT3Bar": 1}, p.H1_lambdaT3_sq),
            "g2^2*lambdaH1": ({"g2": 2, "lambdaH1": 1}, p.H1_g2_linear),
            "gY^2*lambdaH1": ({"gY": 2, "lambdaH1": 1}, p.H1_gY_linear),
            "lambdaH*lambdaH1": ({"lambdaH": 1, "lambdaH1": 1}, p.H1_lambdaH_lambdaH1),
            "lambdaH1^2": ({"lambdaH1": 2}, p.H1_lambdaH1_sq),
            "lambdaH1*lambdaS1": ({"lambdaH1": 1, "lambdaS1": 1}, p.H1_lambdaH1_lambdaS1),
            "lambda12*lambdaH2": ({"lambda12": 1, "lambdaH2": 1}, p.H1_lambda12_lambdaH2),
        },
        "lambdaH2": {
            "Ty2*lambdaH2": ({"Ty2": 1, "lambdaH2": 1}, p.H2_y2_trace),
            "Tyd*lambdaH2": ({"Tyd": 1, "lambdaH2": 1}, p.H2_Yd_trace),
            "Tye*lambdaH2": ({"Tye": 1, "lambdaH2": 1}, p.H2_Ye_trace),
            "Tyu*lambdaH2": ({"Tyu": 1, "lambdaH2": 1}, p.H2_Yu_trace),
            "Ty2Ye": ({"Ty2Ye": 1}, p.H2_y2Ye_mixed),
            "|lambdaT3|^2": ({"lambdaT3": 1, "lambdaT3Bar": 1}, p.H2_lambdaT3_sq),
            "g2^2*lambdaH2": ({"g2": 2, "lambdaH2": 1}, p.H2_g2_linear),
            "gY^2*lambdaH2": ({"gY": 2, "lambdaH2": 1}, p.H2_gY_linear),
            "lambdaH*lambdaH2": ({"lambdaH": 1, "lambdaH2": 1}, p.H2_lambdaH_lambdaH2),
            "lambdaH2^2": ({"lambdaH2": 2}, p.H2_lambdaH2_sq),
            "lambdaH2*lambdaS2": ({"lambdaH2": 1, "lambdaS2": 1}, p.H2_lambdaH2_lambdaS2),
            "lambda12*lambdaH1": ({"lambda12": 1, "lambdaH1": 1}, p.H2_lambda12_lambdaH1),
        },
        "lambda12": {
            "Ty1*lambda12": ({"Ty1": 1, "lambda12": 1}, p.L12_y1_trace),
            "Ty2*lambda12": ({"Ty2": 1, "lambda12": 1}, p.L12_y2_trace),
            "Ty1y2": ({"Ty1y2": 1}, p.L12_y1y2_mixed),
            "|lambdaT3|^2": ({"lambdaT3": 1, "lambdaT3Bar": 1}, p.L12_lambdaT3_sq),
            "g2^2*lambda12": ({"g2": 2, "lambda12": 1}, p.L12_g2_linear),
            "gY^2*lambda12": ({"gY": 2, "lambda12": 1}, p.L12_gY_linear),
            "lambda12^2": ({"lambda12": 2}, p.L12_lambda12_sq),
            "lambda12*lambdaS1": ({"lambda12": 1, "lambdaS1": 1}, p.L12_lambda12_lambdaS1),
            "lambda12*lambdaS2": ({"lambda12": 1, "lambdaS2": 1}, p.L12_lambda12_lambdaS2),
            "lambdaH1*lambdaH2": ({"lambdaH1": 1, "lambdaH2": 1}, p.L12_lambdaH1_lambdaH2),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("output/full/hypercharge"))
    parser.add_argument(
        "--json-output",
        type=Path,
        default=Path("output/group_factors/portal_quartic_backbone_validation.json"),
    )
    args = parser.parse_args()

    paths = sorted(args.root.glob("T3_*_alpha_*/data/uv_rgbeta_rge.json"))
    if not paths:
        raise SystemExit(f"No RGBeta outputs found under {args.root}")

    rows = []
    overall = True

    print("Resolved portal-quartic backbone validation")
    print("model alpha  beta       checks  status")
    print("-" * 45)

    for path in paths:
        model, alpha = model_alpha(path)
        dS1, dS2, dF = MODEL_DIMS[model]
        pred = portal_quartic_backbone(dS1, dS2, dF, alpha)
        betas = load_betas(path)
        checks_by_beta = expected_checks(pred)

        for beta_name, checks in checks_by_beta.items():
            if beta_name not in betas:
                continue
            beta = betas[beta_name]
            details = {}
            ok = True
            for label, (monomial, expected_text) in checks.items():
                got = coeff(beta, monomial)
                expected = F(expected_text)
                residual = got - expected
                passed = residual == 0
                ok &= passed
                details[label] = {
                    "expected": txt(expected),
                    "rgbeta": txt(got),
                    "residual": txt(residual),
                    "match": passed,
                }

            overall &= ok
            print(
                f"T3-{model} {alpha:>5}  {beta_name:<9} "
                f"{len(checks):>6}  {'PASS' if ok else 'FAIL'}"
            )
            rows.append({
                "model": model,
                "alpha": alpha,
                "beta": beta_name,
                "source": str(path),
                "checks": details,
                "match": ok,
            })

    payload = {
        "status": "Success" if overall else "Failed",
        "scope": {
            "included": [
                "Yukawa wave-function trace factors",
                "mixed Yukawa box factors",
                "|lambdaT3|^2 factors",
                "g^2 lambda gauge-linear factors",
                "singlet scalar-quartic products",
            ],
            "deferred": [
                "pure gauge g^4 structures",
                "adjoint-quartic recouplings",
                "T3-E Cross mixing structures",
            ],
        },
        "checks": rows,
    }

    args.json_output.parent.mkdir(parents=True, exist_ok=True)
    args.json_output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print()
    print(f"JSON summary: {args.json_output}")
    print(f"OVERALL:      {payload['status']}")
    return 0 if overall else 1


if __name__ == "__main__":
    raise SystemExit(main())
