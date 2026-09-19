from __future__ import annotations

"""Validate the consolidated portal-quartic group-factor implementation.

This file combines the former resolved-backbone RGBeta validator with the
consolidated representation-level checks.

Validated against saved RGBeta UV outputs:
  * Yukawa wave-function trace factors
  * mixed Yukawa box factors
  * |lambdaT3|^2 factors
  * gauge-linear g^2 lambda factors
  * singlet scalar-quartic products

Validated internally at representation level:
  * adjoint-squared recouplings
  * T3-E Cross/self-adjoint coefficients
  * pure-gauge representation identities over an alpha scan

The independent derivation scripts remain separate evidence:
  * DerivePortalAdjointRecouplings.py
  * DeriveEPortalCrossScalarRecouplings.py
  * DerivePortalGaugeQuartics.py

No A--E model label is used by the production formula itself; labels appear
here only to select the five supported benchmark representation assignments.
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

from RGE.group_factors.PortalQuarticBetaGroupFactors import (
    portal_quartic_beta_group_factors,
)


MODELS = {
    "A": (1, 3, 2),
    "B": (2, 2, 1),
    "C": (2, 2, 3),
    "D": (3, 1, 2),
    "E": (3, 3, 2),
}


def F(value: Any) -> Fraction:
    return Fraction(str(value))


def txt(value: Fraction) -> str:
    return (
        str(value.numerator)
        if value.denominator == 1
        else f"{value.numerator}/{value.denominator}"
    )


def matching_bracket(text: str, open_index: int) -> int:
    depth = 0
    for index in range(open_index, len(text)):
        if text[index] == "[":
            depth += 1
        elif text[index] == "]":
            depth -= 1
            if depth == 0:
                return index
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
    target = {
        sp.Symbol(name): power
        for name, power in monomial.items()
    }

    for symbol, power in target.items():
        out = sp.expand(out).coeff(symbol, power)

    for symbol in list(out.free_symbols):
        if symbol not in target:
            out = out.subs(symbol, 0)

    out = sp.simplify(out)
    return Fraction(int(sp.numer(out)), int(sp.denom(out)))


def model_alpha(path: Path) -> tuple[str, int]:
    name = path.parents[1].name
    match = re.fullmatch(r"T3_([A-E])_alpha_([mp])(\d+)", name)
    if not match:
        raise ValueError(name)

    model, sign, magnitude = match.groups()
    alpha = int(magnitude)
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


def expected_backbone_checks(prediction):
    h1 = prediction.beta_lambdaH1
    h2 = prediction.beta_lambdaH2
    l12 = prediction.beta_lambda12

    return {
        "lambdaH1": {
            "Ty1*lambdaH1": (
                {"Ty1": 1, "lambdaH1": 1},
                h1["lambdaH1*Tr_y1"],
            ),
            "Tyd*lambdaH1": (
                {"Tyd": 1, "lambdaH1": 1},
                h1["lambdaH1*Tr_Yd"],
            ),
            "Tye*lambdaH1": (
                {"Tye": 1, "lambdaH1": 1},
                h1["lambdaH1*Tr_Ye"],
            ),
            "Tyu*lambdaH1": (
                {"Tyu": 1, "lambdaH1": 1},
                h1["lambdaH1*Tr_Yu"],
            ),
            "Ty1Ye": (
                {"Ty1Ye": 1},
                h1["Tr_y1Ye"],
            ),
            "|lambdaT3|^2": (
                {"lambdaT3": 1, "lambdaT3Bar": 1},
                h1["abs_lambdaT3_sq"],
            ),
            "g2^2*lambdaH1": (
                {"g2": 2, "lambdaH1": 1},
                h1["g2_sq*lambdaH1"],
            ),
            "gY^2*lambdaH1": (
                {"gY": 2, "lambdaH1": 1},
                h1["gY_sq*lambdaH1"],
            ),
            "lambdaH*lambdaH1": (
                {"lambdaH": 1, "lambdaH1": 1},
                h1["lambdaH*lambdaH1"],
            ),
            "lambdaH1^2": (
                {"lambdaH1": 2},
                h1["lambdaH1_sq"],
            ),
            "lambdaH1*lambdaS1": (
                {"lambdaH1": 1, "lambdaS1": 1},
                h1["lambdaH1*lambdaS1"],
            ),
            "lambda12*lambdaH2": (
                {"lambda12": 1, "lambdaH2": 1},
                h1["lambda12*lambdaH2"],
            ),
        },
        "lambdaH2": {
            "Ty2*lambdaH2": (
                {"Ty2": 1, "lambdaH2": 1},
                h2["lambdaH2*Tr_y2"],
            ),
            "Tyd*lambdaH2": (
                {"Tyd": 1, "lambdaH2": 1},
                h2["lambdaH2*Tr_Yd"],
            ),
            "Tye*lambdaH2": (
                {"Tye": 1, "lambdaH2": 1},
                h2["lambdaH2*Tr_Ye"],
            ),
            "Tyu*lambdaH2": (
                {"Tyu": 1, "lambdaH2": 1},
                h2["lambdaH2*Tr_Yu"],
            ),
            "Ty2Ye": (
                {"Ty2Ye": 1},
                h2["Tr_y2Ye"],
            ),
            "|lambdaT3|^2": (
                {"lambdaT3": 1, "lambdaT3Bar": 1},
                h2["abs_lambdaT3_sq"],
            ),
            "g2^2*lambdaH2": (
                {"g2": 2, "lambdaH2": 1},
                h2["g2_sq*lambdaH2"],
            ),
            "gY^2*lambdaH2": (
                {"gY": 2, "lambdaH2": 1},
                h2["gY_sq*lambdaH2"],
            ),
            "lambdaH*lambdaH2": (
                {"lambdaH": 1, "lambdaH2": 1},
                h2["lambdaH*lambdaH2"],
            ),
            "lambdaH2^2": (
                {"lambdaH2": 2},
                h2["lambdaH2_sq"],
            ),
            "lambdaH2*lambdaS2": (
                {"lambdaH2": 1, "lambdaS2": 1},
                h2["lambdaH2*lambdaS2"],
            ),
            "lambda12*lambdaH1": (
                {"lambda12": 1, "lambdaH1": 1},
                h2["lambda12*lambdaH1"],
            ),
        },
        "lambda12": {
            "Ty1*lambda12": (
                {"Ty1": 1, "lambda12": 1},
                l12["lambda12*Tr_y1"],
            ),
            "Ty2*lambda12": (
                {"Ty2": 1, "lambda12": 1},
                l12["lambda12*Tr_y2"],
            ),
            "Ty1y2": (
                {"Ty1y2": 1},
                l12["Tr_y1y2"],
            ),
            "|lambdaT3|^2": (
                {"lambdaT3": 1, "lambdaT3Bar": 1},
                l12["abs_lambdaT3_sq"],
            ),
            "g2^2*lambda12": (
                {"g2": 2, "lambda12": 1},
                l12["g2_sq*lambda12"],
            ),
            "gY^2*lambda12": (
                {"gY": 2, "lambda12": 1},
                l12["gY_sq*lambda12"],
            ),
            "lambda12^2": (
                {"lambda12": 2},
                l12["lambda12_sq"],
            ),
            "lambda12*lambdaS1": (
                {"lambda12": 1, "lambdaS1": 1},
                l12["lambda12*lambdaS1"],
            ),
            "lambda12*lambdaS2": (
                {"lambda12": 1, "lambdaS2": 1},
                l12["lambda12*lambdaS2"],
            ),
            "lambdaH1*lambdaH2": (
                {"lambdaH1": 1, "lambdaH2": 1},
                l12["lambdaH1*lambdaH2"],
            ),
        },
    }


def check_rgbeta_backbone(root: Path):
    paths = sorted(root.glob("T3_*_alpha_*/data/uv_rgbeta_rge.json"))
    if not paths:
        return [], None

    rows = []
    overall = True

    print("Resolved portal-quartic backbone vs RGBeta")
    print("model alpha  beta       checks  status")
    print("-" * 45)

    for path in paths:
        model, alpha = model_alpha(path)
        d_s1, d_s2, d_f = MODELS[model]

        prediction = portal_quartic_beta_group_factors(
            d_s1,
            d_s2,
            d_f,
            alpha,
        )
        betas = load_betas(path)
        checks_by_beta = expected_backbone_checks(prediction)

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

            rows.append(
                {
                    "model": model,
                    "alpha": alpha,
                    "beta": beta_name,
                    "source": str(path),
                    "checks": details,
                    "match": ok,
                }
            )

    return rows, overall


def check_representation_identities():
    rows = []
    overall = True

    expected = {
        "A": {"H2Adj2": "2"},
        "B": {
            "H1Adj2": "3/4",
            "H2Adj2": "3/4",
            "L12Adj2": "3/4",
        },
        "C": {
            "H1Adj2": "3/4",
            "H2Adj2": "3/4",
            "L12Adj2": "3/4",
        },
        "D": {"H1Adj2": "2"},
        "E": {
            "H1Adj2": "2",
            "H2Adj2": "2",
            "L12Adj2": "4",
        },
    }

    print()
    print("Representation-level portal identities")

    for model, (d_s1, d_s2, d_f) in MODELS.items():
        prediction = portal_quartic_beta_group_factors(
            d_s1,
            d_s2,
            d_f,
            0,
        )

        actual = {}
        if "lambdaH1Adj_sq" in prediction.beta_lambdaH1:
            actual["H1Adj2"] = prediction.beta_lambdaH1[
                "lambdaH1Adj_sq"
            ]
        if "lambdaH2Adj_sq" in prediction.beta_lambdaH2:
            actual["H2Adj2"] = prediction.beta_lambdaH2[
                "lambdaH2Adj_sq"
            ]
        if "lambda12Adj_sq" in prediction.beta_lambda12:
            actual["L12Adj2"] = prediction.beta_lambda12[
                "lambda12Adj_sq"
            ]

        match = actual == expected[model]
        overall &= match

        rows.append(
            {
                "model": model,
                "actual": actual,
                "expected": expected[model],
                "match": match,
            }
        )
        print(
            f"T3-{model} adjoint-squared "
            f"{'PASS' if match else 'FAIL'} {actual}"
        )

    e = portal_quartic_beta_group_factors(3, 3, 2, 0)

    cross_expected = {
        "H1_cross": "4",
        "H2_cross": "4",
        "L12_cross_sq": "4",
        "L12_cross_S1A": "4",
        "L12_cross_S2A": "4",
        "Cross_beta_g2_four": "6",
        "Cross_beta_adj_sq": "2",
        "Cross_beta_L12_cross": "8",
        "Cross_beta_cross_sq": "10",
        "Cross_beta_cross_S1A": "-2",
        "Cross_beta_cross_S2A": "-2",
    }

    cross_actual = {
        "H1_cross": e.beta_lambdaH1["lambdaH2*lambda12Cross"],
        "H2_cross": e.beta_lambdaH2["lambdaH1*lambda12Cross"],
        "L12_cross_sq": e.beta_lambda12["lambda12Cross_sq"],
        "L12_cross_S1A": e.beta_lambda12[
            "lambda12Cross*lambdaS1Adj"
        ],
        "L12_cross_S2A": e.beta_lambda12[
            "lambda12Cross*lambdaS2Adj"
        ],
        "Cross_beta_g2_four": e.generated_non_singlet[
            "beta_lambda12Cross:g2_four"
        ],
        "Cross_beta_adj_sq": e.generated_non_singlet[
            "beta_lambda12Cross:lambda12Adj_sq"
        ],
        "Cross_beta_L12_cross": e.generated_non_singlet[
            "beta_lambda12Cross:lambda12*lambda12Cross"
        ],
        "Cross_beta_cross_sq": e.generated_non_singlet[
            "beta_lambda12Cross:lambda12Cross_sq"
        ],
        "Cross_beta_cross_S1A": e.generated_non_singlet[
            "beta_lambda12Cross:lambda12Cross*lambdaS1Adj"
        ],
        "Cross_beta_cross_S2A": e.generated_non_singlet[
            "beta_lambda12Cross:lambda12Cross*lambdaS2Adj"
        ],
    }

    match = cross_actual == cross_expected
    overall &= match
    rows.append(
        {
            "model": "E-cross",
            "actual": cross_actual,
            "expected": cross_expected,
            "match": match,
        }
    )
    print(f"T3-E Cross {'PASS' if match else 'FAIL'}")

    return rows, overall


def check_gauge_scan():
    rows = []
    overall = True

    for model, (d_s1, d_s2, d_f) in MODELS.items():
        for alpha in (-2, -1, 0, 1, 2):
            prediction = portal_quartic_beta_group_factors(
                d_s1,
                d_s2,
                d_f,
                alpha,
            )

            h1_total = F(prediction.beta_lambdaH1["g2_four"])
            h2_total = F(prediction.beta_lambdaH2["g2_four"])

            if d_s1 == 1 and h1_total != 0:
                match = False
            elif d_s2 == 1 and h2_total != 0:
                match = False
            else:
                match = True

            match &= F(prediction.beta_lambdaH1["gY_four"]) >= 0
            match &= F(prediction.beta_lambdaH2["gY_four"]) >= 0
            match &= F(prediction.beta_lambda12["gY_four"]) >= 0

            overall &= match
            rows.append(
                {
                    "model": model,
                    "alpha": alpha,
                    "match": match,
                }
            )

    print(
        "25-point gauge identity scan "
        f"{'PASS' if overall else 'FAIL'}"
    )
    return rows, overall


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--root",
        type=Path,
        default=Path("output/full/hypercharge"),
        help="Saved RGBeta UV-output root. Missing data is reported as skipped.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "output/group_factors/portal_quartic_consolidation.json"
        ),
    )
    parser.add_argument(
        "--backbone-output",
        type=Path,
        default=None,
        help=(
            "Optional compatibility artifact containing the former "
            "ValidatePortalQuarticGroupFactors RGBeta-backbone result."
        ),
    )
    args = parser.parse_args()

    backbone_rows, backbone_ok = check_rgbeta_backbone(args.root)
    if backbone_ok is None:
        print(
            "Resolved portal-quartic backbone vs RGBeta: SKIPPED "
            f"(no saved outputs under {args.root})"
        )

    representation_rows, representation_ok = (
        check_representation_identities()
    )
    gauge_rows, gauge_ok = check_gauge_scan()

    overall = (
        representation_ok
        and gauge_ok
        and (backbone_ok is None or backbone_ok)
    )

    if args.backbone_output is not None:
        backbone_payload = {
            "status": (
                "Skipped"
                if backbone_ok is None
                else ("Success" if backbone_ok else "Failed")
            ),
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
            "checks": backbone_rows,
        }
        args.backbone_output.parent.mkdir(
            parents=True,
            exist_ok=True,
        )
        args.backbone_output.write_text(
            json.dumps(backbone_payload, indent=2),
            encoding="utf-8",
        )

    payload = {
        "status": "Success" if overall else "Failed",
        "rgbeta_backbone": {
            "status": (
                "Skipped"
                if backbone_ok is None
                else ("Success" if backbone_ok else "Failed")
            ),
            "checks": backbone_rows,
        },
        "representation_checks": representation_rows,
        "gauge_checks": gauge_rows,
        "independent_derivation_scripts": [
            "DerivePortalAdjointRecouplings.py",
            "DeriveEPortalCrossScalarRecouplings.py",
            "DerivePortalGaugeQuartics.py",
        ],
        "notes": [
            (
                "Production coefficients depend on SU(2) dimensions/Casimirs "
                "and hypercharges, not A-E labels."
            ),
            (
                "The S1-S2 d=3,d=3 sector keeps the independent Cross channel "
                "explicitly."
            ),
        ],
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2),
        encoding="utf-8",
    )

    print()
    if args.backbone_output is not None:
        print(f"Backbone compatibility JSON: {args.backbone_output}")
    print(f"JSON summary: {args.output}")
    print(f"OVERALL:      {payload['status']}")

    return 0 if overall else 1


if __name__ == "__main__":
    raise SystemExit(main())
