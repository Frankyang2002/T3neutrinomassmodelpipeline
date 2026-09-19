from __future__ import annotations

"""Validate the physical T3 scalar-mass group factors against RGBeta.

The analytic scalar-mass formulas that previously lived in
``ScalarMassGroupFactors.py`` are kept directly in this validator so the known
neutral/self-conjugate branch issue remains explicit while reducing one source
file.

Physical branch policy for alpha=-1
-----------------------------------
* dF odd (B,C):
    RGBeta and the physical T3 model are both on the self-conjugate branch.
    Validate the -16 G_Si heavy-mass insertion.

* dF even (A,D,E):
    RGBeta takes a neutral/self-conjugate branch that is not the same field
    content as the current Matchete T3 model.  All common structures are still
    validated, but the heavy-mass insertion is recorded as a branch-mismatch
    diagnostic rather than counted as a physics failure.

No coefficient is changed in this consolidation.  In particular, the possible
factor-4 self-conjugate scalar-mass issue remains visible exactly as before.
"""

import argparse
from dataclasses import dataclass
from fractions import Fraction
import json
from pathlib import Path
import re

import sympy as sp

from RGE.group_factors.RepresentationFactors import (
    su2_quadratic_casimir_from_dimension,
)


MODELS = {
    "A": (1, 3, 2),
    "B": (2, 2, 1),
    "C": (2, 2, 3),
    "D": (3, 1, 2),
    "E": (3, 3, 2),
}


def txt(value: Fraction) -> str:
    return (
        str(value.numerator)
        if value.denominator == 1
        else f"{value.numerator}/{value.denominator}"
    )


def F(value) -> Fraction:
    return Fraction(str(value))


def physical_self_conjugate_f(dF: int, alpha: int) -> bool:
    """Current Matchete-side physical self-conjugate-F criterion."""
    return int(alpha) == -1 and int(dF) % 2 == 1


def _validate_t3_dimensions(d1: int, d2: int, dF: int) -> None:
    if min(d1, d2, dF) < 1:
        raise ValueError("Representation dimensions must be positive.")

    if d1 not in (1, 2, 3) or d2 not in (1, 2, 3) or dF not in (1, 2, 3):
        raise ValueError("Current T3 implementation supports d in {1,2,3}.")

    if abs(d1 - dF) != 1 or abs(d2 - dF) != 1:
        raise ValueError("T3 Yukawa invariance requires dSi=dF+/-1.")


@dataclass(frozen=True)
class ScalarMassGroupFactors:
    dS1: int
    dS2: int
    dF: int
    alpha: int
    self_conjugate_F: bool
    GS1: str
    GS2: str
    C2S1: str
    C2S2: str
    YS1: str
    YS2: str
    beta_mS1Sq: dict[str, str]
    beta_mS2Sq: dict[str, str]


def scalar_mass_group_factors(
    dS1: int,
    dS2: int,
    dF: int,
    alpha: int,
    *,
    self_conjugate_F: bool | None = None,
) -> ScalarMassGroupFactors:
    """Return the preserved physical scalar-mass group factors."""

    d1, d2, dF, alpha = map(int, (dS1, dS2, dF, alpha))
    _validate_t3_dimensions(d1, d2, dF)

    physical_sc = physical_self_conjugate_f(dF, alpha)

    if self_conjugate_F is None:
        self_conjugate_F = physical_sc
    else:
        self_conjugate_F = bool(self_conjugate_F)

    if self_conjugate_F and alpha != -1:
        raise ValueError(
            "A self-conjugate F requires Y_F=0, hence alpha=-1."
        )

    if self_conjugate_F and dF % 2 == 0:
        raise ValueError(
            "Current physical T3 branch does not treat even-dimensional "
            "SU(2) F as self-conjugate."
        )

    C1 = su2_quadratic_casimir_from_dimension(d1)
    C2 = su2_quadratic_casimir_from_dimension(d2)

    Y1 = Fraction(alpha, 2)
    Y2 = Fraction(alpha + 2, 2)

    GS1 = Fraction(max(dF, d1), d1)
    GS2 = Fraction(max(dF, d2), d2)

    # Preserve the existing project formula exactly.
    heavy_factor = Fraction(-16 if self_conjugate_F else -4)

    m1 = {
        "mS1Sq*Tr_y1": txt(2 * GS1),
        "Tr_MF_y1": txt(heavy_factor * GS1),
        "lambdaS1*mS1Sq": txt(Fraction(2 * (d1 + 1))),
        "g2_sq*mS1Sq": txt(-6 * C1),
        "gY_sq*mS1Sq": txt(-6 * Y1 * Y1),
        "lambda12*mS2Sq": txt(Fraction(2 * d2)),
    }

    if d1 == 3:
        m1["lambdaS1Adj*mS1Sq"] = txt(2 * C1)

    if d1 == d2 == 3:
        m1["lambda12Cross*mS2Sq"] = "4"

    m2 = {
        "mS2Sq*Tr_y2": txt(2 * GS2),
        "Tr_MF_y2": txt(heavy_factor * GS2),
        "lambdaS2*mS2Sq": txt(Fraction(2 * (d2 + 1))),
        "g2_sq*mS2Sq": txt(-6 * C2),
        "gY_sq*mS2Sq": txt(-6 * Y2 * Y2),
        "lambda12*mS1Sq": txt(Fraction(2 * d1)),
    }

    if d2 == 3:
        m2["lambdaS2Adj*mS2Sq"] = txt(2 * C2)

    if d1 == d2 == 3:
        m2["lambda12Cross*mS1Sq"] = "4"

    return ScalarMassGroupFactors(
        dS1=d1,
        dS2=d2,
        dF=dF,
        alpha=alpha,
        self_conjugate_F=self_conjugate_F,
        GS1=txt(GS1),
        GS2=txt(GS2),
        C2S1=txt(C1),
        C2S2=txt(C2),
        YS1=txt(Y1),
        YS2=txt(Y2),
        beta_mS1Sq=m1,
        beta_mS2Sq=m2,
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

    raise ValueError("Unbalanced []")


def replace_traces(expr: str) -> str:
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

        if "MF" in inner and "y1" in inner:
            token = "TMy1"
        elif "MF" in inner and "y2" in inner:
            token = "TMy2"
        elif "y1" in inner:
            token = "Ty1"
        elif "y2" in inner:
            token = "Ty2"
        else:
            token = f"Tother{unknown}"
            unknown += 1

        replacements.append((start, close + 1, token))
        pos = close + 1

    for start, end, token in reversed(replacements):
        result = result[:start] + token + result[end:]

    return result


def replace_bar_calls(expr: str) -> str:
    result = expr
    count = 0

    while "Bar[" in result:
        start = result.find("Bar[")
        open_index = start + 3
        close = matching_bracket(result, open_index)
        inner = result[open_index + 1:close]
        safe = re.sub(r"[^A-Za-z0-9_]", "_", inner)

        result = (
            result[:start]
            + f"Bar_{safe}_{count}"
            + result[close + 1:]
        )
        count += 1

    return result


def parse_beta(raw: str) -> sp.Expr:
    text = replace_traces(str(raw))
    text = replace_bar_calls(text)
    text = text.replace("^", "**")

    names = set(
        re.findall(
            r"\b[A-Za-z][A-Za-z0-9_]*\b",
            text,
        )
    )
    locals_ = {
        name: sp.Symbol(name)
        for name in names
    }

    return sp.expand(
        sp.sympify(
            text,
            locals=locals_,
        )
    )


def coeff(
    expr: sp.Expr,
    monomial: dict[str, int],
) -> Fraction:
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

    return Fraction(
        int(sp.numer(out)),
        int(sp.denom(out)),
    )


def model_alpha(path: Path) -> tuple[str, int]:
    name = path.parents[1].name
    match = re.fullmatch(
        r"T3_([A-E])_alpha_([mp])(\d+)",
        name,
    )
    if not match:
        raise ValueError(name)

    model, sign, magnitude = match.groups()
    alpha = int(magnitude)

    return model, (-alpha if sign == "m" else alpha)


def load_betas(path: Path) -> dict[str, sp.Expr]:
    payload = json.loads(
        path.read_text(
            encoding="utf-8-sig",
        )
    )

    if payload.get("status") != "Success":
        raise ValueError(
            f"RGBeta status is not Success: {path}"
        )

    report = payload.get("report_betas", {}) or {}

    return {
        "mS1Sq": parse_beta(report["mS1Sq"]),
        "mS2Sq": parse_beta(report["mS2Sq"]),
    }


def checks_for(
    prediction: ScalarMassGroupFactors,
) -> dict[str, dict[str, tuple[dict[str, int], str]]]:
    m1 = prediction.beta_mS1Sq
    m2 = prediction.beta_mS2Sq

    c1 = {
        "mS1Sq*Tr_y1": (
            {"mS1Sq": 1, "Ty1": 1},
            m1["mS1Sq*Tr_y1"],
        ),
        "Tr_MF_y1": (
            {"TMy1": 1},
            m1["Tr_MF_y1"],
        ),
        "lambdaS1*mS1Sq": (
            {"lambdaS1": 1, "mS1Sq": 1},
            m1["lambdaS1*mS1Sq"],
        ),
        "g2_sq*mS1Sq": (
            {"g2": 2, "mS1Sq": 1},
            m1["g2_sq*mS1Sq"],
        ),
        "gY_sq*mS1Sq": (
            {"gY": 2, "mS1Sq": 1},
            m1["gY_sq*mS1Sq"],
        ),
        "lambda12*mS2Sq": (
            {"lambda12": 1, "mS2Sq": 1},
            m1["lambda12*mS2Sq"],
        ),
    }

    if "lambdaS1Adj*mS1Sq" in m1:
        c1["lambdaS1Adj*mS1Sq"] = (
            {"lambdaS1Adj": 1, "mS1Sq": 1},
            m1["lambdaS1Adj*mS1Sq"],
        )

    if "lambda12Cross*mS2Sq" in m1:
        c1["lambda12Cross*mS2Sq"] = (
            {"lambda12Cross": 1, "mS2Sq": 1},
            m1["lambda12Cross*mS2Sq"],
        )

    c2 = {
        "mS2Sq*Tr_y2": (
            {"mS2Sq": 1, "Ty2": 1},
            m2["mS2Sq*Tr_y2"],
        ),
        "Tr_MF_y2": (
            {"TMy2": 1},
            m2["Tr_MF_y2"],
        ),
        "lambdaS2*mS2Sq": (
            {"lambdaS2": 1, "mS2Sq": 1},
            m2["lambdaS2*mS2Sq"],
        ),
        "g2_sq*mS2Sq": (
            {"g2": 2, "mS2Sq": 1},
            m2["g2_sq*mS2Sq"],
        ),
        "gY_sq*mS2Sq": (
            {"gY": 2, "mS2Sq": 1},
            m2["gY_sq*mS2Sq"],
        ),
        "lambda12*mS1Sq": (
            {"lambda12": 1, "mS1Sq": 1},
            m2["lambda12*mS1Sq"],
        ),
    }

    if "lambdaS2Adj*mS2Sq" in m2:
        c2["lambdaS2Adj*mS2Sq"] = (
            {"lambdaS2Adj": 1, "mS2Sq": 1},
            m2["lambdaS2Adj*mS2Sq"],
        )

    if "lambda12Cross*mS1Sq" in m2:
        c2["lambda12Cross*mS1Sq"] = (
            {"lambda12Cross": 1, "mS1Sq": 1},
            m2["lambda12Cross*mS1Sq"],
        )

    return {
        "mS1Sq": c1,
        "mS2Sq": c2,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(
            "output/full/hypercharge"
        ),
    )
    parser.add_argument(
        "--json-output",
        type=Path,
        default=Path(
            "output/group_factors/scalar_mass_validation.json"
        ),
    )
    args = parser.parse_args()

    paths = sorted(
        args.root.glob(
            "T3_*_alpha_*/data/uv_rgbeta_rge.json"
        )
    )

    if not paths:
        raise SystemExit(
            f"No RGBeta outputs found below {args.root}"
        )

    rows = []
    branch_diagnostics = []
    overall = True

    print(
        "T3 heavy-scalar mass group-factor validation"
    )
    print(
        "model alpha  beta    checks  status"
    )
    print("-" * 46)

    for path in paths:
        model, alpha = model_alpha(path)
        d1, d2, dF = MODELS[model]

        prediction = scalar_mass_group_factors(
            d1,
            d2,
            dF,
            alpha,
        )
        betas = load_betas(path)

        rgbeta_even_neutral_mismatch = (
            alpha == -1
            and dF % 2 == 0
        )

        for beta_name, checks in checks_for(
            prediction
        ).items():
            beta = betas[beta_name]
            details = {}
            row_ok = True
            validated_count = 0

            for label, (
                monomial,
                expected_text,
            ) in checks.items():
                got = coeff(
                    beta,
                    monomial,
                )
                expected = F(
                    expected_text
                )
                residual = (
                    got - expected
                )

                excluded_branch_mismatch = (
                    rgbeta_even_neutral_mismatch
                    and label.startswith(
                        "Tr_MF_"
                    )
                )

                if excluded_branch_mismatch:
                    branch_diagnostics.append(
                        {
                            "model": model,
                            "alpha": alpha,
                            "dF": dF,
                            "beta": beta_name,
                            "structure": label,
                            "physical_expected": txt(
                                expected
                            ),
                            "rgbeta_neutral_branch": txt(
                                got
                            ),
                            "residual": txt(
                                residual
                            ),
                            "reason": (
                                "RGBeta uses its YF=0 neutral/"
                                "self-conjugate branch for even "
                                "dF=2, while the current physical "
                                "Matchete T3 model keeps this "
                                "pseudoreal representation "
                                "non-self-conjugate."
                            ),
                        }
                    )

                    details[label] = {
                        "expected": txt(
                            expected
                        ),
                        "rgbeta": txt(
                            got
                        ),
                        "residual": txt(
                            residual
                        ),
                        "match": None,
                        "status": (
                            "ExcludedBranchMismatch"
                        ),
                    }
                    continue

                validated_count += 1
                passed = residual == 0
                row_ok &= passed

                details[label] = {
                    "expected": txt(
                        expected
                    ),
                    "rgbeta": txt(
                        got
                    ),
                    "residual": txt(
                        residual
                    ),
                    "match": passed,
                    "status": "Validated",
                }

            overall &= row_ok

            suffix = (
                " + branch diagnostic"
                if rgbeta_even_neutral_mismatch
                else ""
            )

            print(
                f"T3-{model} {alpha:>5}  "
                f"{beta_name:<6} "
                f"{validated_count:>6}  "
                f"{'PASS' if row_ok else 'FAIL'}"
                f"{suffix}"
            )

            rows.append(
                {
                    "model": model,
                    "alpha": alpha,
                    "dF": dF,
                    "physical_self_conjugate_F": (
                        physical_self_conjugate_f(
                            dF,
                            alpha,
                        )
                    ),
                    "beta": beta_name,
                    "source": str(
                        path
                    ),
                    "checks": details,
                    "match": row_ok,
                }
            )

    payload = {
        "status": (
            "Success"
            if overall
            else "Failed"
        ),
        "physical_formulas": {
            "mSi_Tr_yi": "2 G_Si",
            "heavy_mass_yukawa_trace_vectorlike": (
                "-4 G_Si"
            ),
            "heavy_mass_yukawa_trace_self_conjugate": (
                "-16 G_Si"
            ),
            "self_conjugate_branch": (
                "alpha=-1 and dF odd"
            ),
            "lambdaSi_mSi": "2(dSi+1)",
            "lambdaSiAdj_mSi_triplet": (
                "2 C2(Si)=4"
            ),
            "g2_sq_mSi": "-6 C2(Si)",
            "gY_sq_mSi": "-6 Yi^2",
            "lambda12_other_mass": (
                "2 d_other"
            ),
            "T3E_cross_other_mass": "4",
        },
        "known_issue": (
            "The possible factor-4 self-conjugate scalar-mass "
            "beta issue remains unresolved. This validator "
            "preserves the existing -16 G_Si physical branch "
            "formula and does not alter it."
        ),
        "branch_diagnostics": (
            branch_diagnostics
        ),
        "checks": rows,
    }

    args.json_output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )
    args.json_output.write_text(
        json.dumps(
            payload,
            indent=2,
        ),
        encoding="utf-8",
    )

    print()
    print(
        "Even-dF alpha=-1 RGBeta "
        "heavy-mass diagnostics excluded: "
        f"{len(branch_diagnostics)}"
    )
    print(
        f"JSON summary: "
        f"{args.json_output}"
    )
    print(
        f"OVERALL:      "
        f"{payload['status']}"
    )

    return 0 if overall else 1


if __name__ == "__main__":
    raise SystemExit(main())
