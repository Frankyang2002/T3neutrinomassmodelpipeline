from __future__ import annotations

"""Fermion- and scalar-mass group-factor validation utilities."""

import argparse
from fractions import Fraction
import json
from pathlib import Path
import re
from typing import Any
import sympy as sp
from RGE.group_factors.RepresentationFactors import su2_quadratic_casimir_from_dimension
from RGE.group_factors.RepresentationFactors import canonical_yukawa_leg_factors
from dataclasses import dataclass

# ---------------------------------------------------------------------------
# Fermion-mass group-factor validation
# Former source: ValidateFermionMassGroupFactors.py
# ---------------------------------------------------------------------------

"""Validate the generic T3 heavy-fermion mass beta function against RGBeta.

The analytic vector-like beta_MF coefficients are small enough to live directly
in this validator.  They are expressed only through the shared representation
helpers consolidated in ``RepresentationFactors.py``:

    GF_i = max(dF,dSi)/dF

    16*pi^2 beta_MF =
        1/2 GF1 (y1^T y1*) MF
      + 1/2 GF2 MF (y2^dagger y2)
      - 6 g2^2 C2(F) MF
      - 6 gY^2 YF^2 MF,

with
    C2(F) = (dF^2-1)/4,
    YF    = (alpha+1)/2.

This validator parses RGBeta's saved ``report_betas["MF"]`` structurally by
replacing each complete ``Matrix[...][heavy[i],heavy[j]]`` object with a
symbolic token and letting SymPy extract exact coefficients.  This also handles
the neutral-F branch, where RGBeta writes both MF and Trans[MF].

The neutral-branch Yukawa tensor structure is deliberately not identified with
the vector-like left/right formula here.
"""





MODEL_DIMS = {
    "A": (1, 3, 2),
    "B": (2, 2, 1),
    "C": (2, 2, 3),
    "D": (3, 1, 2),
    "E": (3, 3, 2),
}


def _txt(value: Fraction | None) -> str | None:
    if value is None:
        return None
    return (
        str(value.numerator)
        if value.denominator == 1
        else f"{value.numerator}/{value.denominator}"
    )


def expected_fermion_mass_coefficients(
    dS1: int,
    dS2: int,
    dF: int,
    alpha: int,
) -> dict[str, Fraction]:
    dS1, dS2, dF, alpha = map(int, (dS1, dS2, dF, alpha))

    if abs(dS1 - dF) != 1 or abs(dS2 - dF) != 1:
        raise ValueError("T3 Yukawa invariance requires dSi=dF+/-1.")

    y1legs = canonical_yukawa_leg_factors(dF, dS1)
    y2legs = canonical_yukawa_leg_factors(dF, dS2)

    c2f = su2_quadratic_casimir_from_dimension(dF)
    yf = Fraction(alpha + 1, 2)

    return {
        "y1_left": Fraction(1, 2) * y1legs.G_heavy,
        "y2_right": Fraction(1, 2) * y2legs.G_heavy,
        "su2_gauge": -6 * c2f,
        "u1_gauge": -6 * yf * yf,
    }


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
    for index in range(open_index, len(text)):
        if text[index] == "[":
            depth += 1
        elif text[index] == "]":
            depth -= 1
            if depth == 0:
                return index

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

    for number, (start, end, inner) in enumerate(objects):
        token = classifier(inner, number)
        if token is None:
            token = f"Z{number}"

        replacements.append((start, end, token))
        labels[token] = inner

    result = expr
    for start, end, token in reversed(replacements):
        result = result[:start] + token + result[end:]

    return result, labels


def _sympy_scalar_expr(
    text: str,
    extra_symbols: tuple[str, ...],
) -> sp.Expr:
    cleaned = text.replace("^", "**")
    locals_ = {
        name: sp.Symbol(name)
        for name in extra_symbols
    }
    locals_.update(
        {
            "g2": sp.Symbol("g2"),
            "gY": sp.Symbol("gY"),
        }
    )

    for name in set(re.findall(r"\bZ\d+\b", cleaned)):
        locals_[name] = sp.Symbol(name)

    return sp.expand(sp.sympify(cleaned, locals=locals_))


def _gauge_coeff(expr: str, gauge: str) -> Fraction:
    """Coefficient of g^2 times the physical MF direction.

    The neutral-F branch contains Matrix[MF] + Matrix[Trans[MF]].
    Both are mapped onto the same physical mass token M, i.e. MF^T=MF.
    """

    def classify(inner: str, _number: int) -> str | None:
        compact = inner.replace(" ", "")
        if compact in {"MF", "Trans[MF]"}:
            return "M"
        return None

    reduced, _ = _replace_matrix_objects(expr, classify)
    parsed = _sympy_scalar_expr(reduced, ("M",))

    gauge_symbol = sp.Symbol(gauge)
    mass_symbol = sp.Symbol("M")

    coefficient = (
        sp.expand(parsed)
        .coeff(mass_symbol, 1)
        .coeff(gauge_symbol, 2)
    )

    for symbol in coefficient.free_symbols:
        if symbol != gauge_symbol:
            coefficient = coefficient.subs(symbol, 0)

    coefficient = sp.simplify(coefficient)

    return Fraction(
        int(sp.numer(coefficient)),
        int(sp.denom(coefficient)),
    )


def _yukawa_coeff(expr: str, target_inner: str) -> Fraction:
    """Coefficient of one exact matrix structure in beta_MF."""

    def classify(inner: str, _number: int) -> str | None:
        if inner.replace(" ", "") == target_inner:
            return "X"
        return None

    reduced, _ = _replace_matrix_objects(expr, classify)
    parsed = _sympy_scalar_expr(reduced, ("X",))
    target = sp.Symbol("X")

    coefficient = sp.expand(parsed).coeff(target, 1)

    for symbol in list(coefficient.free_symbols):
        if symbol != target:
            coefficient = coefficient.subs(symbol, 0)

    coefficient = sp.simplify(coefficient)

    return Fraction(
        int(sp.numer(coefficient)),
        int(sp.denom(coefficient)),
    )


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


def a_main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Validate T3 beta_MF group factors against saved RGBeta results."
        )
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path("output/full/hypercharge"),
    )
    parser.add_argument(
        "--json-output",
        type=Path,
        default=Path(
            "output/group_factors/fermion_mass_validation.json"
        ),
    )
    args = parser.parse_args()

    paths = sorted(
        args.root.glob("T3_*_alpha_*/data/uv_rgbeta_rge.json")
    )
    if not paths:
        raise FileNotFoundError(
            f"No RGBeta files found below {args.root}"
        )

    checks: list[dict[str, Any]] = []
    all_ok = True

    print("Gauge validation across hypercharge scan")
    print("model alpha  gY(exp/rgb/res)   g2(exp/rgb/res)  status")
    print("-" * 70)

    for path in paths:
        model, alpha = _model_alpha(path)
        dS1, dS2, dF = MODEL_DIMS[model]

        expected = expected_fermion_mass_coefficients(
            dS1,
            dS2,
            dF,
            alpha,
        )
        expr = _mf_report(_load(path))

        got_u1 = _gauge_coeff(expr, "gY")
        got_su2 = _gauge_coeff(expr, "g2")
        exp_u1 = expected["u1_gauge"]
        exp_su2 = expected["su2_gauge"]

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

        checks.append(
            {
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
            }
        )

    print()
    print("Vector-like Yukawa-leg validation at alpha=0")
    print("model  y1L(exp/rgb/res)  y2R(exp/rgb/res)  status")
    print("-" * 62)

    for model, (dS1, dS2, dF) in MODEL_DIMS.items():
        path = (
            args.root
            / f"T3_{model}_alpha_p0"
            / "data"
            / "uv_rgbeta_rge.json"
        )
        expr = _mf_report(_load(path))

        expected = expected_fermion_mass_coefficients(
            dS1,
            dS2,
            dF,
            0,
        )

        got1 = _vectorlike_y1_coeff(expr)
        got2 = _vectorlike_y2_coeff(expr)
        exp1 = expected["y1_left"]
        exp2 = expected["y2_right"]

        residual1 = got1 - exp1
        residual2 = got2 - exp2
        ok = residual1 == 0 and residual2 == 0
        all_ok &= ok

        print(
            f"T3-{model:1}    "
            f"{_txt(exp1):>4}/{_txt(got1):>4}/{_txt(residual1):>3}       "
            f"{_txt(exp2):>4}/{_txt(got2):>4}/{_txt(residual2):>3}      "
            f"{'PASS' if ok else 'FAIL'}"
        )

        checks.append(
            {
                "kind": "vectorlike_yukawa",
                "model": model,
                "alpha": 0,
                "expected_y1_left": _txt(exp1),
                "rgbeta_y1_left": _txt(got1),
                "y1_residual": _txt(residual1),
                "expected_y2_right": _txt(exp2),
                "rgbeta_y2_right": _txt(got2),
                "y2_residual": _txt(residual2),
                "match": ok,
                "source": str(path),
            }
        )

    payload = {
        "status": "Success" if all_ok else "Failed",
        "formula_vectorlike_branch": (
            "16*pi^2 beta_MF = 1/2 GF1 (y1^T y1*) MF "
            "+ 1/2 GF2 MF (y2^dagger y2) "
            "- 6 g2^2 C2(F) MF - 6 gY^2 YF^2 MF"
        ),
        "neutral_branch_note": (
            "For YF=0 RGBeta defines the mass as {F,F}; its saved beta "
            "contains MF and Trans[MF]. Gauge validation projects onto "
            "MF^T=MF. The neutral-branch Yukawa tensor structure is not "
            "identified with the vector-like left/right formula in this "
            "validator."
        ),
        "checks": checks,
    }

    args.json_output.parent.mkdir(parents=True, exist_ok=True)
    args.json_output.write_text(
        json.dumps(payload, indent=2),
        encoding="utf-8",
    )

    print()
    print(f"RGBeta files checked: {len(paths)}")
    print(f"JSON summary:         {args.json_output}")
    print(f"OVERALL:              {payload['status']}")

    return 0 if all_ok else 1


# ---------------------------------------------------------------------------
# Scalar-mass group-factor validation
# Former source: ValidateScalarMassGroupFactors.py
# ---------------------------------------------------------------------------

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


def b_main() -> int:
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
