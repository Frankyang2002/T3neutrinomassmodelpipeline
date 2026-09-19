from __future__ import annotations

"""Independent SU(2) recoupling check for beta(lambdaT3).

This script does NOT fit the A--E RGBeta numbers.

It constructs the unique T3 invariant in a canonical SU(2) spherical basis:

  (H H)_{J=1} (S1 \tilde S2^*)_{J=1} -> J=0,

with \tilde S2^* = E_j S2^* the standard SU(2) dual-to-primal map, and then
evaluates the pure-scalar one-loop real-tensor contraction with the physical
adjoint quartics

  (H^\dagger T^A H)(S_i^\dagger T^A S_i),
  (S1^\dagger T^A S1)(S2^\dagger T^A S2).

The real-scalar convention is

  V4 = (1/4!) lambda_abcd phi_a phi_b phi_c phi_d,

for which the scalar-scalar one-loop term is the sum of the three pair
contractions lambda_abef lambda_efcd + permutations.

The overall normalization of the T3 Clebsch tensor cancels in the projection,
so the result is a representation-theory prediction, independent of the
component basis and CG normalization.

If an RGBeta UV JSON is supplied, the script also extracts the corresponding
lambdaT3 * quartic coefficients directly from beta_lambdaT3 and compares them.
A disagreement is therefore a convention/implementation diagnostic, not a fit.
"""

import argparse
from collections import Counter
from itertools import combinations_with_replacement
import json
from math import factorial
from pathlib import Path
import re
import sys

import sympy as sp
from sympy.physics.wigner import clebsch_gordan

PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from RGE.running.eft1.MatcheteParsing import matching_bracket


TARGETS = ("lambdaH1Adj", "lambdaH2Adj", "lambda12Adj")


def _m_values(d: int) -> list[sp.Rational]:
    j = sp.Rational(d - 1, 2)
    return [j - i for i in range(d)]


def _su2_generators(d: int) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    j = sp.Rational(d - 1, 2)
    ms = _m_values(d)
    jp = sp.zeros(d)
    jm = sp.zeros(d)

    for col, m in enumerate(ms):
        mp = m + 1
        if mp in ms:
            row = ms.index(mp)
            jp[row, col] = sp.sqrt((j - m) * (j + m + 1))
        mm = m - 1
        if mm in ms:
            row = ms.index(mm)
            jm[row, col] = sp.sqrt((j + m) * (j - m + 1))

    return (
        sp.simplify((jp + jm) / 2),
        sp.simplify((jp - jm) / (2 * sp.I)),
        sp.diag(*ms),
    )


def _dual_map(d: int) -> sp.Matrix:
    """E_j with E |j,m>^* = (-1)^(j-m) |j,-m>."""
    j = sp.Rational(d - 1, 2)
    ms = _m_values(d)
    E = sp.zeros(d)
    for row, m in enumerate(ms):
        col = ms.index(-m)
        E[row, col] = (-1) ** int(j - m)
    return E


def _complex_block(d: int, first: int):
    z = []
    zb = []
    real_symbols = []
    for c in range(d):
        x = sp.Symbol(f"x{first + 2*c}")
        y = sp.Symbol(f"x{first + 2*c + 1}")
        z.append((x + sp.I*y) / sp.sqrt(2))
        zb.append((x - sp.I*y) / sp.sqrt(2))
        real_symbols.extend((x, y))
    return z, zb, real_symbols


def _t3_polynomial(d1: int, d2: int):
    H, Hb, hv = _complex_block(2, 1)
    S1, S1b, s1v = _complex_block(d1, 5)
    S2, S2b, s2v = _complex_block(d2, 5 + 2*d1)
    variables = hv + s1v + s2v

    mh = _m_values(2)
    m1 = _m_values(d1)
    m2 = _m_values(d2)
    E2 = _dual_map(d2)

    expr = sp.S.Zero
    for i, mi in enumerate(mh):
        for j, mj in enumerate(mh):
            for M in (1, 0, -1):
                c_h = clebsch_gordan(
                    sp.Rational(1, 2), sp.Rational(1, 2), 1,
                    mi, mj, M,
                )
                if c_h == 0:
                    continue

                for b, mb in enumerate(m1):
                    for r, mr in enumerate(m2):
                        for N in (1, 0, -1):
                            c_s = clebsch_gordan(
                                sp.Rational(d1 - 1, 2),
                                sp.Rational(d2 - 1, 2),
                                1,
                                mb, mr, N,
                            )
                            if c_s == 0:
                                continue
                            c_0 = clebsch_gordan(1, 1, 0, M, N, 0)
                            if c_0 == 0:
                                continue

                            for c in range(d2):
                                if E2[r, c] != 0:
                                    expr += (
                                        c_h * c_s * c_0 * E2[r, c]
                                        * H[i] * H[j] * S1[b] * S2b[c]
                                    )

    return variables, sp.expand(expr), (H, Hb, S1, S1b, S2, S2b)


def _bilinear(zb, z, matrix: sp.Matrix):
    return sp.expand(sum(
        zb[i] * matrix[i, j] * z[j]
        for i in range(len(z))
        for j in range(len(z))
    ))


def _adjoint_polynomials(d1: int, d2: int, fields):
    H, Hb, S1, S1b, S2, S2b = fields
    gh = _su2_generators(2)
    out = {}

    if d1 > 1:
        g1 = _su2_generators(d1)
        out["lambdaH1Adj"] = sp.expand(sum(
            _bilinear(Hb, H, gh[a]) * _bilinear(S1b, S1, g1[a])
            for a in range(3)
        ))

    if d2 > 1:
        g2 = _su2_generators(d2)
        out["lambdaH2Adj"] = sp.expand(sum(
            _bilinear(Hb, H, gh[a]) * _bilinear(S2b, S2, g2[a])
            for a in range(3)
        ))

    if d1 > 1 and d2 > 1:
        g1 = _su2_generators(d1)
        g2 = _su2_generators(d2)
        out["lambda12Adj"] = sp.expand(sum(
            _bilinear(S1b, S1, g1[a]) * _bilinear(S2b, S2, g2[a])
            for a in range(3)
        ))

    return out


def _tensor_from_polynomial(expr: sp.Expr, variables: list[sp.Symbol]):
    poly = sp.Poly(sp.expand(expr), *variables)
    entries = {}

    for powers, coefficient in poly.terms():
        if sum(powers) != 4:
            continue
        key = []
        multiplicity = 1
        for index, power in enumerate(powers):
            key.extend([index + 1] * power)
            multiplicity *= factorial(power)
        value = sp.simplify(coefficient * multiplicity)
        if value != 0:
            entries[tuple(sorted(key))] = value

    return entries


def _get(tensor, *indices):
    return tensor.get(tuple(sorted(indices)), sp.S.Zero)


def _ordered_weight(key):
    counts = Counter(key)
    value = factorial(4)
    for count in counts.values():
        value //= factorial(count)
    return value


def _inner(left, right):
    total = sp.S.Zero
    for key in set(left) | set(right):
        total += (
            _ordered_weight(key)
            * sp.conjugate(_get(left, *key))
            * _get(right, *key)
        )
    return sp.simplify(total)


def _cross_beta(target, other, n_real: int):
    generated = {}

    for key in combinations_with_replacement(range(1, n_real + 1), 4):
        a, b, c, d = key
        value = sp.S.Zero

        for (p, q), (r, s) in (
            ((a, b), (c, d)),
            ((a, c), (b, d)),
            ((a, d), (b, c)),
        ):
            for e in range(1, n_real + 1):
                for f in range(1, n_real + 1):
                    value += (
                        _get(target, p, q, e, f) * _get(other, e, f, r, s)
                        + _get(other, p, q, e, f) * _get(target, e, f, r, s)
                    )

        value = sp.simplify(value)
        if value != 0:
            generated[key] = value

    return generated


def _project(generated, target):
    norm = sp.simplify(_inner(target, target))
    coefficient = sp.simplify(_inner(target, generated) / norm)

    residual = 0
    for key in set(generated) | set(target):
        if sp.simplify(
            _get(generated, *key) - coefficient * _get(target, *key)
        ) != 0:
            residual += 1

    return coefficient, residual


def canonical_recouplings(d1: int, d2: int):
    variables, mix, fields = _t3_polynomial(d1, d2)
    target = _tensor_from_polynomial(mix, variables)
    operators = _adjoint_polynomials(d1, d2, fields)

    rows = {}
    for name, polynomial in operators.items():
        other = _tensor_from_polynomial(polynomial, variables)
        generated = _cross_beta(target, other, len(variables))
        coefficient, residual = _project(generated, target)
        rows[name] = {
            "coefficient": sp.sstr(sp.factor(coefficient)),
            "residual_nonzero_components": residual,
            "closure": residual == 0,
        }

    return rows


def _replace_balanced_calls(text: str, head: str, replacement: str = "0") -> str:
    """Replace head[...] calls, including nested brackets, by a scalar token."""
    needle = head + "["
    while True:
        start = text.find(needle)
        if start < 0:
            break

        open_index = start + len(head)
        try:
            close_index = matching_bracket(text, open_index)
        except ValueError as exc:
            raise ValueError(
                f"Unbalanced {head}[...] call in RGBeta expression."
            ) from exc

        text = text[:start] + replacement + text[close_index + 1:]

    return text


def _numeric_cross_coefficient(beta: str, quartic: str):
    """Extract coeff(lambdaT3*quartic) from RGBeta InputForm.

    RGBeta commonly factors lambdaT3 outside a parenthesis, e.g.

        lambdaT3*(lambda12Adj/6 + 2*lambdaH1Adj/3 + ...)

    so term-by-term string matching is insufficient.  Remove trace/function
    structures irrelevant to the scalar-quartic coefficient, zero all other
    symbols, then let SymPy expand the complete expression.
    """
    if not beta or "lambdaT3" not in beta or quartic not in beta:
        return None

    reduced = beta

    # Trace terms can contain Mathematica Dot/Trans/Bar syntax.  They are
    # independent of the quartic coefficient, so remove them structurally.
    reduced = _replace_balanced_calls(reduced, "Tr", "0")
    reduced = _replace_balanced_calls(reduced, "Trans", "0")
    reduced = _replace_balanced_calls(reduced, "Bar", "0")

    # Mathematica InputForm arithmetic -> SymPy arithmetic.
    reduced = reduced.replace("^", "**")
    reduced = reduced.replace(".", "*")

    # Keep only the two symbols whose mixed coefficient is wanted.
    keep = {"lambdaT3", quartic}
    identifiers = set(re.findall(r"\b[A-Za-z]\w*\b", reduced))
    for name in sorted(identifiers, key=len, reverse=True):
        if name not in keep:
            reduced = re.sub(rf"\b{re.escape(name)}\b", "0", reduced)

    try:
        expr = sp.expand(sp.sympify(reduced))
    except Exception as exc:
        raise ValueError(
            f"Could not parse reduced RGBeta beta(lambdaT3): {reduced}"
        ) from exc

    lt3 = sp.Symbol("lambdaT3")
    q = sp.Symbol(quartic)
    coefficient = sp.simplify(expr.coeff(lt3).coeff(q))
    return coefficient if coefficient != 0 else None


def rgbeta_coefficients(path: Path):
    payload = json.loads(path.read_text(encoding="utf-8-sig"))
    betas = payload.get("betas", {})
    beta = (
        betas.get("lambdaT3")
        or payload.get("report_betas", {}).get("lambdaT3")
        or ""
    )
    return {
        name: _numeric_cross_coefficient(beta, name)
        for name in TARGETS
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("dS1", type=int)
    parser.add_argument("dS2", type=int)
    parser.add_argument("--rgbeta-json", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    rows = canonical_recouplings(args.dS1, args.dS2)
    observed = (
        rgbeta_coefficients(args.rgbeta_json)
        if args.rgbeta_json is not None
        else {}
    )

    print(f"dS1={args.dS1}, dS2={args.dS2}")
    print("Canonical SU(2) recoupling, no A--E coefficient fit:")
    print(f"{'quartic':<18} {'canonical':>12} {'residual':>10} {'RGBeta':>12} {'compare':>10}")
    print("-" * 70)

    comparisons = {}
    for name in TARGETS:
        if name not in rows:
            continue
        canonical = sp.sympify(rows[name]["coefficient"])
        rgb = observed.get(name)
        if rgb is None:
            status = "N/A"
            rgb_text = "-"
        else:
            status = "MATCH" if sp.simplify(canonical - rgb) == 0 else "MISMATCH"
            rgb_text = sp.sstr(rgb)
        comparisons[name] = status
        print(
            f"{name:<18} {sp.sstr(canonical):>12} "
            f"{rows[name]['residual_nonzero_components']:>10} "
            f"{rgb_text:>12} {status:>10}"
        )

    result = {
        "dS1": args.dS1,
        "dS2": args.dS2,
        "canonical_recouplings": rows,
        "rgbeta_coefficients": {
            k: None if v is None else sp.sstr(v)
            for k, v in observed.items()
        },
        "comparisons": comparisons,
        "interpretation": (
            "Canonical coefficients are derived from the unique SU(2) invariant "
            "and standard generator normalization Tr(T^A T^B)=T(R) delta^AB. "
            "Any mismatch with RGBeta is a convention/implementation diagnostic."
        ),
    }

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2), encoding="utf-8")
        print()
        print(f"JSON summary: {args.output}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
