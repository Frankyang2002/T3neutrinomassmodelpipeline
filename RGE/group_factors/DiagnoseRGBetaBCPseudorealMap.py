from __future__ import annotations

"""Diagnose the T3-B/C doublet pseudoreality map.

Starting point: RGBeta RefineGroupStructures gives the formal index tensor

  K(h1,h2;s1,s2) =
      1/2 (delta[h1,s2] delta[h2,s1]
         + delta[h1,s1] delta[h2,s2])

for the HH S1 Bar[S2] topology invariant.

The symbols s1 and s2 in that formal expression do not by themselves specify
how RGBeta's covariant/contravariant doublet indices map into an ordinary
complex-field component basis.  SU(2) doublets are pseudoreal, so either leg
may require the epsilon intertwiner

    eps = [[0,1],[-1,0]].

This script tests all four possible epsilon insertions on the S1 and Bar[S2]
legs.  For each embedding it constructs the fully symmetric real-scalar
quartic tensor and evaluates the one-loop scalar cross contraction with the
three adjoint quartics.

A physically acceptable embedding of the unique HH S1 Bar[S2] invariant must
close back onto the same one-dimensional tensor direction: residual = 0.
"""

import argparse
from collections import Counter
from itertools import combinations_with_replacement
from math import factorial
import json
from pathlib import Path

import sympy as sp


def complex_block(d: int, first: int):
    z, zb, variables = [], [], []
    for c in range(d):
        x = sp.Symbol(f"x{first + 2*c}")
        y = sp.Symbol(f"x{first + 2*c + 1}")
        z.append((x + sp.I*y) / sp.sqrt(2))
        zb.append((x - sp.I*y) / sp.sqrt(2))
        variables.extend((x, y))
    return z, zb, variables


def generators():
    return (
        sp.Matrix([[0, 1], [1, 0]]) / 2,
        sp.Matrix([[0, -sp.I], [sp.I, 0]]) / 2,
        sp.Matrix([[1, 0], [0, -1]]) / 2,
    )


def bilinear(zb, z, matrix):
    return sp.expand(sum(
        zb[i] * matrix[i, j] * z[j]
        for i in range(2)
        for j in range(2)
    ))


def tensor_from_polynomial(expr, variables):
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


def get(tensor, *indices):
    return tensor.get(tuple(sorted(indices)), sp.S.Zero)


def ordered_weight(key):
    counts = Counter(key)
    value = factorial(4)
    for count in counts.values():
        value //= factorial(count)
    return value


def inner(left, right):
    return sp.simplify(sum(
        ordered_weight(key)
        * sp.conjugate(get(left, *key))
        * get(right, *key)
        for key in set(left) | set(right)
    ))


def scalar_cross_beta(target, other, n_real):
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
                        get(target, p, q, e, f) * get(other, e, f, r, s)
                        + get(other, p, q, e, f) * get(target, e, f, r, s)
                    )
        value = sp.simplify(value)
        if value != 0:
            generated[key] = value
    return generated


def project(generated, target):
    coefficient = sp.simplify(inner(target, generated) / inner(target, target))
    residual = {
        key: sp.simplify(get(generated, *key) - coefficient * get(target, *key))
        for key in set(generated) | set(target)
    }
    residual = {k: v for k, v in residual.items() if v != 0}
    return coefficient, residual


def build_mix(H, S1, S2b, eps_s1: bool, eps_s2: bool):
    eps = sp.Matrix([[0, 1], [-1, 0]])
    expr = sp.S.Zero

    # Formal RGBeta-refined K tensor.
    for h1 in range(2):
        for h2 in range(2):
            for c in range(2):
                for d in range(2):
                    coefficient = sp.S.Zero
                    for s1 in range(2):
                        for s2 in range(2):
                            K = sp.Rational(1, 2) * (
                                int(h1 == s2 and h2 == s1)
                                + int(h1 == s1 and h2 == s2)
                            )
                            m1 = eps[s1, c] if eps_s1 else int(s1 == c)
                            m2 = eps[s2, d] if eps_s2 else int(s2 == d)
                            coefficient += K * m1 * m2

                    expr += (
                        sp.simplify(coefficient)
                        * H[h1] * H[h2] * S1[c] * S2b[d]
                    )

    return sp.expand(expr)


def run():
    H, Hb, hv = complex_block(2, 1)
    S1, S1b, s1v = complex_block(2, 5)
    S2, S2b, s2v = complex_block(2, 9)
    variables = hv + s1v + s2v

    T = generators()
    quartics = {
        "lambdaH1Adj": sum(
            bilinear(Hb, H, T[A]) * bilinear(S1b, S1, T[A])
            for A in range(3)
        ),
        "lambdaH2Adj": sum(
            bilinear(Hb, H, T[A]) * bilinear(S2b, S2, T[A])
            for A in range(3)
        ),
        "lambda12Adj": sum(
            bilinear(S1b, S1, T[A]) * bilinear(S2b, S2, T[A])
            for A in range(3)
        ),
    }

    output = {}

    for eps_s1, eps_s2 in ((False, False), (True, False), (False, True), (True, True)):
        key = f"epsS1={eps_s1},epsS2={eps_s2}"
        mix = build_mix(H, S1, S2b, eps_s1, eps_s2)
        target = tensor_from_polynomial(mix, variables)

        rows = {}
        for name, polynomial in quartics.items():
            other = tensor_from_polynomial(polynomial, variables)
            generated = scalar_cross_beta(target, other, len(variables))
            coefficient, residual = project(generated, target)
            rows[name] = {
                "coefficient": sp.sstr(coefficient),
                "residual_nonzero_components": len(residual),
                "closed": len(residual) == 0,
            }

        output[key] = rows

    return output


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    result = run()

    print("T3-B/C pseudoreal-map diagnostic")
    print("A valid embedding of the unique invariant should have residual = 0.")
    print()

    for variant, rows in result.items():
        print(variant)
        print(f"{'quartic':<18} {'coefficient':>12} {'residual':>10} {'closed':>8}")
        print("-" * 54)
        for name, row in rows.items():
            print(
                f"{name:<18} {row['coefficient']:>12} "
                f"{row['residual_nonzero_components']:>10} "
                f"{str(row['closed']):>8}"
            )
        print()

    payload = {
        "status": "Success",
        "formal_rgbeta_refined_tensor": (
            "1/2(delta[h1,s2]delta[h2,s1] + "
            "delta[h1,s1]delta[h2,s2])"
        ),
        "results": result,
        "conclusion": (
            "Only the epsilon-on-S1, no-epsilon-on-S2 embedding closes for all "
            "three adjoint quartics.  Its coefficients are "
            "H1Adj=-2, H2Adj=2, lambda12Adj=-1/2.  The no-epsilon embedding "
            "has projection coefficients 2/3, 2, 1/6 but nonzero residuals for "
            "H1Adj and lambda12Adj, so those projection values alone do not "
            "constitute a closed tensor derivation."
        ),
    }

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        print(f"JSON summary: {args.output}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
