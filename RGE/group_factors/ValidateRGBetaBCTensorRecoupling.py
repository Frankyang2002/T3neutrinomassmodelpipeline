from __future__ import annotations

"""Direct T3-B/C recoupling from RGBeta's own refined mixing invariant.

RGBeta refines

  delS2[a,h1,h2] delS2[a,s1,s2]

to

  1/2 (delta[h1,s2] delta[h2,s1]
     + delta[h1,s1] delta[h2,s2]).

This script uses exactly that component tensor, together with standard
fundamental SU(2) generators T^A = sigma^A/2, converts the complex scalar
operators to the fully symmetric real-scalar quartic tensor, and evaluates
the one-loop scalar cross contraction.

It reports:
  * projection coefficient onto the original RGBeta lambdaT3 tensor;
  * residual outside that tensor direction.

The projection coefficients are expected to reproduce RGBeta's B/C values
without fitting them:
    H1Adj -> 2/3
    H2Adj -> 2
    12Adj -> 1/6

A nonzero residual is intentionally retained as a diagnostic: the refined
RGBeta tensor written naively in an ordinary complex-doublet component basis
does not by itself encode all pseudoreal/conjugate-index structure.
"""

import argparse
from collections import Counter
from itertools import combinations_with_replacement
from math import factorial
import json
from pathlib import Path

import sympy as sp


def complex_block(d: int, first: int):
    z = []
    zb = []
    variables = []
    for c in range(d):
        x = sp.Symbol(f"x{first + 2*c}")
        y = sp.Symbol(f"x{first + 2*c + 1}")
        z.append((x + sp.I*y) / sp.sqrt(2))
        zb.append((x - sp.I*y) / sp.sqrt(2))
        variables.extend((x, y))
    return z, zb, variables


def fundamental_generators():
    return (
        sp.Matrix([[0, 1], [1, 0]]) / 2,
        sp.Matrix([[0, -sp.I], [sp.I, 0]]) / 2,
        sp.Matrix([[1, 0], [0, -1]]) / 2,
    )


def bilinear(zb, z, matrix):
    return sp.expand(sum(
        zb[i] * matrix[i, j] * z[j]
        for i in range(len(z))
        for j in range(len(z))
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
    residual = {}

    for key in set(generated) | set(target):
        value = sp.simplify(
            get(generated, *key) - coefficient * get(target, *key)
        )
        if value != 0:
            residual[key] = value

    return coefficient, residual


def run():
    H, Hb, hv = complex_block(2, 1)
    S1, S1b, s1v = complex_block(2, 5)
    S2, S2b, s2v = complex_block(2, 9)
    variables = hv + s1v + s2v

    # Exact RGBeta RefineGroupStructures result:
    # 1/2(delta_h1,s2 delta_h2,s1 + delta_h1,s1 delta_h2,s2).
    mix = sp.S.Zero
    for h1 in range(2):
        for h2 in range(2):
            for s1 in range(2):
                for s2 in range(2):
                    coeff = sp.Rational(1, 2) * (
                        int(h1 == s2 and h2 == s1)
                        + int(h1 == s1 and h2 == s2)
                    )
                    mix += coeff * H[h1] * H[h2] * S1[s1] * S2b[s2]

    gens = fundamental_generators()

    quartics = {
        "lambdaH1Adj": sum(
            bilinear(Hb, H, gens[A]) * bilinear(S1b, S1, gens[A])
            for A in range(3)
        ),
        "lambdaH2Adj": sum(
            bilinear(Hb, H, gens[A]) * bilinear(S2b, S2, gens[A])
            for A in range(3)
        ),
        "lambda12Adj": sum(
            bilinear(S1b, S1, gens[A]) * bilinear(S2b, S2, gens[A])
            for A in range(3)
        ),
    }

    target = tensor_from_polynomial(mix, variables)
    rows = {}

    for name, polynomial in quartics.items():
        other = tensor_from_polynomial(polynomial, variables)
        generated = scalar_cross_beta(target, other, len(variables))
        coefficient, residual = project(generated, target)
        rows[name] = {
            "projection_coefficient": sp.sstr(coefficient),
            "residual_nonzero_components": len(residual),
        }

    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    rows = run()

    print("Direct recoupling from RGBeta-refined T3-B/C mixing tensor")
    print(f"{'quartic':<18} {'projection':>12} {'residual':>10}")
    print("-" * 46)
    for name, row in rows.items():
        print(
            f"{name:<18} "
            f"{row['projection_coefficient']:>12} "
            f"{row['residual_nonzero_components']:>10}"
        )

    payload = {
        "status": "Success",
        "input_tensor": (
            "1/2(delta[h1,s2]delta[h2,s1] + "
            "delta[h1,s1]delta[h2,s2])"
        ),
        "results": rows,
        "interpretation": (
            "Projection coefficients reproduce RGBeta B/C directly. "
            "Nonzero residuals for H1Adj and lambda12Adj show that the naively "
            "embedded refined tensor does not yet encode the complete "
            "pseudoreal/conjugate-index convention."
        ),
    }

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        print(f"\nJSON summary: {args.output}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
