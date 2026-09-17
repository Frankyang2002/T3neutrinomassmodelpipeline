from __future__ import annotations

"""Canonical T3-E lambda12Cross recoupling.

This reconstructs the exact T3-E Cross operator independently of RGBeta's
standalone del[S2,...] reduction, which returns zero in RGBeta v1.2.0.

Basis
-----
The triplets are represented in the orthonormal spherical spin-1 basis
m=(+1,0,-1).  The invariant metric for two same-orientation triplet indices is

    C_{m m'} = (-1)^(1-m) delta_{m,-m'}.

For the RGBeta Cross field ordering {Bar[S1], S1, Bar[S2], S2},

    delta(BarS1,BarS2) delta(S1,S2)
  + delta(BarS1,S2)    delta(S1,BarS2)

becomes

    (BarS1^T C BarS2)(S1^T C S2)
  + (BarS1^T S2)(S1^T BarS2).

The T3-E mixing tensor is the same gauge-covariant embedding selected by the
Ward-identity test in ValidateRGBetaEInvariantRecoupling.py: charge-conjugation
map C on the Bar[S2] triplet leg.

The one-loop scalar cross contraction uses the standard fully symmetric real
quartic tensor convention.
"""

import argparse
import json
from collections import Counter
from itertools import combinations_with_replacement
from math import factorial
from pathlib import Path

import sympy as sp


EPS2 = sp.Matrix([[0, 1], [-1, 0]])
C1 = sp.Matrix([
    [0, 0, 1],
    [0, -1, 0],
    [1, 0, 0],
])

D = {}
for a in range(3):
    for i in range(2):
        for j in range(2):
            D[a, i, j] = sp.S.Zero
D[0, 0, 0] = 1
D[1, 0, 1] = 1 / sp.sqrt(2)
D[1, 1, 0] = 1 / sp.sqrt(2)
D[2, 1, 1] = 1


def refined_mix_pair(h1, h2, p, q, r, s):
    terms = (
        EPS2[h1, s] * EPS2[h2, q] * EPS2[p, r],
        EPS2[h1, q] * EPS2[h2, s] * EPS2[p, r],
        EPS2[h1, r] * EPS2[h2, q] * EPS2[p, s],
        EPS2[h1, q] * EPS2[h2, r] * EPS2[p, s],
        EPS2[h1, s] * EPS2[h2, p] * EPS2[q, r],
        EPS2[h1, p] * EPS2[h2, s] * EPS2[q, r],
        EPS2[h1, r] * EPS2[h2, p] * EPS2[q, s],
        EPS2[h1, p] * EPS2[h2, r] * EPS2[q, s],
    )
    return sp.simplify(sum(terms) / 8)


def mix_triplet_raw(h1, h2, a, b):
    return sp.simplify(sum(
        D[a, p, q] * D[b, r, s] * refined_mix_pair(h1, h2, p, q, r, s)
        for p in range(2)
        for q in range(2)
        for r in range(2)
        for s in range(2)
    ))


def mix_gauge_covariant(h1, h2, a, b):
    # Ward test selected C1 on the Bar[S2] triplet leg.
    return sp.simplify(sum(
        C1[v, b] * mix_triplet_raw(h1, h2, a, v)
        for v in range(3)
    ))


def complex_block(d, first):
    z, zb, variables = [], [], []
    for c in range(d):
        x = sp.Symbol(f"x{first + 2*c}")
        y = sp.Symbol(f"x{first + 2*c + 1}")
        z.append((x + sp.I*y) / sp.sqrt(2))
        zb.append((x - sp.I*y) / sp.sqrt(2))
        variables.extend((x, y))
    return z, zb, variables


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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    H, Hb, hv = complex_block(2, 1)
    S1, S1b, s1v = complex_block(3, 5)
    S2, S2b, s2v = complex_block(3, 11)
    variables = hv + s1v + s2v

    mix_poly = sp.expand(sum(
        mix_gauge_covariant(h1, h2, a, b)
        * H[h1] * H[h2] * S1[a] * S2b[b]
        for h1 in range(2)
        for h2 in range(2)
        for a in range(3)
        for b in range(3)
    ))

    # Exact Cross operator in spherical basis.
    same_orientation = (
        sum(S1b[a] * C1[a, b] * S2b[b] for a in range(3) for b in range(3))
        * sum(S1[a] * C1[a, b] * S2[b] for a in range(3) for b in range(3))
    )
    opposite_orientation = (
        sum(S1b[a] * S2[a] for a in range(3))
        * sum(S1[a] * S2b[a] for a in range(3))
    )
    cross_poly = sp.expand(same_orientation + opposite_orientation)

    target = tensor_from_polynomial(mix_poly, variables)
    cross = tensor_from_polynomial(cross_poly, variables)
    generated = scalar_cross_beta(target, cross, len(variables))
    coefficient, residual = project(generated, target)

    print("T3-E exact Cross recoupling")
    print(f"projection coefficient : {coefficient}")
    print(f"nonzero residual       : {len(residual)}")
    print(f"closed                 : {len(residual) == 0}")

    payload = {
        "status": "Success" if len(residual) == 0 else "ResidualNonzero",
        "coefficient": sp.sstr(coefficient),
        "residual_nonzero_components": len(residual),
        "closed": len(residual) == 0,
        "basis": "orthonormal spherical spin-1",
        "cross_operator": (
            "(BarS1^T C BarS2)(S1^T C S2) "
            "+ (BarS1^T S2)(S1^T BarS2)"
        ),
    }

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        print(f"JSON summary: {args.output}")

    return 0 if len(residual) == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
