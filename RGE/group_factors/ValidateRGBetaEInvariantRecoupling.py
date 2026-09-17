from __future__ import annotations

"""Validate the T3-E RGBeta mixing invariant in an ordinary component basis.

Starting point
--------------
RGBeta RefineGroupStructures gives the T3-E mixing invariant as an 8-term
expression built from fundamental SU(2) epsilon tensors.  The external S1 and
Bar[S2] triplet indices are represented as symmetric pairs of fundamental
indices.

This script:

1. Maps the symmetric-pair representation to an orthonormal j=1 basis
   |1,+1>, |1,0>, |1,-1>.
2. Tests all four possible spin-1 charge-conjugation maps on S1 and Bar[S2].
3. Applies the physical Ward identity for fields {H,H,S1,Bar[S2]}:
       +T_H + T_H + T_S1 - T_S2^T.
4. Selects the gauge-covariant embedding.
5. Converts that invariant to a fully symmetric real-scalar quartic tensor.
6. Evaluates the standard one-loop scalar cross contraction with
   lambdaH1Adj, lambdaH2Adj and lambda12Adj.
7. Projects back onto the mixing tensor and reports the exact residual.

No A--E numerical beta coefficients are fitted.

Important:
RGBeta's T3-E Cross invariant is not validated here.  The direct Wolfram
probe showed RefineGroupStructures[Cross] == 0 in RGBeta v1.2.0 because the
standalone del[S2,...] structure is routed through twoIndexRepDelta, whose
current definition returns 0.  Therefore Cross needs a separate treatment.
"""

import argparse
import json
from collections import Counter
from itertools import combinations_with_replacement
from math import factorial
from pathlib import Path

import sympy as sp


EPS2 = sp.Matrix([[0, 1], [-1, 0]])

# Orthonormal symmetric-pair map for j=1, ordered as m=+1,0,-1.
D = {}
for a in range(3):
    for i in range(2):
        for j in range(2):
            D[a, i, j] = sp.S.Zero
D[0, 0, 0] = 1
D[1, 0, 1] = 1 / sp.sqrt(2)
D[1, 1, 0] = 1 / sp.sqrt(2)
D[2, 1, 1] = 1

# Spin-1 charge-conjugation / real-structure intertwiner
# C_{m,m'} = (-1)^(1-m) delta_{m,-m'} in the (+1,0,-1) basis.
C1 = sp.Matrix([
    [0, 0, 1],
    [0, -1, 0],
    [1, 0, 0],
])


def fundamental_generators():
    return (
        sp.Matrix([[0, 1], [1, 0]]) / 2,
        sp.Matrix([[0, -sp.I], [sp.I, 0]]) / 2,
        sp.Matrix([[1, 0], [0, -1]]) / 2,
    )


def triplet_generators():
    jp = sp.Matrix([
        [0, sp.sqrt(2), 0],
        [0, 0, sp.sqrt(2)],
        [0, 0, 0],
    ])
    jm = jp.T
    return (
        (jp + jm) / 2,
        (jp - jm) / (2 * sp.I),
        sp.diag(1, 0, -1),
    )


def refined_mix_pair(h1, h2, p, q, r, s):
    """Exact 8-term RGBeta-refined T3-E tensor in symmetric-pair indices."""
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


def mix_with_maps(h1, h2, a, b, map_s1=False, map_s2bar=False):
    return sp.simplify(sum(
        (C1[u, a] if map_s1 else int(u == a))
        * (C1[v, b] if map_s2bar else int(v == b))
        * mix_triplet_raw(h1, h2, u, v)
        for u in range(3)
        for v in range(3)
    ))


def ward_residual_count(map_s1=False, map_s2bar=False):
    tf = fundamental_generators()
    tt = triplet_generators()
    total = 0
    by_generator = []

    for A in range(3):
        nonzero = 0
        for h1 in range(2):
            for h2 in range(2):
                for a in range(3):
                    for b in range(3):
                        value = sp.S.Zero

                        for p in range(2):
                            value += (
                                tf[A][h1, p]
                                * mix_with_maps(p, h2, a, b, map_s1, map_s2bar)
                            )
                            value += (
                                tf[A][h2, p]
                                * mix_with_maps(h1, p, a, b, map_s1, map_s2bar)
                            )

                        for c in range(3):
                            value += (
                                tt[A][a, c]
                                * mix_with_maps(h1, h2, c, b, map_s1, map_s2bar)
                            )
                            # Bar[S2] transforms in the conjugate representation.
                            value -= (
                                tt[A][c, b]
                                * mix_with_maps(h1, h2, a, c, map_s1, map_s2bar)
                            )

                        if sp.simplify(value) != 0:
                            nonzero += 1

        by_generator.append(nonzero)
        total += nonzero

    return by_generator, total


def complex_block(d, first):
    z, zb, variables = [], [], []
    for c in range(d):
        x = sp.Symbol(f"x{first + 2*c}")
        y = sp.Symbol(f"x{first + 2*c + 1}")
        z.append((x + sp.I * y) / sp.sqrt(2))
        zb.append((x - sp.I * y) / sp.sqrt(2))
        variables.extend((x, y))
    return z, zb, variables


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
    residual = {
        key: sp.simplify(
            get(generated, *key) - coefficient * get(target, *key)
        )
        for key in set(generated) | set(target)
    }
    residual = {key: value for key, value in residual.items() if value != 0}
    return coefficient, len(residual)


def recoupling_for_embedding(map_s1, map_s2bar):
    H, Hb, hv = complex_block(2, 1)
    S1, S1b, s1v = complex_block(3, 5)
    S2, S2b, s2v = complex_block(3, 11)
    variables = hv + s1v + s2v

    mix_poly = sp.expand(sum(
        mix_with_maps(h1, h2, a, b, map_s1, map_s2bar)
        * H[h1] * H[h2] * S1[a] * S2b[b]
        for h1 in range(2)
        for h2 in range(2)
        for a in range(3)
        for b in range(3)
    ))
    target = tensor_from_polynomial(mix_poly, variables)

    tf = fundamental_generators()
    tt = triplet_generators()

    quartics = {
        "lambdaH1Adj": sum(
            bilinear(Hb, H, tf[A]) * bilinear(S1b, S1, tt[A])
            for A in range(3)
        ),
        "lambdaH2Adj": sum(
            bilinear(Hb, H, tf[A]) * bilinear(S2b, S2, tt[A])
            for A in range(3)
        ),
        "lambda12Adj": sum(
            bilinear(S1b, S1, tt[A]) * bilinear(S2b, S2, tt[A])
            for A in range(3)
        ),
    }

    rows = {}
    for name, polynomial in quartics.items():
        other = tensor_from_polynomial(polynomial, variables)
        generated = scalar_cross_beta(target, other, len(variables))
        coefficient, residual = project(generated, target)
        rows[name] = {
            "coefficient": sp.sstr(coefficient),
            "residual_nonzero_components": residual,
            "closed": residual == 0,
        }

    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    print("T3-E RGBeta invariant: Ward identity + scalar recoupling")
    print()

    ward = {}
    selected = None

    for map_s1, map_s2bar in (
        (False, False),
        (True, False),
        (False, True),
        (True, True),
    ):
        label = f"C1_on_S1={map_s1},C1_on_BarS2={map_s2bar}"
        by_gen, total = ward_residual_count(map_s1, map_s2bar)
        ward[label] = {
            "by_generator": by_gen,
            "total": total,
            "gauge_invariant": total == 0,
        }
        print(
            f"{label:<39} "
            f"A1={by_gen[0]:>2} A2={by_gen[1]:>2} A3={by_gen[2]:>2} "
            f"total={total:>2} invariant={total == 0}"
        )
        if total == 0:
            selected = (map_s1, map_s2bar, label)

    if selected is None:
        print("\nNo gauge-covariant embedding found.")
        return 1

    map_s1, map_s2bar, selected_label = selected
    print(f"\nSelected gauge-covariant embedding: {selected_label}")
    print()

    recoupling = recoupling_for_embedding(map_s1, map_s2bar)

    print(f"{'quartic':<18} {'coefficient':>12} {'residual':>10} {'closed':>8}")
    print("-" * 54)
    for name, row in recoupling.items():
        print(
            f"{name:<18} {row['coefficient']:>12} "
            f"{row['residual_nonzero_components']:>10} "
            f"{str(row['closed']):>8}"
        )

    payload = {
        "status": "Success",
        "ward_identity": ward,
        "selected_embedding": selected_label,
        "recoupling": recoupling,
        "cross_status": (
            "Unresolved here: direct RGBeta v1.2.0 probe gave "
            "RefineGroupStructures[lambda12Cross invariant] == 0."
        ),
    }

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        print(f"\nJSON summary: {args.output}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
