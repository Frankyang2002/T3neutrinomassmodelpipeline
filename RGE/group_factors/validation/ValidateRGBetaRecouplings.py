from __future__ import annotations

"""Independent RGBeta recoupling diagnostics and validations.

Each former calculation remains a separate function set; this module only consolidates their source-file packaging. Use the first positional argument to select the former script."""

import sys
import argparse
from collections import Counter
from itertools import combinations_with_replacement
from math import factorial
import json
from pathlib import Path
import sympy as sp

# ---------------------------------------------------------------------------
# T3-B/C pseudoreal-map diagnostic
# Former source: DiagnoseRGBetaBCPseudorealMap.py
# ---------------------------------------------------------------------------

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




def bc_diagonal_complex_block(d: int, first: int):
    z, zb, variables = [], [], []
    for c in range(d):
        x = sp.Symbol(f"x{first + 2*c}")
        y = sp.Symbol(f"x{first + 2*c + 1}")
        z.append((x + sp.I*y) / sp.sqrt(2))
        zb.append((x - sp.I*y) / sp.sqrt(2))
        variables.extend((x, y))
    return z, zb, variables


def bc_diagonal_generators():
    return (
        sp.Matrix([[0, 1], [1, 0]]) / 2,
        sp.Matrix([[0, -sp.I], [sp.I, 0]]) / 2,
        sp.Matrix([[1, 0], [0, -1]]) / 2,
    )


def bc_diagonal_bilinear(zb, z, matrix):
    return sp.expand(sum(
        zb[i] * matrix[i, j] * z[j]
        for i in range(2)
        for j in range(2)
    ))


def bc_diagonal_tensor_from_polynomial(expr, variables):
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


def bc_diagonal_get(tensor, *indices):
    return tensor.get(tuple(sorted(indices)), sp.S.Zero)


def bc_diagonal_ordered_weight(key):
    counts = Counter(key)
    value = factorial(4)
    for count in counts.values():
        value //= factorial(count)
    return value


def bc_diagonal_inner(left, right):
    return sp.simplify(sum(
        bc_diagonal_ordered_weight(key)
        * sp.conjugate(bc_diagonal_get(left, *key))
        * bc_diagonal_get(right, *key)
        for key in set(left) | set(right)
    ))


def bc_diagonal_scalar_cross_beta(target, other, n_real):
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
                        bc_diagonal_get(target, p, q, e, f) * bc_diagonal_get(other, e, f, r, s)
                        + bc_diagonal_get(other, p, q, e, f) * bc_diagonal_get(target, e, f, r, s)
                    )
        value = sp.simplify(value)
        if value != 0:
            generated[key] = value
    return generated


def bc_diagonal_project(generated, target):
    coefficient = sp.simplify(bc_diagonal_inner(target, generated) / bc_diagonal_inner(target, target))
    residual = {
        key: sp.simplify(bc_diagonal_get(generated, *key) - coefficient * bc_diagonal_get(target, *key))
        for key in set(generated) | set(target)
    }
    residual = {k: v for k, v in residual.items() if v != 0}
    return coefficient, residual


def build_bc_diagonal_mix(H, S1, S2b, eps_s1: bool, eps_s2: bool):
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


def bc_diagonal_run():
    H, Hb, hv = bc_diagonal_complex_block(2, 1)
    S1, S1b, s1v = bc_diagonal_complex_block(2, 5)
    S2, S2b, s2v = bc_diagonal_complex_block(2, 9)
    variables = hv + s1v + s2v

    T = bc_diagonal_generators()
    quartics = {
        "lambdaH1Adj": sum(
            bc_diagonal_bilinear(Hb, H, T[A]) * bc_diagonal_bilinear(S1b, S1, T[A])
            for A in range(3)
        ),
        "lambdaH2Adj": sum(
            bc_diagonal_bilinear(Hb, H, T[A]) * bc_diagonal_bilinear(S2b, S2, T[A])
            for A in range(3)
        ),
        "lambda12Adj": sum(
            bc_diagonal_bilinear(S1b, S1, T[A]) * bc_diagonal_bilinear(S2b, S2, T[A])
            for A in range(3)
        ),
    }

    output = {}

    for eps_s1, eps_s2 in ((False, False), (True, False), (False, True), (True, True)):
        key = f"epsS1={eps_s1},epsS2={eps_s2}"
        mix = build_bc_diagonal_mix(H, S1, S2b, eps_s1, eps_s2)
        target = bc_diagonal_tensor_from_polynomial(mix, variables)

        rows = {}
        for name, polynomial in quartics.items():
            other = bc_diagonal_tensor_from_polynomial(polynomial, variables)
            generated = bc_diagonal_scalar_cross_beta(target, other, len(variables))
            coefficient, residual = bc_diagonal_project(generated, target)
            rows[name] = {
                "coefficient": sp.sstr(coefficient),
                "residual_nonzero_components": len(residual),
                "closed": len(residual) == 0,
            }

        output[key] = rows

    return output


def run_bc_diagonal_diagnostic_cli():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    result = bc_diagonal_run()

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
            "has projection coefficients 2/3, 2, 1/6 but nonzero bc_ward_residuals for "
            "H1Adj and lambda12Adj, so those projection values alone do not "
            "constitute a closed tensor derivation."
        ),
    }

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        print(f"JSON summary: {args.output}")

    return 0


# ---------------------------------------------------------------------------
# T3-B/C mixing-invariant Ward-identity validation
# Former source: ValidateRGBetaBCMixInvariantWardIdentity.py
# ---------------------------------------------------------------------------

"""Ward-identity test for the T3-B/C RGBeta mixing invariant.

RGBeta defines, for dS1=dS2=2,

    I(h1,h2,s1,s2bar)
      = delS2[a,h1,h2] delS2[a,s1,s2bar]

and RefineGroupStructures gives the formal component tensor

    K = 1/2 (delta[h1,s2bar] delta[h2,s1]
           + delta[h1,s1]    delta[h2,s2bar]).

The physical field content is {H,H,S1,Bar[S2]}, so the infinitesimal SU(2)
action is

    +T on H,
    +T on H,
    +T on S1,
    -T^T on Bar[S2].

Because SU(2) doublets are pseudoreal, the formal RGBeta index on the S1 leg
must be mapped into the ordinary complex-component convention with epsilon.
This is exactly the same leg identified independently by the scalar-loop
closure diagnostic.

This script compares:
  1. the raw RGBeta-refined formal delta tensor;
  2. epsilon inserted on S1;
  3. epsilon inserted on Bar[S2];
  4. epsilon inserted on both scalar legs.

A correct ordinary-component embedding should satisfy the Ward identity
component-by-component.
"""



def bc_ward_generators():
    return (
        sp.Matrix([[0, 1], [1, 0]]) / 2,
        sp.Matrix([[0, -sp.I], [sp.I, 0]]) / 2,
        sp.Matrix([[1, 0], [0, -1]]) / 2,
    )


EPS = sp.Matrix([[0, 1], [-1, 0]])


def build_bc_ward_raw_tensor(h1, h2, s1, sb2):
    return sp.Rational(1, 2) * (
        int(h1 == sb2 and h2 == s1)
        + int(h1 == s1 and h2 == sb2)
    )


def eps_on_s1(h1, h2, s1, sb2):
    return sp.simplify(sum(
        build_bc_ward_raw_tensor(h1, h2, u, sb2) * EPS[u, s1]
        for u in range(2)
    ))


def eps_on_s2bar(h1, h2, s1, sb2):
    return sp.simplify(sum(
        build_bc_ward_raw_tensor(h1, h2, s1, u) * EPS[u, sb2]
        for u in range(2)
    ))


def eps_on_both(h1, h2, s1, sb2):
    return sp.simplify(sum(
        build_bc_ward_raw_tensor(h1, h2, u, v) * EPS[u, s1] * EPS[v, sb2]
        for u in range(2)
        for v in range(2)
    ))


def bc_ward_component(tensor, T, h1, h2, s1, sb2):
    out = sp.S.Zero

    for p in range(2):
        out += T[h1, p] * tensor(p, h2, s1, sb2)
        out += T[h2, p] * tensor(h1, p, s1, sb2)
        out += T[s1, p] * tensor(h1, h2, p, sb2)

        # Bar[S2] transforms in the anti-fundamental:
        out -= T[p, sb2] * tensor(h1, h2, s1, p)

    return sp.simplify(out)


def bc_ward_residuals(tensor):
    rows = []
    for A, T in enumerate(bc_ward_generators(), start=1):
        nonzero = {}
        for h1 in range(2):
            for h2 in range(2):
                for s1 in range(2):
                    for sb2 in range(2):
                        r = bc_ward_component(tensor, T, h1, h2, s1, sb2)
                        if r != 0:
                            nonzero[(h1, h2, s1, sb2)] = r
        rows.append((A, nonzero))
    return rows


def print_bc_ward_result(name, rows):
    print(name)
    total = 0
    for A, nonzero in rows:
        total += len(nonzero)
        print(f"  generator A={A}: nonzero Ward components = {len(nonzero)}")
        for key, value in list(nonzero.items())[:6]:
            print(f"    {key}: {value}")
        if len(nonzero) > 6:
            print("    ...")
    print(f"  total nonzero Ward components = {total}")
    print(f"  gauge invariant = {total == 0}")
    print()


def run_bc_ward_validation_cli():
    print("T3-B/C RGBeta mixing-invariant Ward-identity test")
    print("field ordering: H, H, S1, Bar[S2]")
    print()

    cases = [
        ("Raw RGBeta-refined tensor", build_bc_ward_raw_tensor),
        ("Epsilon on S1", eps_on_s1),
        ("Epsilon on Bar[S2]", eps_on_s2bar),
        ("Epsilon on both scalar legs", eps_on_both),
    ]

    totals = {}
    for name, tensor in cases:
        rows = bc_ward_residuals(tensor)
        print_bc_ward_result(name, rows)
        totals[name] = sum(len(x) for _, x in rows)

    ok = (
        totals["Raw RGBeta-refined tensor"] != 0
        and totals["Epsilon on S1"] == 0
        and totals["Epsilon on Bar[S2]"] != 0
    )

    if ok:
        print(
            "CONCLUSION: in an ordinary complex-component basis the RGBeta "
            "formal tensor becomes gauge invariant only when the SU(2) "
            "pseudoreality intertwiner is placed on the S1 leg. This matches "
            "the independent scalar-loop closure diagnostic."
        )
        return 0

    print("CONCLUSION: expected pseudoreal-map pattern was not obtained.")
    return 1


# ---------------------------------------------------------------------------
# T3-B/C tensor recoupling validation
# Former source: ValidateRGBetaBCTensorRecoupling.py
# ---------------------------------------------------------------------------

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




def bc_tensor_complex_block(d: int, first: int):
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


def bc_tensor_fundamental_generators():
    return (
        sp.Matrix([[0, 1], [1, 0]]) / 2,
        sp.Matrix([[0, -sp.I], [sp.I, 0]]) / 2,
        sp.Matrix([[1, 0], [0, -1]]) / 2,
    )


def bc_tensor_bilinear(zb, z, matrix):
    return sp.expand(sum(
        zb[i] * matrix[i, j] * z[j]
        for i in range(len(z))
        for j in range(len(z))
    ))


def bc_tensor_tensor_from_polynomial(expr, variables):
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


def bc_tensor_get(tensor, *indices):
    return tensor.get(tuple(sorted(indices)), sp.S.Zero)


def bc_tensor_ordered_weight(key):
    counts = Counter(key)
    value = factorial(4)
    for count in counts.values():
        value //= factorial(count)
    return value


def bc_tensor_inner(left, right):
    return sp.simplify(sum(
        bc_tensor_ordered_weight(key)
        * sp.conjugate(bc_tensor_get(left, *key))
        * bc_tensor_get(right, *key)
        for key in set(left) | set(right)
    ))


def bc_tensor_scalar_cross_beta(target, other, n_real):
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
                        bc_tensor_get(target, p, q, e, f) * bc_tensor_get(other, e, f, r, s)
                        + bc_tensor_get(other, p, q, e, f) * bc_tensor_get(target, e, f, r, s)
                    )

        value = sp.simplify(value)
        if value != 0:
            generated[key] = value

    return generated


def bc_tensor_project(generated, target):
    coefficient = sp.simplify(bc_tensor_inner(target, generated) / bc_tensor_inner(target, target))
    residual = {}

    for key in set(generated) | set(target):
        value = sp.simplify(
            bc_tensor_get(generated, *key) - coefficient * bc_tensor_get(target, *key)
        )
        if value != 0:
            residual[key] = value

    return coefficient, residual


def bc_tensor_run():
    H, Hb, hv = bc_tensor_complex_block(2, 1)
    S1, S1b, s1v = bc_tensor_complex_block(2, 5)
    S2, S2b, s2v = bc_tensor_complex_block(2, 9)
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

    gens = bc_tensor_fundamental_generators()

    quartics = {
        "lambdaH1Adj": sum(
            bc_tensor_bilinear(Hb, H, gens[A]) * bc_tensor_bilinear(S1b, S1, gens[A])
            for A in range(3)
        ),
        "lambdaH2Adj": sum(
            bc_tensor_bilinear(Hb, H, gens[A]) * bc_tensor_bilinear(S2b, S2, gens[A])
            for A in range(3)
        ),
        "lambda12Adj": sum(
            bc_tensor_bilinear(S1b, S1, gens[A]) * bc_tensor_bilinear(S2b, S2, gens[A])
            for A in range(3)
        ),
    }

    target = bc_tensor_tensor_from_polynomial(mix, variables)
    rows = {}

    for name, polynomial in quartics.items():
        other = bc_tensor_tensor_from_polynomial(polynomial, variables)
        generated = bc_tensor_scalar_cross_beta(target, other, len(variables))
        coefficient, residual = bc_tensor_project(generated, target)
        rows[name] = {
            "projection_coefficient": sp.sstr(coefficient),
            "residual_nonzero_components": len(residual),
        }

    return rows


def run_bc_tensor_validation_cli():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    rows = bc_tensor_run()

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
            "Nonzero bc_ward_residuals for H1Adj and lambda12Adj show that the naively "
            "embedded refined tensor does not yet encode the complete "
            "pseudoreal/conjugate-index convention."
        ),
    }

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        print(f"\nJSON summary: {args.output}")

    return 0


# ---------------------------------------------------------------------------
# T3-E invariant recoupling validation
# Former source: ValidateRGBetaEInvariantRecoupling.py
# ---------------------------------------------------------------------------

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




e_invariant_EPS2 = sp.Matrix([[0, 1], [-1, 0]])

# Orthonormal symmetric-pair map for j=1, ordered as m=+1,0,-1.
e_invariant_D = {}
for a in range(3):
    for i in range(2):
        for j in range(2):
            e_invariant_D[a, i, j] = sp.S.Zero
e_invariant_D[0, 0, 0] = 1
e_invariant_D[1, 0, 1] = 1 / sp.sqrt(2)
e_invariant_D[1, 1, 0] = 1 / sp.sqrt(2)
e_invariant_D[2, 1, 1] = 1

# Spin-1 charge-conjugation / real-structure intertwiner
# C_{m,m'} = (-1)^(1-m) delta_{m,-m'} in the (+1,0,-1) basis.
e_invariant_C1 = sp.Matrix([
    [0, 0, 1],
    [0, -1, 0],
    [1, 0, 0],
])


def e_invariant_fundamental_generators():
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


def e_invariant_refined_mix_pair(h1, h2, p, q, r, s):
    """Exact 8-term RGBeta-refined T3-E tensor in symmetric-pair indices."""
    terms = (
        e_invariant_EPS2[h1, s] * e_invariant_EPS2[h2, q] * e_invariant_EPS2[p, r],
        e_invariant_EPS2[h1, q] * e_invariant_EPS2[h2, s] * e_invariant_EPS2[p, r],
        e_invariant_EPS2[h1, r] * e_invariant_EPS2[h2, q] * e_invariant_EPS2[p, s],
        e_invariant_EPS2[h1, q] * e_invariant_EPS2[h2, r] * e_invariant_EPS2[p, s],
        e_invariant_EPS2[h1, s] * e_invariant_EPS2[h2, p] * e_invariant_EPS2[q, r],
        e_invariant_EPS2[h1, p] * e_invariant_EPS2[h2, s] * e_invariant_EPS2[q, r],
        e_invariant_EPS2[h1, r] * e_invariant_EPS2[h2, p] * e_invariant_EPS2[q, s],
        e_invariant_EPS2[h1, p] * e_invariant_EPS2[h2, r] * e_invariant_EPS2[q, s],
    )
    return sp.simplify(sum(terms) / 8)


def e_invariant_mix_triplet_raw(h1, h2, a, b):
    return sp.simplify(sum(
        e_invariant_D[a, p, q] * e_invariant_D[b, r, s] * e_invariant_refined_mix_pair(h1, h2, p, q, r, s)
        for p in range(2)
        for q in range(2)
        for r in range(2)
        for s in range(2)
    ))


def build_e_invariant_mix_with_maps(h1, h2, a, b, map_s1=False, map_s2bar=False):
    return sp.simplify(sum(
        (e_invariant_C1[u, a] if map_s1 else int(u == a))
        * (e_invariant_C1[v, b] if map_s2bar else int(v == b))
        * e_invariant_mix_triplet_raw(h1, h2, u, v)
        for u in range(3)
        for v in range(3)
    ))


def count_e_invariant_ward_residuals(map_s1=False, map_s2bar=False):
    tf = e_invariant_fundamental_generators()
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
                                * build_e_invariant_mix_with_maps(p, h2, a, b, map_s1, map_s2bar)
                            )
                            value += (
                                tf[A][h2, p]
                                * build_e_invariant_mix_with_maps(h1, p, a, b, map_s1, map_s2bar)
                            )

                        for c in range(3):
                            value += (
                                tt[A][a, c]
                                * build_e_invariant_mix_with_maps(h1, h2, c, b, map_s1, map_s2bar)
                            )
                            # Bar[S2] transforms in the conjugate representation.
                            value -= (
                                tt[A][c, b]
                                * build_e_invariant_mix_with_maps(h1, h2, a, c, map_s1, map_s2bar)
                            )

                        if sp.simplify(value) != 0:
                            nonzero += 1

        by_generator.append(nonzero)
        total += nonzero

    return by_generator, total


def e_invariant_complex_block(d, first):
    z, zb, variables = [], [], []
    for c in range(d):
        x = sp.Symbol(f"x{first + 2*c}")
        y = sp.Symbol(f"x{first + 2*c + 1}")
        z.append((x + sp.I * y) / sp.sqrt(2))
        zb.append((x - sp.I * y) / sp.sqrt(2))
        variables.extend((x, y))
    return z, zb, variables


def e_invariant_bilinear(zb, z, matrix):
    return sp.expand(sum(
        zb[i] * matrix[i, j] * z[j]
        for i in range(len(z))
        for j in range(len(z))
    ))


def e_invariant_tensor_from_polynomial(expr, variables):
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


def e_invariant_get(tensor, *indices):
    return tensor.get(tuple(sorted(indices)), sp.S.Zero)


def e_invariant_ordered_weight(key):
    counts = Counter(key)
    value = factorial(4)
    for count in counts.values():
        value //= factorial(count)
    return value


def e_invariant_inner(left, right):
    return sp.simplify(sum(
        e_invariant_ordered_weight(key)
        * sp.conjugate(e_invariant_get(left, *key))
        * e_invariant_get(right, *key)
        for key in set(left) | set(right)
    ))


def e_invariant_scalar_cross_beta(target, other, n_real):
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
                        e_invariant_get(target, p, q, e, f) * e_invariant_get(other, e, f, r, s)
                        + e_invariant_get(other, p, q, e, f) * e_invariant_get(target, e, f, r, s)
                    )

        value = sp.simplify(value)
        if value != 0:
            generated[key] = value

    return generated


def e_invariant_project(generated, target):
    coefficient = sp.simplify(e_invariant_inner(target, generated) / e_invariant_inner(target, target))
    residual = {
        key: sp.simplify(
            e_invariant_get(generated, *key) - coefficient * e_invariant_get(target, *key)
        )
        for key in set(generated) | set(target)
    }
    residual = {key: value for key, value in residual.items() if value != 0}
    return coefficient, len(residual)


def recoupling_for_embedding(map_s1, map_s2bar):
    H, Hb, hv = e_invariant_complex_block(2, 1)
    S1, S1b, s1v = e_invariant_complex_block(3, 5)
    S2, S2b, s2v = e_invariant_complex_block(3, 11)
    variables = hv + s1v + s2v

    mix_poly = sp.expand(sum(
        build_e_invariant_mix_with_maps(h1, h2, a, b, map_s1, map_s2bar)
        * H[h1] * H[h2] * S1[a] * S2b[b]
        for h1 in range(2)
        for h2 in range(2)
        for a in range(3)
        for b in range(3)
    ))
    target = e_invariant_tensor_from_polynomial(mix_poly, variables)

    tf = e_invariant_fundamental_generators()
    tt = triplet_generators()

    quartics = {
        "lambdaH1Adj": sum(
            e_invariant_bilinear(Hb, H, tf[A]) * e_invariant_bilinear(S1b, S1, tt[A])
            for A in range(3)
        ),
        "lambdaH2Adj": sum(
            e_invariant_bilinear(Hb, H, tf[A]) * e_invariant_bilinear(S2b, S2, tt[A])
            for A in range(3)
        ),
        "lambda12Adj": sum(
            e_invariant_bilinear(S1b, S1, tt[A]) * e_invariant_bilinear(S2b, S2, tt[A])
            for A in range(3)
        ),
    }

    rows = {}
    for name, polynomial in quartics.items():
        other = e_invariant_tensor_from_polynomial(polynomial, variables)
        generated = e_invariant_scalar_cross_beta(target, other, len(variables))
        coefficient, residual = e_invariant_project(generated, target)
        rows[name] = {
            "coefficient": sp.sstr(coefficient),
            "residual_nonzero_components": residual,
            "closed": residual == 0,
        }

    return rows


def run_e_invariant_validation_cli():
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
        by_gen, total = count_e_invariant_ward_residuals(map_s1, map_s2bar)
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


# ---------------------------------------------------------------------------
# T3-E cross-recoupling validation
# Former source: ValidateT3ECrossRecoupling.py
# ---------------------------------------------------------------------------

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




e_cross_EPS2 = sp.Matrix([[0, 1], [-1, 0]])
e_cross_C1 = sp.Matrix([
    [0, 0, 1],
    [0, -1, 0],
    [1, 0, 0],
])

e_cross_D = {}
for a in range(3):
    for i in range(2):
        for j in range(2):
            e_cross_D[a, i, j] = sp.S.Zero
e_cross_D[0, 0, 0] = 1
e_cross_D[1, 0, 1] = 1 / sp.sqrt(2)
e_cross_D[1, 1, 0] = 1 / sp.sqrt(2)
e_cross_D[2, 1, 1] = 1


def e_cross_refined_mix_pair(h1, h2, p, q, r, s):
    terms = (
        e_cross_EPS2[h1, s] * e_cross_EPS2[h2, q] * e_cross_EPS2[p, r],
        e_cross_EPS2[h1, q] * e_cross_EPS2[h2, s] * e_cross_EPS2[p, r],
        e_cross_EPS2[h1, r] * e_cross_EPS2[h2, q] * e_cross_EPS2[p, s],
        e_cross_EPS2[h1, q] * e_cross_EPS2[h2, r] * e_cross_EPS2[p, s],
        e_cross_EPS2[h1, s] * e_cross_EPS2[h2, p] * e_cross_EPS2[q, r],
        e_cross_EPS2[h1, p] * e_cross_EPS2[h2, s] * e_cross_EPS2[q, r],
        e_cross_EPS2[h1, r] * e_cross_EPS2[h2, p] * e_cross_EPS2[q, s],
        e_cross_EPS2[h1, p] * e_cross_EPS2[h2, r] * e_cross_EPS2[q, s],
    )
    return sp.simplify(sum(terms) / 8)


def e_cross_mix_triplet_raw(h1, h2, a, b):
    return sp.simplify(sum(
        e_cross_D[a, p, q] * e_cross_D[b, r, s] * e_cross_refined_mix_pair(h1, h2, p, q, r, s)
        for p in range(2)
        for q in range(2)
        for r in range(2)
        for s in range(2)
    ))


def build_e_cross_gauge_covariant_mix(h1, h2, a, b):
    # Ward test selected C1 on the Bar[S2] triplet leg.
    return sp.simplify(sum(
        e_cross_C1[v, b] * e_cross_mix_triplet_raw(h1, h2, a, v)
        for v in range(3)
    ))


def e_cross_complex_block(d, first):
    z, zb, variables = [], [], []
    for c in range(d):
        x = sp.Symbol(f"x{first + 2*c}")
        y = sp.Symbol(f"x{first + 2*c + 1}")
        z.append((x + sp.I*y) / sp.sqrt(2))
        zb.append((x - sp.I*y) / sp.sqrt(2))
        variables.extend((x, y))
    return z, zb, variables


def e_cross_tensor_from_polynomial(expr, variables):
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


def e_cross_get(tensor, *indices):
    return tensor.get(tuple(sorted(indices)), sp.S.Zero)


def e_cross_ordered_weight(key):
    counts = Counter(key)
    value = factorial(4)
    for count in counts.values():
        value //= factorial(count)
    return value


def e_cross_inner(left, right):
    return sp.simplify(sum(
        e_cross_ordered_weight(key)
        * sp.conjugate(e_cross_get(left, *key))
        * e_cross_get(right, *key)
        for key in set(left) | set(right)
    ))


def e_cross_scalar_cross_beta(target, other, n_real):
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
                        e_cross_get(target, p, q, e, f) * e_cross_get(other, e, f, r, s)
                        + e_cross_get(other, p, q, e, f) * e_cross_get(target, e, f, r, s)
                    )
        value = sp.simplify(value)
        if value != 0:
            generated[key] = value
    return generated


def e_cross_project(generated, target):
    coefficient = sp.simplify(e_cross_inner(target, generated) / e_cross_inner(target, target))
    residual = {
        key: sp.simplify(e_cross_get(generated, *key) - coefficient * e_cross_get(target, *key))
        for key in set(generated) | set(target)
    }
    residual = {k: v for k, v in residual.items() if v != 0}
    return coefficient, residual


def run_e_cross_validation_cli():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    H, Hb, hv = e_cross_complex_block(2, 1)
    S1, S1b, s1v = e_cross_complex_block(3, 5)
    S2, S2b, s2v = e_cross_complex_block(3, 11)
    variables = hv + s1v + s2v

    mix_poly = sp.expand(sum(
        build_e_cross_gauge_covariant_mix(h1, h2, a, b)
        * H[h1] * H[h2] * S1[a] * S2b[b]
        for h1 in range(2)
        for h2 in range(2)
        for a in range(3)
        for b in range(3)
    ))

    # Exact Cross operator in spherical basis.
    same_orientation = (
        sum(S1b[a] * e_cross_C1[a, b] * S2b[b] for a in range(3) for b in range(3))
        * sum(S1[a] * e_cross_C1[a, b] * S2[b] for a in range(3) for b in range(3))
    )
    opposite_orientation = (
        sum(S1b[a] * S2[a] for a in range(3))
        * sum(S1[a] * S2b[a] for a in range(3))
    )
    cross_poly = sp.expand(same_orientation + opposite_orientation)

    target = e_cross_tensor_from_polynomial(mix_poly, variables)
    cross = e_cross_tensor_from_polynomial(cross_poly, variables)
    generated = e_cross_scalar_cross_beta(target, cross, len(variables))
    coefficient, residual = e_cross_project(generated, target)

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


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Consolidated group-factor research/validation driver."
    )
    parser.add_argument("mode", choices=[
        "diagnose-bc-pseudoreal",
        "bc-ward",
        "bc-tensor",
        "e-invariant",
        "e-cross",
    ])
    args, remaining = parser.parse_known_args()
    sys.argv = [sys.argv[0], *remaining]

    if args.mode == "diagnose-bc-pseudoreal":
        run_bc_diagonal_diagnostic_cli()
    elif args.mode == "bc-ward":
        run_bc_ward_validation_cli()
    elif args.mode == "bc-tensor":
        run_bc_tensor_validation_cli()
    elif args.mode == "e-invariant":
        run_e_invariant_validation_cli()
    elif args.mode == "e-cross":
        run_e_cross_validation_cli()


if __name__ == "__main__":
    main()
