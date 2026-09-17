from __future__ import annotations

"""Derive the remaining T3-E scalar recouplings in beta_lambdaH1/H2/lambda12.

This is the companion to DerivePortalAdjointRecouplings.py.

It derives, in the physical gauge-covariant component basis, the scalar-loop
contributions involving the T3-E Cross quartic and triplet self-adjoint
quartics.  No RGBeta coefficient is used as input.

External operator bases
------------------------
H1 sector:
    O_H1    = (H†H)(S1†S1)
    O_H1Adj = sum_A (H† t^A H)(S1† T^A S1)

H2 sector:
    O_H2    = (H†H)(S2†S2)
    O_H2Adj = sum_A (H† t^A H)(S2† T^A S2)

12 sector:
    O_12      = (S1†S1)(S2†S2)
    O_12Adj   = sum_A (S1† T^A S1)(S2† T^A S2)
    O_12Cross = (BarS1^T C BarS2)(S1^T C S2)
              + (BarS1^T S2)(S1^T BarS2)

Triplet self-adjoint quartics:
    O_SiAdj = 1/2 sum_A (Si† T^A Si)^2

The one-loop scalar quartic tensor convention is
    V4 = (1/4!) lambda_abcd phi_a phi_b phi_c phi_d
with
    beta_abcd|_{lambda^2}
      = lambda_abef lambda_efcd
      + lambda_acef lambda_efbd
      + lambda_adef lambda_efbc.

For distinct couplings X,Y both orderings contribute.
"""

import argparse
from collections import Counter
from itertools import combinations_with_replacement
from math import factorial
from pathlib import Path
import json

import sympy as sp


def su2_generators(d: int):
    j = sp.Rational(d - 1, 2)
    ms = [j - i for i in range(d)]
    jp = sp.zeros(d)

    for col, m in enumerate(ms):
        mp = m + 1
        if mp in ms:
            row = ms.index(mp)
            jp[row, col] = sp.sqrt((j - m) * (j + m + 1))

    jm = jp.T
    return (
        sp.simplify((jp + jm) / 2),
        sp.simplify((jp - jm) / (2 * sp.I)),
        sp.diag(*ms),
    )


def charge_conjugation_metric(d: int):
    j = sp.Rational(d - 1, 2)
    ms = [j - i for i in range(d)]
    C = sp.zeros(d)
    for i, m in enumerate(ms):
        for k, mp in enumerate(ms):
            if mp == -m:
                C[i, k] = (-1) ** int(j - m)
    return C


def complex_block(d: int, prefix: str):
    z, zb, variables = [], [], []
    for c in range(d):
        x = sp.Symbol(f"{prefix}R{c}")
        y = sp.Symbol(f"{prefix}I{c}")
        z.append((x + sp.I * y) / sp.sqrt(2))
        zb.append((x - sp.I * y) / sp.sqrt(2))
        variables.extend((x, y))
    return z, zb, variables


def dot(a, b):
    return sp.expand(sum(x * y for x, y in zip(a, b)))


def bilinear(zb, z, M):
    return sp.expand(sum(
        zb[i] * M[i, j] * z[j]
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
        for i, power in enumerate(powers):
            key.extend([i] * power)
            multiplicity *= factorial(power)

        value = sp.simplify(coefficient * multiplicity)
        if value != 0:
            entries[tuple(sorted(key))] = value

    return entries


def tget(tensor, a, b, c, d):
    return tensor.get(tuple(sorted((a, b, c, d))), sp.S.Zero)


def ordered_weight(key):
    counts = Counter(key)
    result = factorial(4)
    for n in counts.values():
        result //= factorial(n)
    return result


def inner(A, B):
    return sp.simplify(sum(
        ordered_weight(key)
        * sp.conjugate(A.get(key, 0))
        * B.get(key, 0)
        for key in set(A) | set(B)
    ))


def pair_maps(tensor, n):
    out = {}
    for a in range(n):
        for b in range(a, n):
            vec = {}
            for e in range(n):
                for f in range(n):
                    value = tget(tensor, a, b, e, f)
                    if value != 0:
                        vec[(e, f)] = value
            out[(a, b)] = vec
    return out


def sparse_dot(left, right):
    if len(left) > len(right):
        left, right = right, left
    return sp.simplify(sum(
        value * right.get(key, 0)
        for key, value in left.items()
    ))


def scalar_cross_tensor(A, B, n, same=False):
    pa = pair_maps(A, n)
    pb = pa if same else pair_maps(B, n)

    def pk(a, b):
        return (a, b) if a <= b else (b, a)

    out = {}

    for key in combinations_with_replacement(range(n), 4):
        a, b, c, d = key

        channels = (
            ((a, b), (c, d)),
            ((a, c), (b, d)),
            ((a, d), (b, c)),
        )

        value = 0
        for left_pair, right_pair in channels:
            lp = pk(*left_pair)
            rp = pk(*right_pair)
            value += sparse_dot(pa[lp], pb[rp])
            if not same:
                # Distinct couplings require the second ordering explicitly.
                # It is NOT generally equal term-by-term to the first one.
                value += sparse_dot(pb[lp], pa[rp])

        value = sp.simplify(value)
        if value != 0:
            out[key] = value

    return out


def restrict_sector(tensor, groups, required):
    out = {}
    for key, value in tensor.items():
        counts = {name: 0 for name in groups}
        for idx in key:
            for name, indices in groups.items():
                if idx in indices:
                    counts[name] += 1
                    break
        if counts == required:
            out[key] = value
    return out


def decompose(generated, basis):
    names = list(basis)
    tensors = [basis[name] for name in names]

    gram = sp.Matrix([
        [inner(A, B) for B in tensors]
        for A in tensors
    ])
    rhs = sp.Matrix([inner(A, generated) for A in tensors])
    coefficients = gram.LUsolve(rhs)

    residual = {}
    keys = set(generated)
    for tensor in tensors:
        keys |= set(tensor)

    for key in keys:
        value = generated.get(key, 0)
        for coefficient, tensor in zip(coefficients, tensors):
            value -= coefficient * tensor.get(key, 0)
        value = sp.simplify(value)
        if value != 0:
            residual[key] = value

    return (
        {name: sp.simplify(value) for name, value in zip(names, coefficients)},
        residual,
    )


def build_e_operators():
    H, Hb, hv = complex_block(2, "H")
    S1, S1b, s1v = complex_block(3, "S1")
    S2, S2b, s2v = complex_block(3, "S2")
    variables = hv + s1v + s2v

    groups = {
        "H": set(range(0, len(hv))),
        "S1": set(range(len(hv), len(hv) + len(s1v))),
        "S2": set(range(len(hv) + len(s1v), len(variables))),
    }

    tH = su2_generators(2)
    t1 = su2_generators(3)
    t2 = su2_generators(3)
    C = charge_conjugation_metric(3)

    polynomials = {
        "H1": dot(Hb, H) * dot(S1b, S1),
        "H2": dot(Hb, H) * dot(S2b, S2),
        "12": dot(S1b, S1) * dot(S2b, S2),
        "H1Adj": sum(
            bilinear(Hb, H, tH[A]) * bilinear(S1b, S1, t1[A])
            for A in range(3)
        ),
        "H2Adj": sum(
            bilinear(Hb, H, tH[A]) * bilinear(S2b, S2, t2[A])
            for A in range(3)
        ),
        "12Adj": sum(
            bilinear(S1b, S1, t1[A]) * bilinear(S2b, S2, t2[A])
            for A in range(3)
        ),
        "S1Adj": sp.Rational(1, 2) * sum(
            bilinear(S1b, S1, t1[A]) ** 2
            for A in range(3)
        ),
        "S2Adj": sp.Rational(1, 2) * sum(
            bilinear(S2b, S2, t2[A]) ** 2
            for A in range(3)
        ),
    }

    polynomials["12Cross"] = sp.expand(
        (
            sum(
                S1b[a] * C[a, b] * S2b[b]
                for a in range(3) for b in range(3)
            )
            * sum(
                S1[a] * C[a, b] * S2[b]
                for a in range(3) for b in range(3)
            )
        )
        + dot(S1b, S2) * dot(S1, S2b)
    )

    tensors = {
        name: tensor_from_polynomial(expr, variables)
        for name, expr in polynomials.items()
    }

    bases = {
        "H1": {
            "lambdaH1": tensors["H1"],
            "lambdaH1Adj": tensors["H1Adj"],
        },
        "H2": {
            "lambdaH2": tensors["H2"],
            "lambdaH2Adj": tensors["H2Adj"],
        },
        "12": {
            "lambda12": tensors["12"],
            "lambda12Adj": tensors["12Adj"],
            "lambda12Cross": tensors["12Cross"],
        },
    }

    return variables, groups, tensors, bases


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "output/group_factors/e_portal_cross_scalar_recouplings.json"
        ),
    )
    args = parser.parse_args()

    variables, groups, tensors, bases = build_e_operators()
    n = len(variables)

    jobs = [
        ("beta_lambdaH1", "lambdaH2*lambda12Cross",
         "H2", "12Cross", False, "H1",
         {"H": 2, "S1": 2, "S2": 0}),
        ("beta_lambdaH2", "lambdaH1*lambda12Cross",
         "H1", "12Cross", False, "H2",
         {"H": 2, "S1": 0, "S2": 2}),
        ("beta_lambda12", "lambda12*lambda12Cross",
         "12", "12Cross", False, "12",
         {"H": 0, "S1": 2, "S2": 2}),
        ("beta_lambda12", "lambda12Cross^2",
         "12Cross", "12Cross", True, "12",
         {"H": 0, "S1": 2, "S2": 2}),
        ("beta_lambda12", "lambda12Cross*lambdaS1Adj",
         "12Cross", "S1Adj", False, "12",
         {"H": 0, "S1": 2, "S2": 2}),
        ("beta_lambda12", "lambda12Cross*lambdaS2Adj",
         "12Cross", "S2Adj", False, "12",
         {"H": 0, "S1": 2, "S2": 2}),
        ("beta_lambdaH1", "lambdaH1*lambdaS1Adj",
         "H1", "S1Adj", False, "H1",
         {"H": 2, "S1": 2, "S2": 0}),
        ("beta_lambdaH2", "lambdaH2*lambdaS2Adj",
         "H2", "S2Adj", False, "H2",
         {"H": 2, "S1": 0, "S2": 2}),
        ("beta_lambda12", "lambda12*lambdaS1Adj",
         "12", "S1Adj", False, "12",
         {"H": 0, "S1": 2, "S2": 2}),
        ("beta_lambda12", "lambda12*lambdaS2Adj",
         "12", "S2Adj", False, "12",
         {"H": 0, "S1": 2, "S2": 2}),
    ]

    rows = []
    overall = True

    print("T3-E remaining scalar recouplings")

    for beta, source, left, right, same, sector, required in jobs:
        generated = scalar_cross_tensor(
            tensors[left], tensors[right], n, same=same
        )
        generated = restrict_sector(generated, groups, required)
        coefficients, residual = decompose(generated, bases[sector])

        row = {
            "beta": beta,
            "source": source,
            "coefficients": {
                name: sp.sstr(value)
                for name, value in coefficients.items()
            },
            "residual_nonzero_components": len(residual),
            "closed": len(residual) == 0,
        }
        rows.append(row)
        overall &= row["closed"]

        coeffs = ", ".join(
            f"{name}={value}"
            for name, value in row["coefficients"].items()
        )
        print(
            f"  {beta:<15} {source:<32} "
            f"{coeffs}; residual={len(residual)}"
        )

    payload = {
        "status": "Success" if overall else "ResidualNonzero",
        "model": "E",
        "rows": rows,
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print()
    print(f"JSON summary: {args.output}")
    print(f"OVERALL:      {payload['status']}")
    return 0 if overall else 1


if __name__ == "__main__":
    raise SystemExit(main())
