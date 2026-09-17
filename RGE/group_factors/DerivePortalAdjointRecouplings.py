from __future__ import annotations

"""Derive physical adjoint-squared recouplings in portal-singlet beta functions.

Targets:
    beta_lambdaH1  <- (lambdaH1Adj)^2
    beta_lambdaH2  <- (lambdaH2Adj)^2
    beta_lambda12  <- (lambda12Adj)^2

No A--E coefficient is used as input.

The calculation is performed in an orthonormal spherical SU(2) basis, then
converted to the fully symmetric real-scalar quartic tensor convention

    V4 = (1/4!) lambda_abcd phi_a phi_b phi_c phi_d.

For each external field sector, the generated one-loop tensor is decomposed
onto the COMPLETE invariant basis:
  H-S: {singlet, adjoint}
  S1-S2 for 2x2: {singlet, adjoint}
  S1-S2 for 3x3: {singlet, adjoint, Cross}

Only components with the target external field content are retained before
the decomposition; diagrams generating H^4, S^4, etc. belong to different
beta functions and are intentionally excluded from the residual.
"""

import argparse
from collections import Counter
from itertools import combinations_with_replacement
from math import factorial
from pathlib import Path
import json

import sympy as sp


MODEL_DIMS = {
    "A": (1, 3, 2),
    "B": (2, 2, 1),
    "C": (2, 2, 3),
    "D": (3, 1, 2),
    "E": (3, 3, 2),
}


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
    value = factorial(4)
    for count in counts.values():
        value //= factorial(count)
    return value


def inner(A, B):
    return sp.simplify(sum(
        ordered_weight(key)
        * sp.conjugate(A.get(key, 0))
        * B.get(key, 0)
        for key in set(A) | set(B)
    ))


def pair_maps(tensor, n):
    maps = {}
    for a in range(n):
        for b in range(a, n):
            vec = {}
            for e in range(n):
                for f in range(n):
                    value = tget(tensor, a, b, e, f)
                    if value != 0:
                        vec[(e, f)] = value
            maps[(a, b)] = vec
    return maps


def sparse_dot(left, right):
    if len(left) > len(right):
        left, right = right, left
    return sp.simplify(sum(
        value * right.get(key, 0)
        for key, value in left.items()
    ))


def square_loop_tensor(A, n):
    """One-loop scalar tensor proportional to one coupling squared."""
    pairs = pair_maps(A, n)

    def pk(a, b):
        return (a, b) if a <= b else (b, a)

    out = {}
    for key in combinations_with_replacement(range(n), 4):
        a, b, c, d = key
        value = (
            sparse_dot(pairs[pk(a, b)], pairs[pk(c, d)])
            + sparse_dot(pairs[pk(a, c)], pairs[pk(b, d)])
            + sparse_dot(pairs[pk(a, d)], pairs[pk(b, c)])
        )
        value = sp.simplify(value)
        if value != 0:
            out[key] = value
    return out


def restrict_sector(tensor, groups, required):
    """Keep entries with exactly the requested number of real legs per field."""
    out = {}
    for key, value in tensor.items():
        counts = {name: 0 for name in groups}
        valid = True
        for idx in key:
            found = False
            for name, indices in groups.items():
                if idx in indices:
                    counts[name] += 1
                    found = True
                    break
            if not found:
                valid = False
                break
        if valid and counts == required:
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
        {
            name: sp.simplify(value)
            for name, value in zip(names, coefficients)
        },
        residual,
    )


def build_model(d1, d2):
    H, Hb, hv = complex_block(2, "H")
    S1, S1b, s1v = complex_block(d1, "S1")
    S2, S2b, s2v = complex_block(d2, "S2")
    variables = hv + s1v + s2v

    groups = {
        "H": set(range(0, len(hv))),
        "S1": set(range(len(hv), len(hv) + len(s1v))),
        "S2": set(range(len(hv) + len(s1v), len(variables))),
    }

    tH = su2_generators(2)
    t1 = su2_generators(d1)
    t2 = su2_generators(d2)

    poly = {
        "H1": dot(Hb, H) * dot(S1b, S1),
        "H2": dot(Hb, H) * dot(S2b, S2),
        "12": dot(S1b, S1) * dot(S2b, S2),
    }

    if d1 > 1:
        poly["H1Adj"] = sum(
            bilinear(Hb, H, tH[A]) * bilinear(S1b, S1, t1[A])
            for A in range(3)
        )
    if d2 > 1:
        poly["H2Adj"] = sum(
            bilinear(Hb, H, tH[A]) * bilinear(S2b, S2, t2[A])
            for A in range(3)
        )
    if d1 > 1 and d2 > 1:
        poly["12Adj"] = sum(
            bilinear(S1b, S1, t1[A]) * bilinear(S2b, S2, t2[A])
            for A in range(3)
        )

    if d1 == 3 and d2 == 3:
        C = charge_conjugation_metric(3)
        poly["12Cross"] = sp.expand(
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

    tensor = {
        name: tensor_from_polynomial(expr, variables)
        for name, expr in poly.items()
    }

    return variables, groups, tensor


def derive_model(model):
    d1, d2, _ = MODEL_DIMS[model]
    variables, groups, tensor = build_model(d1, d2)
    n = len(variables)
    rows = []

    jobs = []
    if "H1Adj" in tensor:
        basis = {"H1": tensor["H1"], "H1Adj": tensor["H1Adj"]}
        jobs.append((
            "beta_lambdaH1",
            "lambdaH1Adj^2",
            "H1Adj",
            basis,
            {"H": 2, "S1": 2, "S2": 0},
        ))
    if "H2Adj" in tensor:
        basis = {"H2": tensor["H2"], "H2Adj": tensor["H2Adj"]}
        jobs.append((
            "beta_lambdaH2",
            "lambdaH2Adj^2",
            "H2Adj",
            basis,
            {"H": 2, "S1": 0, "S2": 2},
        ))
    if "12Adj" in tensor:
        basis = {"lambda12": tensor["12"], "lambda12Adj": tensor["12Adj"]}
        if "12Cross" in tensor:
            basis["lambda12Cross"] = tensor["12Cross"]
        jobs.append((
            "beta_lambda12",
            "lambda12Adj^2",
            "12Adj",
            basis,
            {"H": 0, "S1": 2, "S2": 2},
        ))

    for beta, source, source_tensor, basis, sector in jobs:
        generated = square_loop_tensor(tensor[source_tensor], n)
        generated = restrict_sector(generated, groups, sector)
        coefficients, residual = decompose(generated, basis)
        rows.append({
            "beta": beta,
            "source": source,
            "coefficients": {
                name: sp.sstr(value)
                for name, value in coefficients.items()
            },
            "residual_nonzero_components": len(residual),
            "closed": len(residual) == 0,
        })

    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "output/group_factors/portal_adjoint_recouplings.json"
        ),
    )
    args = parser.parse_args()

    payload = {}
    overall = True

    for model in MODEL_DIMS:
        rows = derive_model(model)
        payload[model] = rows
        print(f"T3-{model}")
        if not rows:
            print("  no adjoint portal operator in this model")
            continue

        for row in rows:
            coeffs = ", ".join(
                f"{name}={value}"
                for name, value in row["coefficients"].items()
            )
            print(
                f"  {row['beta']:<15} {row['source']:<18} "
                f"{coeffs}; residual={row['residual_nonzero_components']}"
            )
            overall &= row["closed"]

    result = {
        "status": "Success" if overall else "ResidualNonzero",
        "models": payload,
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2), encoding="utf-8")

    print()
    print(f"JSON summary: {args.output}")
    print(f"OVERALL:      {result['status']}")
    return 0 if overall else 1


if __name__ == "__main__":
    raise SystemExit(main())
