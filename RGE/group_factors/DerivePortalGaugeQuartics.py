from __future__ import annotations

"""Derive the pure-gauge g^4 contributions to beta_lambdaH1/H2/lambda12.

The calculation uses the standard real-scalar quartic gauge structure

    beta_abcd |_(g^4)
      = 3 * sum_{A,B} [
          {Theta^A,Theta^B}_{ab} {Theta^A,Theta^B}_{cd}
        + {Theta^A,Theta^B}_{ac} {Theta^A,Theta^B}_{bd}
        + {Theta^A,Theta^B}_{ad} {Theta^A,Theta^B}_{bc}
      ],

where Theta are the real antisymmetric scalar generators INCLUDING the gauge
coupling.  This normalization reproduces the familiar 9/4 g2^4 contribution
to a doublet-doublet singlet portal.

Gauge group:
    SU(2)_L x U(1)_Y

Conventions:
    Q = T3 + Y
    Y(H)=1/2
    Y(S1)=alpha/2
    Y(S2)=(alpha+2)/2

For each portal sector the generated tensor is decomposed onto the complete
operator basis:
    H-S:   {singlet, adjoint}
    S1-S2: {singlet, adjoint[, Cross for 3x3]}

No RGBeta coefficients are used as inputs.
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

g2, gY = sp.symbols("g2 gY", real=True)


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


def real_generator(T: sp.Matrix) -> sp.Matrix:
    """Complex Hermitian T -> real antisymmetric generator in (R...,I...) order."""
    d = T.rows
    Re = T.applyfunc(sp.re)
    Im = T.applyfunc(sp.im)
    out = sp.zeros(2 * d)
    out[:d, :d] = -Im
    out[:d, d:] = -Re
    out[d:, :d] = Re
    out[d:, d:] = -Im
    return out


def block_diag(*blocks):
    return sp.diag(*blocks)


def charge_conjugation_metric(d: int):
    j = sp.Rational(d - 1, 2)
    ms = [j - i for i in range(d)]
    C = sp.zeros(d)
    for i, m in enumerate(ms):
        for k, mp in enumerate(ms):
            if mp == -m:
                C[i, k] = (-1) ** int(j - m)
    return C


def complex_block(d: int, symbols):
    z, zb = [], []
    R = symbols[:d]
    I = symbols[d:]
    for i in range(d):
        z.append((R[i] + sp.I * I[i]) / sp.sqrt(2))
        zb.append((R[i] - sp.I * I[i]) / sp.sqrt(2))
    return z, zb


def dot(a, b):
    return sp.expand(sum(x * y for x, y in zip(a, b)))


def bilinear(zb, z, T):
    return sp.expand(sum(
        zb[i] * T[i, j] * z[j]
        for i in range(len(z))
        for j in range(len(z))
    ))


def tensor_from_polynomial(expr, variables):
    poly = sp.Poly(sp.expand(expr), *variables)
    out = {}
    for powers, coefficient in poly.terms():
        if sum(powers) != 4:
            continue
        key = []
        mult = 1
        for i, power in enumerate(powers):
            key.extend([i] * power)
            mult *= factorial(power)
        value = sp.simplify(coefficient * mult)
        if value != 0:
            out[tuple(sorted(key))] = value
    return out


def ordered_weight(key):
    counts = Counter(key)
    value = factorial(4)
    for n in counts.values():
        value //= factorial(n)
    return value


def inner(A, B):
    return sp.simplify(sum(
        ordered_weight(k) * sp.conjugate(A.get(k, 0)) * B.get(k, 0)
        for k in set(A) | set(B)
    ))


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
    gram = sp.Matrix([[inner(A, B) for B in tensors] for A in tensors])
    rhs = sp.Matrix([inner(A, generated) for A in tensors])
    coeffs = gram.LUsolve(rhs)

    residual = {}
    keys = set(generated)
    for T in tensors:
        keys |= set(T)

    for key in keys:
        value = generated.get(key, 0)
        for c, T in zip(coeffs, tensors):
            value -= c * T.get(key, 0)
        value = sp.simplify(value)
        if value != 0:
            residual[key] = value

    return (
        {name: sp.factor(sp.simplify(c)) for name, c in zip(names, coeffs)},
        residual,
    )


def gauge_tensor(generators, n):
    """Pure gauge quartic tensor with weighted generators."""
    anticom = {}
    for i, A in enumerate(generators):
        for j, B in enumerate(generators):
            anticom[i, j] = sp.expand(A * B + B * A)

    out = {}
    for key in combinations_with_replacement(range(n), 4):
        a, b, c, d = key
        value = 0
        for i in range(len(generators)):
            for j in range(len(generators)):
                X = anticom[i, j]
                value += (
                    X[a, b] * X[c, d]
                    + X[a, c] * X[b, d]
                    + X[a, d] * X[b, c]
                )
        value = sp.factor(3 * sp.simplify(value))
        if value != 0:
            out[key] = value
    return out


def build_model(d1, d2, alpha):
    nH, n1, n2 = 4, 2 * d1, 2 * d2
    n = nH + n1 + n2
    variables = sp.symbols(f"x0:{n}", real=True)

    Hvars = list(variables[:nH])
    S1vars = list(variables[nH:nH+n1])
    S2vars = list(variables[nH+n1:])

    H, Hb = complex_block(2, Hvars)
    S1, S1b = complex_block(d1, S1vars)
    S2, S2b = complex_block(d2, S2vars)

    tH = su2_generators(2)
    t1 = su2_generators(d1)
    t2 = su2_generators(d2)

    zeroH = sp.zeros(nH)
    zero1 = sp.zeros(n1)
    zero2 = sp.zeros(n2)

    weighted = []

    for A in range(3):
        weighted.append(
            g2 * block_diag(
                real_generator(tH[A]),
                real_generator(t1[A]),
                real_generator(t2[A]),
            )
        )

    YH = sp.Rational(1, 2)
    Y1 = sp.Rational(alpha, 2)
    Y2 = sp.Rational(alpha + 2, 2)
    weighted.append(
        gY * block_diag(
            real_generator(YH * sp.eye(2)),
            real_generator(Y1 * sp.eye(d1)),
            real_generator(Y2 * sp.eye(d2)),
        )
    )

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
                sum(S1b[a] * C[a, b] * S2b[b] for a in range(3) for b in range(3))
                * sum(S1[a] * C[a, b] * S2[b] for a in range(3) for b in range(3))
            )
            + dot(S1b, S2) * dot(S1, S2b)
        )

    tensors = {
        name: tensor_from_polynomial(expr, variables)
        for name, expr in poly.items()
    }

    bases = {
        "H1": {"lambdaH1": tensors["H1"]},
        "H2": {"lambdaH2": tensors["H2"]},
        "12": {"lambda12": tensors["12"]},
    }
    if "H1Adj" in tensors:
        bases["H1"]["lambdaH1Adj"] = tensors["H1Adj"]
    if "H2Adj" in tensors:
        bases["H2"]["lambdaH2Adj"] = tensors["H2Adj"]
    if "12Adj" in tensors:
        bases["12"]["lambda12Adj"] = tensors["12Adj"]
    if "12Cross" in tensors:
        bases["12"]["lambda12Cross"] = tensors["12Cross"]

    groups = {
        "H": set(range(0, nH)),
        "S1": set(range(nH, nH+n1)),
        "S2": set(range(nH+n1, n)),
    }

    return n, weighted, tensors, bases, groups


def derive(model, alpha):
    d1, d2, _ = MODEL_DIMS[model]
    n, weighted, tensors, bases, groups = build_model(d1, d2, alpha)
    generated = gauge_tensor(weighted, n)

    sectors = {
        "beta_lambdaH1": ("H1", {"H": 2, "S1": 2, "S2": 0}),
        "beta_lambdaH2": ("H2", {"H": 2, "S1": 0, "S2": 2}),
        "beta_lambda12": ("12", {"H": 0, "S1": 2, "S2": 2}),
    }

    rows = {}
    for beta, (sector, required) in sectors.items():
        restricted = restrict_sector(generated, groups, required)
        coeffs, residual = decompose(restricted, bases[sector])
        rows[beta] = {
            "coefficients": {
                name: sp.sstr(sp.factor(value))
                for name, value in coeffs.items()
            },
            "residual_nonzero_components": len(residual),
            "closed": len(residual) == 0,
        }
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/group_factors/portal_gauge_quartics.json"),
    )
    args = parser.parse_args()

    payload = {}
    overall = True

    for model in MODEL_DIMS:
        print(f"T3-{model}")
        payload[model] = {}
        for alpha in (-2, -1, 0, 1, 2):
            rows = derive(model, alpha)
            payload[model][str(alpha)] = rows

            print(f"  alpha={alpha:+d}")
            for beta, row in rows.items():
                coeffs = ", ".join(
                    f"{name}={value}"
                    for name, value in row["coefficients"].items()
                )
                print(
                    f"    {beta:<15} {coeffs}; "
                    f"residual={row['residual_nonzero_components']}"
                )
                overall &= row["closed"]

    result = {
        "status": "Success" if overall else "ResidualNonzero",
        "formula": (
            "3 sum_AB [ {ThetaA,ThetaB}_ab {ThetaA,ThetaB}_cd + 2 perms ]"
        ),
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
