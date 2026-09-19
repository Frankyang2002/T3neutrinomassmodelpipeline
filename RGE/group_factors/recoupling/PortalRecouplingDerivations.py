from __future__ import annotations

"""Portal-quartic recoupling derivations.

The three derivations remain logically independent and retain their original calculation bodies; only their source packaging is combined."""

import argparse
from itertools import combinations_with_replacement
from pathlib import Path
import sys
import json
import sympy as sp
from RGE.group_factors.core.QuarticTensorAlgebra import pair_maps, restrict_sector, sparse_dot, tensor_from_polynomial, tensor_inner
from RGE.group_factors.core.QuarticTensorAlgebra import restrict_sector, tensor_from_polynomial, tensor_inner

# ---------------------------------------------------------------------------
# T3-E portal cross-scalar recoupling derivation
# Former source: DeriveEPortalCrossScalarRecouplings.py
# ---------------------------------------------------------------------------

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



ecrossportal_PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(ecrossportal_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(ecrossportal_PROJECT_ROOT))



def ecrossportal_su2_generators(d: int):
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


def ecrossportal_charge_conjugation_metric(d: int):
    j = sp.Rational(d - 1, 2)
    ms = [j - i for i in range(d)]
    C = sp.zeros(d)
    for i, m in enumerate(ms):
        for k, mp in enumerate(ms):
            if mp == -m:
                C[i, k] = (-1) ** int(j - m)
    return C


def ecrossportal_complex_block(d: int, prefix: str):
    z, zb, variables = [], [], []
    for c in range(d):
        x = sp.Symbol(f"{prefix}R{c}")
        y = sp.Symbol(f"{prefix}I{c}")
        z.append((x + sp.I * y) / sp.sqrt(2))
        zb.append((x - sp.I * y) / sp.sqrt(2))
        variables.extend((x, y))
    return z, zb, variables


def ecrossportal_dot(a, b):
    return sp.expand(sum(x * y for x, y in zip(a, b)))


def ecrossportal_bilinear(zb, z, M):
    return sp.expand(sum(
        zb[i] * M[i, j] * z[j]
        for i in range(len(z))
        for j in range(len(z))
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


def ecrossportal_decompose(generated, basis):
    names = list(basis)
    tensors = [basis[name] for name in names]

    gram = sp.Matrix([
        [tensor_inner(A, B) for B in tensors]
        for A in tensors
    ])
    rhs = sp.Matrix([tensor_inner(A, generated) for A in tensors])
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
    H, Hb, hv = ecrossportal_complex_block(2, "H")
    S1, S1b, s1v = ecrossportal_complex_block(3, "S1")
    S2, S2b, s2v = ecrossportal_complex_block(3, "S2")
    variables = hv + s1v + s2v

    groups = {
        "H": set(range(0, len(hv))),
        "S1": set(range(len(hv), len(hv) + len(s1v))),
        "S2": set(range(len(hv) + len(s1v), len(variables))),
    }

    tH = ecrossportal_su2_generators(2)
    t1 = ecrossportal_su2_generators(3)
    t2 = ecrossportal_su2_generators(3)
    C = ecrossportal_charge_conjugation_metric(3)

    polynomials = {
        "H1": ecrossportal_dot(Hb, H) * ecrossportal_dot(S1b, S1),
        "H2": ecrossportal_dot(Hb, H) * ecrossportal_dot(S2b, S2),
        "12": ecrossportal_dot(S1b, S1) * ecrossportal_dot(S2b, S2),
        "H1Adj": sum(
            ecrossportal_bilinear(Hb, H, tH[A]) * ecrossportal_bilinear(S1b, S1, t1[A])
            for A in range(3)
        ),
        "H2Adj": sum(
            ecrossportal_bilinear(Hb, H, tH[A]) * ecrossportal_bilinear(S2b, S2, t2[A])
            for A in range(3)
        ),
        "12Adj": sum(
            ecrossportal_bilinear(S1b, S1, t1[A]) * ecrossportal_bilinear(S2b, S2, t2[A])
            for A in range(3)
        ),
        "S1Adj": sp.Rational(1, 2) * sum(
            ecrossportal_bilinear(S1b, S1, t1[A]) ** 2
            for A in range(3)
        ),
        "S2Adj": sp.Rational(1, 2) * sum(
            ecrossportal_bilinear(S2b, S2, t2[A]) ** 2
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
        + ecrossportal_dot(S1b, S2) * ecrossportal_dot(S1, S2b)
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


def ecrossportal_main():
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
        coefficients, residual = ecrossportal_decompose(generated, bases[sector])

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


# ---------------------------------------------------------------------------
# Portal adjoint recoupling derivation
# Former source: DerivePortalAdjointRecouplings.py
# ---------------------------------------------------------------------------

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



adjportal_PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(adjportal_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(adjportal_PROJECT_ROOT))



adjportal_MODEL_DIMS = {
    "A": (1, 3, 2),
    "B": (2, 2, 1),
    "C": (2, 2, 3),
    "D": (3, 1, 2),
    "E": (3, 3, 2),
}


def adjportal_su2_generators(d: int):
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


def adjportal_charge_conjugation_metric(d: int):
    j = sp.Rational(d - 1, 2)
    ms = [j - i for i in range(d)]
    C = sp.zeros(d)
    for i, m in enumerate(ms):
        for k, mp in enumerate(ms):
            if mp == -m:
                C[i, k] = (-1) ** int(j - m)
    return C


def adjportal_complex_block(d: int, prefix: str):
    z, zb, variables = [], [], []
    for c in range(d):
        x = sp.Symbol(f"{prefix}R{c}")
        y = sp.Symbol(f"{prefix}I{c}")
        z.append((x + sp.I * y) / sp.sqrt(2))
        zb.append((x - sp.I * y) / sp.sqrt(2))
        variables.extend((x, y))
    return z, zb, variables


def adjportal_dot(a, b):
    return sp.expand(sum(x * y for x, y in zip(a, b)))


def adjportal_bilinear(zb, z, M):
    return sp.expand(sum(
        zb[i] * M[i, j] * z[j]
        for i in range(len(z))
        for j in range(len(z))
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


def adjportal_decompose(generated, basis):
    names = list(basis)
    tensors = [basis[name] for name in names]

    gram = sp.Matrix([
        [tensor_inner(A, B) for B in tensors]
        for A in tensors
    ])
    rhs = sp.Matrix([tensor_inner(A, generated) for A in tensors])
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


def adjportal_build_model(d1, d2):
    H, Hb, hv = adjportal_complex_block(2, "H")
    S1, S1b, s1v = adjportal_complex_block(d1, "S1")
    S2, S2b, s2v = adjportal_complex_block(d2, "S2")
    variables = hv + s1v + s2v

    groups = {
        "H": set(range(0, len(hv))),
        "S1": set(range(len(hv), len(hv) + len(s1v))),
        "S2": set(range(len(hv) + len(s1v), len(variables))),
    }

    tH = adjportal_su2_generators(2)
    t1 = adjportal_su2_generators(d1)
    t2 = adjportal_su2_generators(d2)

    poly = {
        "H1": adjportal_dot(Hb, H) * adjportal_dot(S1b, S1),
        "H2": adjportal_dot(Hb, H) * adjportal_dot(S2b, S2),
        "12": adjportal_dot(S1b, S1) * adjportal_dot(S2b, S2),
    }

    if d1 > 1:
        poly["H1Adj"] = sum(
            adjportal_bilinear(Hb, H, tH[A]) * adjportal_bilinear(S1b, S1, t1[A])
            for A in range(3)
        )
    if d2 > 1:
        poly["H2Adj"] = sum(
            adjportal_bilinear(Hb, H, tH[A]) * adjportal_bilinear(S2b, S2, t2[A])
            for A in range(3)
        )
    if d1 > 1 and d2 > 1:
        poly["12Adj"] = sum(
            adjportal_bilinear(S1b, S1, t1[A]) * adjportal_bilinear(S2b, S2, t2[A])
            for A in range(3)
        )

    if d1 == 3 and d2 == 3:
        C = adjportal_charge_conjugation_metric(3)
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
            + adjportal_dot(S1b, S2) * adjportal_dot(S1, S2b)
        )

    tensor = {
        name: tensor_from_polynomial(expr, variables)
        for name, expr in poly.items()
    }

    return variables, groups, tensor


def derive_model(model):
    d1, d2, _ = adjportal_MODEL_DIMS[model]
    variables, groups, tensor = adjportal_build_model(d1, d2)
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
        generated = restrict_sector(
            generated,
            groups,
            sector,
            reject_unassigned=True,
        )
        coefficients, residual = adjportal_decompose(generated, basis)
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


def adjportal_main():
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

    for model in adjportal_MODEL_DIMS:
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


# ---------------------------------------------------------------------------
# Portal gauge-quartic derivation
# Former source: DerivePortalGaugeQuartics.py
# ---------------------------------------------------------------------------

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



gaugeportal_PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(gaugeportal_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(gaugeportal_PROJECT_ROOT))



gaugeportal_MODEL_DIMS = {
    "A": (1, 3, 2),
    "B": (2, 2, 1),
    "C": (2, 2, 3),
    "D": (3, 1, 2),
    "E": (3, 3, 2),
}

g2, gY = sp.symbols("g2 gY", real=True)


def gaugeportal_su2_generators(d: int):
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


def gaugeportal_charge_conjugation_metric(d: int):
    j = sp.Rational(d - 1, 2)
    ms = [j - i for i in range(d)]
    C = sp.zeros(d)
    for i, m in enumerate(ms):
        for k, mp in enumerate(ms):
            if mp == -m:
                C[i, k] = (-1) ** int(j - m)
    return C


def gaugeportal_complex_block(d: int, symbols):
    z, zb = [], []
    R = symbols[:d]
    I = symbols[d:]
    for i in range(d):
        z.append((R[i] + sp.I * I[i]) / sp.sqrt(2))
        zb.append((R[i] - sp.I * I[i]) / sp.sqrt(2))
    return z, zb


def gaugeportal_dot(a, b):
    return sp.expand(sum(x * y for x, y in zip(a, b)))


def gaugeportal_bilinear(zb, z, T):
    return sp.expand(sum(
        zb[i] * T[i, j] * z[j]
        for i in range(len(z))
        for j in range(len(z))
    ))


def gaugeportal_decompose(generated, basis):
    names = list(basis)
    tensors = [basis[name] for name in names]
    gram = sp.Matrix([[tensor_inner(A, B) for B in tensors] for A in tensors])
    rhs = sp.Matrix([tensor_inner(A, generated) for A in tensors])
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


def gaugeportal_build_model(d1, d2, alpha):
    nH, n1, n2 = 4, 2 * d1, 2 * d2
    n = nH + n1 + n2
    variables = sp.symbols(f"x0:{n}", real=True)

    Hvars = list(variables[:nH])
    S1vars = list(variables[nH:nH+n1])
    S2vars = list(variables[nH+n1:])

    H, Hb = gaugeportal_complex_block(2, Hvars)
    S1, S1b = gaugeportal_complex_block(d1, S1vars)
    S2, S2b = gaugeportal_complex_block(d2, S2vars)

    tH = gaugeportal_su2_generators(2)
    t1 = gaugeportal_su2_generators(d1)
    t2 = gaugeportal_su2_generators(d2)

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
        "H1": gaugeportal_dot(Hb, H) * gaugeportal_dot(S1b, S1),
        "H2": gaugeportal_dot(Hb, H) * gaugeportal_dot(S2b, S2),
        "12": gaugeportal_dot(S1b, S1) * gaugeportal_dot(S2b, S2),
    }

    if d1 > 1:
        poly["H1Adj"] = sum(
            gaugeportal_bilinear(Hb, H, tH[A]) * gaugeportal_bilinear(S1b, S1, t1[A])
            for A in range(3)
        )
    if d2 > 1:
        poly["H2Adj"] = sum(
            gaugeportal_bilinear(Hb, H, tH[A]) * gaugeportal_bilinear(S2b, S2, t2[A])
            for A in range(3)
        )
    if d1 > 1 and d2 > 1:
        poly["12Adj"] = sum(
            gaugeportal_bilinear(S1b, S1, t1[A]) * gaugeportal_bilinear(S2b, S2, t2[A])
            for A in range(3)
        )

    if d1 == 3 and d2 == 3:
        C = gaugeportal_charge_conjugation_metric(3)
        poly["12Cross"] = sp.expand(
            (
                sum(S1b[a] * C[a, b] * S2b[b] for a in range(3) for b in range(3))
                * sum(S1[a] * C[a, b] * S2[b] for a in range(3) for b in range(3))
            )
            + gaugeportal_dot(S1b, S2) * gaugeportal_dot(S1, S2b)
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
    d1, d2, _ = gaugeportal_MODEL_DIMS[model]
    n, weighted, tensors, bases, groups = gaugeportal_build_model(d1, d2, alpha)
    generated = gauge_tensor(weighted, n)

    sectors = {
        "beta_lambdaH1": ("H1", {"H": 2, "S1": 2, "S2": 0}),
        "beta_lambdaH2": ("H2", {"H": 2, "S1": 0, "S2": 2}),
        "beta_lambda12": ("12", {"H": 0, "S1": 2, "S2": 2}),
    }

    rows = {}
    for beta, (sector, required) in sectors.items():
        restricted = restrict_sector(generated, groups, required)
        coeffs, residual = gaugeportal_decompose(restricted, bases[sector])
        rows[beta] = {
            "coefficients": {
                name: sp.sstr(sp.factor(value))
                for name, value in coeffs.items()
            },
            "residual_nonzero_components": len(residual),
            "closed": len(residual) == 0,
        }
    return rows


def gaugeportal_main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/group_factors/portal_gauge_quartics.json"),
    )
    args = parser.parse_args()

    payload = {}
    overall = True

    for model in gaugeportal_MODEL_DIMS:
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


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Consolidated group-factor research/validation driver."
    )
    parser.add_argument("mode", choices=[
        "e-cross",
        "adjoint",
        "gauge",
    ])
    args, remaining = parser.parse_known_args()
    sys.argv = [sys.argv[0], *remaining]

    if args.mode == "e-cross":
        ecrossportal_main()
    elif args.mode == "adjoint":
        adjportal_main()
    elif args.mode == "gauge":
        gaugeportal_main()


if __name__ == "__main__":
    main()
