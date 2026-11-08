"""Independent exact A(-4) scalar-only one-loop beta polynomial.

Mathematics: with 12 canonically normalised real scalar fields phi and
V_4=sum_k x_k P_k(phi), 16*pi^2 beta(V_4)|scalar = 1/2 Tr[(d^2 V_4/dphi^2)^2].
The ten P_k come directly from the *declared* SU(2) invariant potential.
This does not read RGBeta tensors, their quartic projection, or its beta
expressions, and does not certify a Matchete-to-RGBeta coupling map.

SU(2) generators: Pauli/2 for H, spin-1 weight basis for S2.
Complex scalars (H,S1,S2) = (Re+i Im)/sqrt(2), with the same 12-field
order as Numerical/running/Aminus4EFT1Couplings.py.
"""
from __future__ import annotations

from functools import lru_cache
from itertools import combinations_with_replacement

import sympy as sp

COORDINATES = (
    'lambdaH', 'lambdaS1', 'lambdaS2', 'lambdaS2Adj',
    'lambdaH1', 'lambdaH2', 'lambdaH2Adj', 'lambda12',
    'lambdaT3.real', 'lambdaT3.imag',
)


def _model():
    phi = sp.symbols('phi0:12', real=True)
    z = [(phi[2 * k] + sp.I * phi[2 * k + 1]) / sp.sqrt(2)
         for k in range(6)]
    h = sp.Matrix(z[:2]); s1 = z[2]; s2 = sp.Matrix(z[3:])
    hd = h.conjugate(); td = s2.conjugate()
    half = sp.Rational(1, 2)
    t_half = (
        sp.Matrix([[0, half], [half, 0]]),
        sp.Matrix([[0, -sp.I * half], [sp.I * half, 0]]),
        sp.diag(half, -half),
    )
    q = 1 / sp.sqrt(2)
    t_one = (
        sp.Matrix([[0, q, 0], [q, 0, q], [0, q, 0]]),
        sp.Matrix([[0, -sp.I*q, 0], [sp.I*q, 0, -sp.I*q],
                   [0, sp.I*q, 0]]),
        sp.diag(1, 0, -1),
    )
    hh = (hd.T * h)[0]
    s1bar_s1 = sp.conjugate(s1)*s1
    s2bar_s2 = (td.T * s2)[0]
    triplet_adj = sum(((td.T * t * s2)[0])**2 for t in t_one)
    portal_adj = sum((hd.T * t * h)[0] * (td.T * u * s2)[0]
                     for t, u in zip(t_half, t_one))
    # 2x2 -> 3 symmetric CG coefficients (1, 1/sqrt(2), 1);
    # the pair of off-diagonal HH components contributes sqrt(2).
    mix = (h[0]**2 * td[0] + sp.sqrt(2)*h[0]*h[1]*td[1]
           + h[1]**2*td[2])*s1
    p = (
        hh**2 / 2, s1bar_s1**2 / 2, s2bar_s2**2 / 2,
        triplet_adj / 2, hh * s1bar_s1, hh * s2bar_s2,
        portal_adj, s1bar_s1*s2bar_s2,
        mix + sp.conjugate(mix),
        sp.I*(mix - sp.conjugate(mix)),
    )
    return phi, tuple(sp.expand(x) for x in p), t_half, t_one, s2, td


def group_identities() -> dict:
    """Exact Pauli Fierz and complex spin-1 adjoint magnitude identity."""
    phi, ps, fundamental, spin1, t, td = _model()
    def delta(a, b): return int(a == b)
    max_nonzero = []
    for a in range(2):
        for b in range(2):
            for c in range(2):
                for d in range(2):
                    lhs = sum(x[a, b]*x[c, d] for x in fundamental)
                    rhs = sp.Rational(1, 2)*delta(a, d)*delta(b, c) - (
                        sp.Rational(1, 4)*delta(a, b)*delta(c, d))
                    max_nonzero.append(sp.simplify(lhs-rhs))
    triplet_spin = sum((td.T * g * t)[0]**2 for g in spin1)
    singlet = 2*t[0]*t[2] - t[1]**2
    expected = (td.T*t)[0]**2 - singlet*sp.conjugate(singlet)
    return {
        'fundamental_Fierz_exact': all(v == 0 for v in max_nonzero),
        'spin_one_adjoint_square_exact': sp.expand(triplet_spin-expected) == 0,
        'identity_spin_one': 'sum_A (S2_dagger T_A S2)^2 = (S2_dagger S2)^2 - abs(2*S2_1*S2_3 - S2_2^2)^2',
    }


@lru_cache(maxsize=1)
def beta_monomials_exact() -> dict:
    """Return all 55 exact scalar-loop beta monomials by Hessian projection.

    Constructs the 10 invariant quartic polynomials, selects 10 algebraically
    independent coefficients of the 12-variable polynomials and solves the
    resulting rational linear system. Every projected identity is *also*
    checked as a full 12-variable polynomial, not merely at ten sample points.
    """
    phi, basis, *_ = _model()
    polynomials = [sp.Poly(v, *phi) for v in basis]
    monomials = sorted(set().union(*(set(p.monoms()) for p in polynomials)))
    coefficient_matrix = sp.Matrix([
        [p.coeff_monomial(m) for p in polynomials] for m in monomials
    ])
    rows = []
    independent = sp.zeros(0, len(COORDINATES))
    for idx in range(coefficient_matrix.rows):
        candidate = independent.col_join(coefficient_matrix[idx, :])
        if candidate.rank() > len(rows):
            rows.append(idx)
            independent = candidate
            if len(rows) == len(COORDINATES):
                break
    if len(rows) != len(COORDINATES) or independent.det() == 0:
        raise ValueError('Declared invariant basis is not rank 10')
    inverse = independent.inv()
    hessians = [sp.hessian(v, phi) for v in basis]
    answer = {}
    for a, b in combinations_with_replacement(range(len(COORDINATES)), 2):
        # Coefficients for x_a*x_b, including the factor 1/2 for squares.
        loop = sp.Rational(1, 2) if a == b else sp.Integer(1)
        loop *= sum(hessians[a][i, j]*hessians[b][j, i]
                    for i in range(12) for j in range(12))
        loop = sp.expand(loop)
        coeff = sp.Poly(loop, *phi)
        projected = inverse * sp.Matrix([
            coeff.coeff_monomial(monomials[idx]) for idx in rows
        ])
        # The full exact polynomial closure is mandatory.
        rem = sp.Poly(sp.expand(loop - sum(projected[k]*basis[k]
                          for k in range(len(COORDINATES)))), *phi)
        if not rem.is_zero:
            raise ValueError(f'Scalar quartic beta leaves invariant basis: {a},{b}')
        result = {}
        for key, value in zip(COORDINATES, projected):
            if value != 0:
                exact = sp.simplify(value)
                if exact.is_Rational is not True:
                    raise ValueError(f'Unexpected nonrational beta {a},{b}->{key}: {exact}')
                result[key] = exact
        answer[(COORDINATES[a], COORDINATES[b])] = result
    return answer


def monomial_json(monomials: dict) -> list[dict]:
    return [
        {'couplings': list(pair), 'beta_coefficients': {
            key: str(val) for key, val in row.items()}}
        for pair, row in monomials.items()
    ]


def evaluated_beta(monomials: dict, point: dict) -> dict:
    """Evaluate quadratic beta in arbitrary real coordinates exactly or float."""
    if set(point) != set(COORDINATES):
        raise ValueError('Missing/extra quartic coordinate')
    result = {key: sp.Integer(0) for key in COORDINATES}
    for (a, b), row in monomials.items():
        for key, coefficient in row.items():
            result[key] += coefficient * point[a] * point[b]
    return result


def restricted_coefficients(monomials: dict) -> dict:
    """Five disputed channels, exact reference coefficients, no RGBeta use."""
    selected = (
        ('lambdaS2Adj', 'lambdaS2Adj'),
        ('lambdaS2Adj', 'lambdaH2Adj'),
        ('lambdaH2Adj', 'lambdaH2Adj'),
        ('lambdaT3.real', 'lambdaT3.real'),
        ('lambdaT3.imag', 'lambdaT3.imag'),
    )
    return {pair: monomials[pair] for pair in selected}
