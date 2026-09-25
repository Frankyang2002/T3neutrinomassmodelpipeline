from __future__ import annotations

"""
We get common operations on quartic tensors here, note that T3 is lambda5 scalar

"""

from collections import Counter
from itertools import combinations_with_replacement
from math import factorial

import sympy as sp

from RGE.running.intermediate.ScalarOnlyTensorAdapters import SparseQuarticTensor


def basis_tensor(
    full: SparseQuarticTensor,
    coupling_name: str,
    *,
    identify_conjugate: bool,
) -> SparseQuarticTensor:
    """Extract the coefficient tensor of one coupling from the full tensor.
    So if we have Lambda_abcd = Lambda_a T^a + Lambda_b T_b
    if we do basistensor(...,"lambda_a"), we get T^a"""
    symbol = sp.Symbol(coupling_name)
    conjugate = sp.conjugate(symbol)

    entries: dict[tuple[int, int, int, int], sp.Expr] = {}
    for key, raw in full.nonzero_items():
        expr = sp.expand(raw)
        if identify_conjugate:
            expr = sp.expand(expr.xreplace({conjugate: symbol}))
        coefficient = sp.simplify(expr.coeff(symbol))
        if coefficient != 0:
            entries[key] = coefficient

    return SparseQuarticTensor(entries)


def scalar_dimension(full: SparseQuarticTensor) -> int:
    """get number of real scalars in tensor"""
    return max((max(key) for key in full.entries), default=0)


def multiplicity_weight(key: tuple[int, int, int, int]) -> int:
    """Calculate distinct permutations, eg: 1234 = 4!, 1112 = 4!/3!"""
    counts = Counter(key)
    weight = factorial(4)
    for count in counts.values():
        weight //= factorial(count)
    return weight


def tensor_inner_product(
    left: SparseQuarticTensor,
    right: SparseQuarticTensor,
) -> sp.Expr:
    """Defines inner product, <A,B> = sum_abcd A*_abcd B_abcd."""
    total = sp.S.Zero
    for key in set(left.entries) | set(right.entries):
        total += (
            multiplicity_weight(key)
            * sp.conjugate(left[key])
            * right[key]
        )
    return sp.simplify(total)


def cross_component(
    t3: SparseQuarticTensor,
    other: SparseQuarticTensor,
    a: int,
    b: int,
    c: int,
    d: int,
    scalar_dimension_value: int,
) -> sp.Expr:
    """
    We sum over (ab)(cd), (ac)(bd), (ad)(bc), summing our internal scalar indices
    This gives us, lambdaT3*lambdaX(T^T3 T^X + T^X T^T3) for each of the pairing (..)(..)
    We get this from expanding lambda*lambda
    Then we can get our coefficient of it
    This is for the RGE of the quartics themselves"""
    total = sp.S.Zero

    pairings = (
        ((a, b), (c, d)),
        ((a, c), (b, d)),
        ((a, d), (b, c)),
    )

    for (p, q), (r, s) in pairings:
        for e in range(1, scalar_dimension_value + 1):
            for f in range(1, scalar_dimension_value + 1):
                total += (
                    t3[p, q, e, f] * other[e, f, r, s]
                    + other[p, q, e, f] * t3[e, f, r, s]
                )

    return sp.simplify(total)


def build_cross_tensor(
    t3: SparseQuarticTensor,
    other: SparseQuarticTensor,
    scalar_dimension_value: int,
) -> SparseQuarticTensor:
    """Build the complete symmetric lambdaT3*X scalar-loop tensor."""
    entries: dict[tuple[int, int, int, int], sp.Expr] = {}

    for key in combinations_with_replacement(
        range(1, scalar_dimension_value + 1),
        4,
    ):
        value = cross_component(
            t3,
            other,
            *key,
            scalar_dimension_value,
        )
        if value != 0:
            entries[key] = value

    return SparseQuarticTensor(entries)


def project_onto_direction(
    generated: SparseQuarticTensor,
    target: SparseQuarticTensor,
) -> tuple[sp.Expr, SparseQuarticTensor, sp.Expr]:
    """We project our tensor onto our target tensor
    c = <T,G>/<T,T>
    We get residuals R=G-cT
    Where if R=0, G=ct and is directly proportional
    If not, then our group factor isnt there, and so we just return the residue
    """
    norm = sp.simplify(tensor_inner_product(target, target))
    if norm == 0:
        raise ValueError("target quartic tensor has zero norm.")

    coefficient = sp.simplify(
        tensor_inner_product(target, generated) / norm
    )

    residual_entries: dict[tuple[int, int, int, int], sp.Expr] = {}
    for key in set(generated.entries) | set(target.entries):
        residual = sp.simplify(
            generated[key] - coefficient * target[key]
        )
        if residual != 0:
            residual_entries[key] = residual

    return coefficient, SparseQuarticTensor(residual_entries), norm
