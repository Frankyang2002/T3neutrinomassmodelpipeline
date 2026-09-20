from __future__ import annotations

"""Shared tensor algebra for exact lambdaT3 scalar-quartic recouplings.

This module operates on the one-based ``SparseQuarticTensor`` convention used
by ``EFT1TensorAdapters``.  It contains only common basis extraction,
cross-contraction, inner-product, and projection machinery.  It does not define
which quartic couplings are physical targets or assign any recoupling
coefficient by hand.
"""

from collections import Counter
from itertools import combinations_with_replacement
from math import factorial

import sympy as sp

from RGE.running.eft1.EFT1TensorAdapters import SparseQuarticTensor


def basis_tensor(
    full: SparseQuarticTensor,
    coupling_name: str,
    *,
    identify_conjugate: bool,
) -> SparseQuarticTensor:
    """Extract the coefficient tensor of one coupling from the full tensor."""
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
    """Return the largest one-based real-scalar index present in the tensor."""
    return max((max(key) for key in full.entries), default=0)


def multiplicity_weight(key: tuple[int, int, int, int]) -> int:
    """Number of ordered index tuples represented by one sorted quartic key."""
    counts = Counter(key)
    weight = factorial(4)
    for count in counts.values():
        weight //= factorial(count)
    return weight


def tensor_inner_product(
    left: SparseQuarticTensor,
    right: SparseQuarticTensor,
) -> sp.Expr:
    """Hermitian tensor inner product summed over all ordered quartic indices."""
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
    """Coefficient of lambdaT3*X in one beta_abcd component."""
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
    """Project ``generated`` onto ``target`` and return coefficient/residual/norm."""
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
