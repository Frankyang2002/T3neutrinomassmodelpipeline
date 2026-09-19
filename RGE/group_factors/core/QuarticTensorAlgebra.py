from __future__ import annotations

"""Shared sparse quartic-tensor bookkeeping for portal recoupling derivations.

The tensors represented here use sorted zero-based real-scalar index tuples and
the fully symmetric convention

    V4 = (1/4!) lambda_abcd phi_a phi_b phi_c phi_d.

This module contains only generic tensor bookkeeping. It does not define any
portal operator basis, SU(2) generator, gauge normalization, or recoupling
coefficient.
"""

from collections import Counter
from math import factorial

import sympy as sp


def tensor_from_polynomial(expr, variables):
    """Convert a quartic polynomial into a sparse fully symmetric tensor."""
    poly = sp.Poly(sp.expand(expr), *variables)
    entries = {}

    for powers, coefficient in poly.terms():
        if sum(powers) != 4:
            continue

        key = []
        multiplicity = 1
        for index, power in enumerate(powers):
            key.extend([index] * power)
            multiplicity *= factorial(power)

        value = sp.simplify(coefficient * multiplicity)
        if value != 0:
            entries[tuple(sorted(key))] = value

    return entries


def tensor_get(tensor, a, b, c, d):
    """Read one component of a fully symmetric sparse quartic tensor."""
    return tensor.get(tuple(sorted((a, b, c, d))), sp.S.Zero)


def ordered_weight(key):
    """Number of ordered index tuples represented by one sorted quartic key."""
    counts = Counter(key)
    value = factorial(4)
    for count in counts.values():
        value //= factorial(count)
    return value


def tensor_inner(left, right):
    """Hermitian inner product over all ordered quartic-index tuples."""
    return sp.simplify(sum(
        ordered_weight(key)
        * sp.conjugate(left.get(key, 0))
        * right.get(key, 0)
        for key in set(left) | set(right)
    ))


def pair_maps(tensor, dimension):
    """Precompute two-index slices used in scalar-loop quartic contractions."""
    maps = {}

    for a in range(dimension):
        for b in range(a, dimension):
            vector = {}
            for e in range(dimension):
                for f in range(dimension):
                    value = tensor_get(tensor, a, b, e, f)
                    if value != 0:
                        vector[(e, f)] = value
            maps[(a, b)] = vector

    return maps


def sparse_dot(left, right):
    """Sparse bilinear contraction without complex conjugation."""
    if len(left) > len(right):
        left, right = right, left

    return sp.simplify(sum(
        value * right.get(key, 0)
        for key, value in left.items()
    ))


def restrict_sector(
    tensor,
    groups,
    required,
    *,
    reject_unassigned: bool = False,
):
    """Keep entries with the requested number of real legs in each field group."""
    out = {}

    for key, value in tensor.items():
        counts = {name: 0 for name in groups}
        valid = True

        for index in key:
            found = False
            for name, indices in groups.items():
                if index in indices:
                    counts[name] += 1
                    found = True
                    break

            if reject_unassigned and not found:
                valid = False
                break

        if valid and counts == required:
            out[key] = value

    return out
