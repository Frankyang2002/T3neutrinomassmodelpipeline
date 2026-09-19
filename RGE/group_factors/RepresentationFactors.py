from __future__ import annotations

"""Exact representation-theory utilities shared by T3 group-factor modules.

This module contains only representation identities that are independent of
the particular beta function or matching calculation.
"""

from fractions import Fraction


def su2_quadratic_casimir_from_dimension(dimension: int) -> Fraction:
    """Return C2(j)=j(j+1) for the SU(2) irrep of dimension d=2j+1."""
    dimension = int(dimension)

    if dimension < 1:
        raise ValueError("SU(2) representation dimension must be positive.")

    return Fraction(dimension * dimension - 1, 4)
