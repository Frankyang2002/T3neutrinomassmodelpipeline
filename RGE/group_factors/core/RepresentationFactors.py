from __future__ import annotations

"""Exact representation-theory utilities shared by T3 group-factor modules.

This module contains only representation identities that are independent of
the particular beta function or matching calculation.
"""

from dataclasses import dataclass
from fractions import Fraction


def su2_quadratic_casimir_from_dimension(dimension: int) -> Fraction:
    """Return C2(j)=j(j+1) for the SU(2) irrep of dimension d=2j+1."""
    dimension = int(dimension)

    if dimension < 1:
        raise ValueError("SU(2) representation dimension must be positive.")

    return Fraction(dimension * dimension - 1, 4)


@dataclass(frozen=True)
class CanonicalYukawaLegFactors:
    dF: int
    dS: int
    d_larger: int
    G_scalar: Fraction
    G_lepton: Fraction
    G_heavy: Fraction
    self_coefficient: Fraction


def canonical_yukawa_leg_factors(
    dF: int,
    dS: int,
) -> CanonicalYukawaLegFactors:
    """Return canonical exact SU(2) leg factors for one T3 Yukawa invariant."""
    dF = int(dF)
    dS = int(dS)
    if dF < 1 or dS < 1:
        raise ValueError("SU(2) representation dimensions must be positive.")
    if abs(dS - dF) != 1:
        raise ValueError(
            f"T3 Yukawa requires dS=dF+/-1, got dF={dF}, dS={dS}."
        )
    d_larger = max(dF, dS)
    g_scalar = Fraction(d_larger, dS)
    g_lepton = Fraction(d_larger, 2)
    g_heavy = Fraction(d_larger, dF)
    return CanonicalYukawaLegFactors(
        dF=dF,
        dS=dS,
        d_larger=d_larger,
        G_scalar=g_scalar,
        G_lepton=g_lepton,
        G_heavy=g_heavy,
        self_coefficient=Fraction(1, 2) * (g_lepton + g_heavy),
    )
