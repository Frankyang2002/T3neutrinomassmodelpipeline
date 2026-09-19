from __future__ import annotations

"""Shared exact representation factors for T3 Yukawa vertices.

For a T3 Yukawa invariant involving the lepton doublet L, a heavy fermion F,
and a scalar S, SU(2) invariance requires dS = dF +/- 1. In the canonical
orthonormal Clebsch--Gordan convention define

    d_> = max(dF, dS),
    G_S = d_> / dS,
    G_L = d_> / 2,
    G_F = d_> / dF.

This module owns only these representation-theory identities and the SU(2)
quadratic Casimir. It does not contain beta-function normalisations or
branch-specific physics.
"""

from dataclasses import dataclass
from fractions import Fraction

from RGE.group_factors.RepresentationFactors import (
    su2_quadratic_casimir_from_dimension,
)


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
    """Return the canonical exact SU(2) leg factors for one T3 Yukawa."""
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


