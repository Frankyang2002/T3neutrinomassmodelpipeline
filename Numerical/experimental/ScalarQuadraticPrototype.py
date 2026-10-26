"""Isolated analytic/numerical prototypes for extra T3 scalar mass terms.

NO integration with the RG flow, Matchete propagators, matching, or thresholds.
Use this module only to explore whether quadratic mixing is numerically
relevant before deciding to refactor the production mass/threshold pipeline.

Conventions for ordinary T3-B/C at alpha=-1:

  V2 = m1sq S1^dag S1 + m2sq S2^dag S2
       + (mu12sq S1^T eps S2 + h.c.)
  eps = [[0,1],[-1,0]],  tilde(S1) = eps S1*
  Phi = (tilde(S1), S2)
  V2 = Phi^dag [[m1sq, -mu12sq],[-mu12sq*,m2sq]] Phi

Eigenvalues are the gauge-symmetric doublet masses-squared (before EWSB,
quartic-induced splitting, and scale-dependent running). A unitary U satisfies
U^dag M2 U = diag(eigenvalues) and Phi = U Phi_mass.

For a complex Y=0 singlet/triplet, with an SU(2)-invariant symmetric real
contraction, V2 = mSq S^dag S + 1/2(muSq S.S + h.c.) splits the two real
multiplets into mass-squared mSq +/- abs(muSq). The triplet components remain
threefold degenerate in the unbroken gauge phase.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class MixedDoubletSpectrum:
    mass_squared: tuple[float, float]
    mass: tuple[float, float]
    # Column eigenvectors of the mass matrix; U^dag M U = diag(m_-,m_+).
    unitary: np.ndarray
    mass_matrix: np.ndarray


def mixed_t3_bc_doublet_spectrum(
    m_s1_sq: float,
    m_s2_sq: float,
    mu12_sq: complex,
    *,
    require_positive: bool = True,
) -> MixedDoubletSpectrum:
    """Diagonalise the gauge-invariant B/C alpha=-1 quadratic doublet sector.

    Inputs are mass-squared parameters in GeV^2; complex mu12_sq includes its
    phase. Reject tachyonic eigenvalues by default because the inert-vacuum
    benchmark assumes no scalar field develops a VEV.
    """
    m1 = float(m_s1_sq)
    m2 = float(m_s2_sq)
    mu = complex(mu12_sq)
    if not (np.isfinite(m1) and np.isfinite(m2) and
            np.isfinite(mu.real) and np.isfinite(mu.imag)):
        raise ValueError("All quadratic mass-squared parameters must be finite")
    matrix = np.array([[m1, -mu], [-mu.conjugate(), m2]], dtype=np.complex128)
    eigenvalues, unitary = np.linalg.eigh(matrix)
    if require_positive and np.any(eigenvalues <= 0):
        raise ValueError("Quadratic mass matrix is not strictly positive definite")
    masses = tuple(float(np.sqrt(v)) if v >= 0 else float("nan")
                   for v in eigenvalues)
    return MixedDoubletSpectrum(
        mass_squared=tuple(float(v) for v in eigenvalues),
        mass=masses,
        unitary=unitary,
        mass_matrix=matrix,
    )


@dataclass(frozen=True)
class NeutralComplexScalarSplitting:
    su2_dimension: int
    mass_squared: tuple[float, float]
    mass: tuple[float, float]
    real_multiplet_degeneracy: tuple[int, int]


def neutral_complex_scalar_splitting(
    m_s_sq: float,
    mu_sq: complex,
    *,
    su2_dimension: int,
    require_positive: bool = True,
) -> NeutralComplexScalarSplitting:
    """Split complex Y=0 singlet or triplet into two real SU(2) multiplets.

    Applies to complex scalars in real representations only, d=1 or d=3.
    It does NOT represent the complete post-EWSB component spectrum.
    """
    if su2_dimension not in (1, 3):
        raise ValueError("Prototype supports only Y=0 singlet/triplet (d=1 or 3)")
    m = float(m_s_sq)
    mu = complex(mu_sq)
    if not (np.isfinite(m) and np.isfinite(mu.real) and np.isfinite(mu.imag)):
        raise ValueError("Mass-squared parameters must be finite")
    eigenvalues = (m - abs(mu), m + abs(mu))
    if require_positive and eigenvalues[0] <= 0:
        raise ValueError("Complex neutral scalar has a nonpositive quadratic eigenvalue")
    masses = tuple(float(np.sqrt(v)) if v >= 0 else float("nan")
                   for v in eigenvalues)
    return NeutralComplexScalarSplitting(
        su2_dimension=su2_dimension,
        mass_squared=eigenvalues,
        mass=masses,
        real_multiplet_degeneracy=(su2_dimension, su2_dimension),
    )
