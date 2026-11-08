"""Independent gaugeless one-loop psi^2 phi^2 reference from Fonseca et al.

Reference: R. M. Fonseca, P. Olgoso, J. Santiago, JHEP 07 (2026) 094,
  arXiv:2512.15866v1, eqs. (5), (6), (17), (18), (22).

Assumptions: only a^(5)_{psi phi^2}, Yukawas and scalar quartics are nonzero;
  all gauge couplings are zero. No dimension-5 dipoles, bosonic dimension-5
  operators, masses, dimension-6 insertions or trilinear sources are included.

The reference uses L5 = (1/4) a^(5)_ijab psi_i^T C psi_j phi_a phi_b + hc,
whereas Batch34 uses L5 = C_ijab psi_i psi_j phi_a phi_b + hc.
Thus a^(5) = 4*C and da^(5)/dt = 4*dC/dt. Since the retained beta is
real-linear in a^(5), this module can evaluate eq.(22) using C directly.

This is intentionally NOT implemented by calling the repository master beta.
The contraction follows the eight separate terms explicitly printed in eq.(22)
with independent external swaps (i<->j), (a<->b). It is a numerical
cross-check, NOT an alternative gauge-complete anomalous dimension.
"""
from __future__ import annotations

from itertools import permutations
from typing import Mapping

import numpy as np

REFERENCES = (
    'Fonseca-Olgoso-Santiago arXiv:2512.15866v1 Eq.(22) JHEP07(2026)094',
    'Aebischer-Bresciani-Selimovic arXiv:2512.16890v1 Eq.(4.85)',
)
CHANNELS = (
    'scalar_anomalous_dimension', 'yukawa_wavefunction',
    'fermion_anomalous_dimension', 'mixed_yukawa_gauge',
    'scalar_pair', 'conjugate_coefficient_yukawa',
)


def symmetrize_permutation_sum(raw: np.ndarray) -> np.ndarray:
    """Sum four independent external permutations, with repetitions on diagonals."""
    return (raw + raw.swapaxes(0,1) + raw.swapaxes(2,3)
            + raw.swapaxes(0,1).swapaxes(2,3))


def equation22_gaugeless(c: np.ndarray, yukawa: np.ndarray,
                          quartic: np.ndarray) -> dict[str,np.ndarray]:
    """Eight terms of the 2026 independent published formula, grouped as six.

    Inputs: full complex separately pair-symmetric C[i,j,a,b], symmetric
    complex y[i,j,a], fully S4-symmetric real lambda[a,b,c,d].
    Returns 16*pi^2 dC/d(ln mu), independently for each source channel.
    The output only covers the g=0, other-dim5=0 limit.
    """
    c=np.asarray(c,dtype=np.complex128)
    y=np.asarray(yukawa,dtype=np.complex128)
    lam=np.asarray(quartic,dtype=np.float64)
    if c.ndim!=4 or c.shape[0]!=c.shape[1] or c.shape[2]!=c.shape[3]:
        raise ValueError('Expected C[i,j,a,b] with equal pair dimensions')
    nf, _, ns, _=c.shape
    if y.shape!=(nf,nf,ns) or lam.shape!=(ns,)*4:
        raise ValueError('Yukawa or quartic tensor dimension mismatch')
    if not np.isfinite(c).all() or not np.isfinite(y).all() or not np.isfinite(lam).all():
        raise ValueError('Nonfinite reference input')
    if not np.allclose(c,c.swapaxes(0,1),atol=1e-13,rtol=1e-12) or not np.allclose(c,c.swapaxes(2,3),atol=1e-13,rtol=1e-12):
        raise ValueError('Wilson must have separate fermion/scalar exchange symmetries')
    if not np.allclose(y,y.swapaxes(0,1),atol=1e-13,rtol=1e-12):
        raise ValueError('Yukawa must be fermion symmetric')
    if any(not np.allclose(lam,lam.transpose(p),atol=1e-13,rtol=1e-12)
           for p in permutations(range(4))):
        raise ValueError('Scalar quartic must be fully S4 symmetric')

    raw={label:np.zeros_like(c) for label in CHANNELS}
    for a in range(ns):
        Ya=y[:,:,a]
        for b in range(ns):
            Yb=y[:,:,b]
            Cab=c[:,:,a,b]
            Cab_dag=Cab.conj().T
            # lambda term appears ONCE under the sum over the pair permutations,
            # with no sum over the dummy Yukawa scalar index d.
            raw['scalar_pair'][:,:,a,b] = .25 * np.einsum(
                'de,ijde->ij',lam[a,b],c,optimize=True)
            for d in range(ns):
                Yd=y[:,:,d]
                Yd_dag=Yd.conj().T
                Cad=c[:,:,a,d]
                Cbd=c[:,:,b,d]
                raw['scalar_anomalous_dimension'][:,:,a,b] += (
                    .25*Cbd*np.trace(Ya@Yd_dag+Yd@Ya.conj().T))
                raw['yukawa_wavefunction'][:,:,a,b] += (
                    .25*np.trace(Cab@Yd_dag)*Yd)
                raw['fermion_anomalous_dimension'][:,:,a,b] += (
                    .25*Cab@Yd_dag@Yd)
                raw['mixed_yukawa_gauge'][:,:,a,b] += (
                    Cad@Yd_dag@Yb + 2*Cad@Yb.conj().T@Yd)
                raw['conjugate_coefficient_yukawa'][:,:,a,b] += (
                    .25*np.trace(Cab_dag@Yd)*Yd
                    +.5*Yd@Cab_dag@Yd)
    return {name:symmetrize_permutation_sum(value)
            for name,value in raw.items()}


def sum_reference(channels: Mapping[str,np.ndarray]) -> np.ndarray:
    if set(channels)!=set(CHANNELS):
        raise ValueError('Expected exactly the six non-crossed reference channels')
    return sum((channels[key] for key in CHANNELS), np.zeros_like(channels[CHANNELS[0]]))
