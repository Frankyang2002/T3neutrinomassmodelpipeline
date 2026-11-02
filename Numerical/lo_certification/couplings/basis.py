"""Exact coordinate changes for real scalar components, independent of T3 models.

Conventions: s_original = R_F s_F = R_S s_S.  Thus the muS Hessian
in the fixed muF Wilson basis is R_F.T M_original(muS) R_F.  The
instantaneous eigenbasis change is s_F = (R_F.T R_S) s_S.  Nothing in
this module assumes equality of eigenvectors at muF and muS.
"""
from __future__ import annotations
import numpy as np


def orthogonal_matrix(value, n: int, label: str = 'rotation') -> np.ndarray:
    r = np.asarray(value)
    if (not np.isrealobj(r) or r.shape != (n,n) or
            not np.isfinite(r).all() or
            not np.allclose(r.T @ r, np.eye(n), rtol=1e-10, atol=1e-10)):
        raise ValueError(f'{label}: expected finite orthogonal real {n}x{n} matrix')
    return np.asarray(r, dtype=float)


def mass2_matrix(value, n: int, label: str = 'mass2') -> np.ndarray:
    a = np.asarray(value)
    if (not np.isrealobj(a) or a.shape != (n,n) or
            not np.isfinite(a).all() or
            not np.allclose(a, a.T, rtol=1e-10, atol=1e-8)):
        raise ValueError(f'{label}: real symmetric finite mass matrix required')
    a = np.asarray(a, float)
    if np.linalg.eigvalsh(a)[0] <= 0:
        raise ValueError(f'{label}: nonpositive scalar eigenvalue; no clipping')
    return a


def transport_mass2(m_original, r_original_from_target):
    n = len(m_original)
    r = orthogonal_matrix(r_original_from_target, n)
    m = mass2_matrix(m_original, n, 'original mass2')
    return mass2_matrix(r.T @ m @ r, n, 'transported mass2')


def transport_hhss(channels_original, r_original_from_target):
    """q_target[A,B,a,b] = R[i,a] R[j,b] q_original[A,B,i,j]."""
    q = np.asarray(channels_original)
    if q.ndim not in (4,5):
        raise ValueError('HHSS tensor must have rank 4 or channel-resolved rank 5')
    n = q.shape[-1]
    r = orthogonal_matrix(r_original_from_target, n)
    if q.shape[-4:] != (4,4,n,n) or not np.isrealobj(q) or not np.isfinite(q).all():
        raise ValueError('Invalid real HHSS channel shape/values')
    if not np.allclose(q, q.swapaxes(-1,-2), atol=1e-10, rtol=1e-10):
        raise ValueError('Scalar legs of HHSS are not symmetric')
    if not np.allclose(q, q.swapaxes(-3,-4), atol=1e-10, rtol=1e-10):
        raise ValueError('Higgs legs of HHSS are not symmetric')
    out = np.einsum('...ABij,ia,jb->...ABab', q, r, r, optimize=True)
    if not np.allclose(out, out.swapaxes(-1,-2), atol=1e-10, rtol=1e-10):
        raise ValueError('HHSS rotation corrupted scalar symmetry')
    return out


def transport_llss(c_original, r_original_from_target):
    c = np.asarray(c_original)
    n = c.shape[-1]
    r = orthogonal_matrix(r_original_from_target, n)
    if c.ndim < 2 or c.shape[-2] != n or not np.isfinite(c).all():
        raise ValueError('Invalid LLSS tensor')
    return np.einsum('...ij,ia,jb->...ab', c, r, r, optimize=True)


def relative_basis_rotation(r_original_from_F, r_original_from_S):
    n = len(r_original_from_F)
    rf = orthogonal_matrix(r_original_from_F, n, 'muF basis')
    rs = orthogonal_matrix(r_original_from_S, n, 'muS basis')
    return orthogonal_matrix(rf.T @ rs, n, 'relative muF/muS rotation')


def self_bilinear_hessian(dim: int, mass2: float, mu: complex) -> np.ndarray:
    """Matchete-native fixed Minus/Plus basis: V=m2 S†S-[mu(S†)^2+h.c.].

    This is the *source-convention* Hessian, not a universal BSM potential.
    """
    if dim not in (1,3):
        raise ValueError('Only source-native self-bilinear dimensions 1,3 are supported')
    z = complex(mu)
    if not np.isfinite([mass2, z.real, z.imag]).all():
        raise ValueError('Nonfinite self-bilinear parameters')
    i = np.eye(dim)
    return np.block([[(mass2-2*z.real)*i, -2*z.imag*i],
                     [-2*z.imag*i, (mass2+2*z.real)*i]])


def checked_eigenbasis(m):
    """Diagonalization for diagnostics; degenerate eigenvectors are not canonical."""
    m = mass2_matrix(m, len(m))
    values, u = np.linalg.eigh(m)
    if not np.allclose(u.T @ m @ u, np.diag(values), rtol=1e-10, atol=1e-7):
        raise ValueError('Eigenbasis validation failure')
    return values, u
