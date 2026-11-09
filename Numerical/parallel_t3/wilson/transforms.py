"""Explicit field-coordinate transforms; NEVER imply their physical certification.

psi_native[p] = U[p,i] psi_canonical[i] (unitary),
phi_native[a] = R[a,A] phi_canonical[A] (real orthogonal).
In L5=C psi psi phi phi the tensors transform with U, U, R, R,
not conjugate U. In kinetic generators T_can=U^dag T_native U.
"""
from __future__ import annotations
import numpy as np


def check_maps(U,R,tol=1e-10):
    U=np.asarray(U,complex);R=np.asarray(R,float)
    if U.ndim!=2 or R.ndim!=2 or U.shape[0]!=U.shape[1] or R.shape[0]!=R.shape[1]:
        raise ValueError('Coordinate transforms must be square matrices')
    if not np.isfinite(U).all() or not np.isfinite(R).all():raise ValueError('Nonfinite transform')
    if not np.allclose(U.conj().T@U,np.eye(len(U)),atol=tol) or not np.allclose(R.T@R,np.eye(len(R)),atol=tol):
        raise ValueError('Require unitary Weyl and real-orthogonal scalar maps')
    return U,R


def wilson_to_canonical(C,U,R):
    U,R=check_maps(U,R)
    return np.einsum('pi,qj,ca,db,pqcd->ijab',U,U,R,R,np.asarray(C,complex),optimize=True)


def yukawa_to_canonical(Y,U,R):
    U,R=check_maps(U,R)
    return np.einsum('pi,qj,ca,pqc->ija',U,U,R,np.asarray(Y,complex),optimize=True)


def quartic_to_canonical(lam,R):
    n=len(R);_,R=check_maps(np.eye(1),R)
    out=np.asarray(lam,float)
    if out.shape!=(n,)*4:raise ValueError('lambda4 axes do not match scalar transform')
    for axis in range(4):
        out=np.tensordot(out,R,axes=([axis],[0]))
        out=np.moveaxis(out,-1,axis)
    return out


def scalar_mass_to_canonical(m2,R):
    _,R=check_maps(np.eye(1),R)
    M=np.asarray(m2,complex)
    if M.shape!=(len(R),len(R)):raise ValueError('Scalar mass mismatch')
    return R.T@M@R


def fermion_mass_to_canonical(m,U):
    U,_=check_maps(U,np.eye(1));M=np.asarray(m,complex)
    if M.shape!=U.shape:raise ValueError('Weyl mass mismatch')
    return U.T@M@U


def generator_residual(native_f,native_s,canonical_f,canonical_s,U,R):
    U,R=check_maps(U,R)
    return max(float(np.max(np.abs(U.conj().T@native_f@U-canonical_f))),
               float(np.max(np.abs(R.T@native_s@R-canonical_s))))
