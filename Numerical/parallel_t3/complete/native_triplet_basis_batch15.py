"""Source-exact Matchete T3MixCG triplet -> highest-weight spin-1 basis.

Matchete's T3BSMd3 source components in the verified A(-4)/D(-2) CG
are Cartesian-like.  The reference master ADM and RGBeta columns use
SU(2) spin-1 |m=+1,0,-1> components.  This module derives the unitary
change of basis from the *saved* CG, never from a guessed permutation.

The source tensor C and all Matchete source quartics must be transformed
with this map. RGBeta reference quartics already in highest-weight basis
must NOT be transformed again.  This SU2 basis result does not certify
Matchete-Weyl spinor normalization or full finite matching.
"""
from __future__ import annotations
import math
from typing import Any, Mapping
import numpy as np


def _require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError('Batch15 SU2 basis: '+message)


def fundamental() -> tuple[np.ndarray,...]:
    return (np.array([[0,.5],[.5,0]],complex),
            np.array([[0,-.5j],[.5j,0]],complex),
            np.diag([.5,-.5]).astype(complex))


def spherical_spin_one() -> tuple[np.ndarray,...]:
    tp=np.zeros((3,3),complex)
    tp[0,1]=tp[1,2]=math.sqrt(2)
    return ((tp+tp.conj().T)/2,(tp-tp.conj().T)/(2j),
            np.diag([1,0,-1]).astype(complex))


def triplet_cg_from_seed(seed: Mapping[str,Any]) -> np.ndarray:
    """Use the actual Matchete SparseArray parser; no hardcoded CG values."""
    from RGE.running.intermediate.ScalarOnlyTensorAdapters import load_wilson_cg_registry
    reg=load_wilson_cg_registry({'CGRegistry':seed.get('CGRegistry',[])})
    _require('T3MixCG' in reg,'T3MixCG not found in native source')
    obj=reg['T3MixCG']
    _require(tuple(obj.reps)==('SU2L[fund]','SU2L[fund]','T3BSMd3') or
             tuple(obj.reps)==('Bar[SU2L[fund]]','Bar[SU2L[fund]]','T3BSMd3'),
             'unexpected T3MixCG representation labels')
    cg=np.asarray(obj.matrix.tolist(),dtype=np.complex128)
    _require(cg.shape==(2,2,3) and np.isfinite(cg).all(),
             'expected finite rank-three 2x2x3 native CG')
    return cg


def generators_from_cg(cg: np.ndarray) -> tuple[np.ndarray,...]:
    """Solve (t⊗1+1⊗t) CG + CG.T_source = 0 for source triplet T.

    The source has two *dual* fundamental slots and one fundamental triplet
    slot. This equation follows from covariance of the exact invariant CG.
    """
    x=np.asarray(cg,complex)
    _require(x.shape==(2,2,3) and np.isfinite(x).all(),
             'native triplet CG must be complex 2x2x3')
    a=x.reshape(4,3)
    _require(np.linalg.matrix_rank(a)==3,'CG does not span a triplet')
    out=[]
    for t in fundamental():
        b=np.einsum('ki,kjm->ijm',t,x)+np.einsum('kj,ikm->ijm',t,x)
        trip=np.linalg.lstsq(a,-b.reshape(4,3),rcond=None)[0]
        err=np.linalg.norm(a@trip+b.reshape(4,3)) / max(np.linalg.norm(b),1e-100)
        _require(err<1e-10,'native CG fails triplet intertwining identity')
        _require(np.linalg.norm(trip-trip.conj().T)<1e-10,
                 'CG-derived triplet generator not Hermitian')
        out.append(trip)
    t1,t2,t3=out
    _require(np.linalg.norm(t1@t2-t2@t1-1j*t3)<1e-10,
             'native CG-derived generators do not close su(2)')
    _require(np.linalg.norm(t2@t3-t3@t2-1j*t1)<1e-10 and
             np.linalg.norm(t3@t1-t1@t3-1j*t2)<1e-10,
             'native triplet generator su(2) axes inconsistent')
    return tuple(out)


def source_from_spherical(cg: np.ndarray) -> tuple[np.ndarray,dict]:
    """Return U with z_source = U z_spherical, U^dagger T_source U = T_spin1.

    Fix the otherwise arbitrary eigenvector phases by demanding positive
    spin-1 raising matrix elements for (+1,0) and (0,-1).
    """
    gens=generators_from_cg(cg)
    eig,vec=np.linalg.eigh(gens[2])
    order=np.argsort(eig)[::-1]
    eig=eig[order]; vec=vec[:,order]
    _require(np.max(abs(eig-np.array([1.,0.,-1.])))<1e-10,
             'native triplet has wrong T3 eigenvalues')
    tp=vec.conj().T@(gens[0]+1j*gens[1])@vec
    _require(abs(tp[0,1])>1e-12 and abs(tp[1,2])>1e-12,
             'triplet raising operator basis degenerate')
    p1=-np.angle(tp[0,1]);p2=p1-np.angle(tp[1,2])
    U=vec@np.diag(np.exp(1j*np.asarray([0.,p1,p2])))
    errs=[float(np.linalg.norm(U.conj().T @ a @ U - b))
          for a,b in zip(gens,spherical_spin_one())]
    _require(max(errs)<1e-10 and
             np.linalg.norm(U.conj().T@U-np.eye(3))<1e-10,
             'native/spherical SU2 unitary intertwiner failed')
    return U,{'cg_triplet_spherical_intertwining_residual_max':max(errs),
              'native_source_triplet_generator_type':'Cartesian-like Matchete T3BSMd3',
              'target_triplet_generator_type':'SU2 highest-weight spin-1 m=+1,0,-1',
              'source_CG_intertwiner_unique_up_to_global_complex_phase':True,
              'native_to_RGBeta_quartic_normalization_certified':False}


def complex_to_interleaved_real(U: np.ndarray) -> np.ndarray:
    """z_source=U z_spherical -> phi_source=R phi_spherical (interleaved R,I)."""
    z=np.asarray(U,complex)
    _require(z.shape==(3,3) and np.linalg.norm(z.conj().T@z-np.eye(3))<1e-10,
             'expected unitary 3x3 complex basis transform')
    R=np.zeros((6,6),float)
    for i in range(3):
        for j in range(3):
            a,b=z[i,j].real,z[i,j].imag
            R[2*i:2*i+2,2*j:2*j+2]=[[a,-b],[b,a]]
    _require(np.linalg.norm(R.T@R-np.eye(6))<1e-10,
             'native-to-spherical real transform not orthogonal')
    return R


def heavy_source_from_reference(model: str, seed: Mapping[str,Any]) -> tuple[np.ndarray,dict]:
    """Heavy source -> canonical spherical reference, in 8 real coordinates.

    A: (S1 singlet 2, S2 triplet 6).
    D: (S1 triplet 6, S2 complex singlet 2), before the D singlet split.
    """
    _require(model in ('A_-4','D_-2'),'only A(-4) and D(-2) supported')
    U,meta=source_from_spherical(triplet_cg_from_seed(seed))
    R=complex_to_interleaved_real(U)
    out=np.eye(8)
    offset=2 if model=='A_-4' else 0
    out[offset:offset+6,offset:offset+6]=R
    return out,meta


def relative_operator_ward(c:np.ndarray,t:np.ndarray,theta:np.ndarray) -> float:
    """Variation of C_ijab psi_i psi_j phi_a phi_b in the selected basis."""
    v=(np.einsum('ki,kjab->ijab',t,c,optimize=True)
       +np.einsum('kj,ikab->ijab',t,c,optimize=True)
       +np.einsum('ca,ijcb->ijab',theta,c,optimize=True)
       +np.einsum('cb,ijac->ijab',theta,c,optimize=True))
    return float(np.max(abs(v))/max(np.max(abs(c)),1e-100))
