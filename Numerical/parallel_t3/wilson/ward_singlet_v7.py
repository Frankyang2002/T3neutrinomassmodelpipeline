"""A(-4) Ward-singlet *diagnostic* boundary, never native matching.

The result is the orthogonal Frobenius projection of the already-supplied
LLSS block onto the canonical SU(2)xU(1) invariant template span computed by
`synthetic_invariants.synthetic_llss`. No Matchete <-> RGBeta map is inferred.
This is a controlled change of the numerical EFT1 boundary, NOT a repair of
the original UV fermion matching or a physical matching certificate.

Fermion-first EFT ordering remains UV -> all F at MF -> both S at MS -> SMEFT.
"""
from __future__ import annotations

import hashlib
from pathlib import Path
import numpy as np

from .basis import build_basis
from .coordinates import symmetry_residual
from .synthetic_invariants import synthetic_llss

STATUS='WardProjectedBoundaryDiagnosticOnly'


def sha256(path: Path) -> str:
    hasher=hashlib.sha256()
    with Path(path).open('rb') as fh:
        for block in iter(lambda:fh.read(1024*1024),b''):
            hasher.update(block)
    return hasher.hexdigest()


def _gauge_norm(c:np.ndarray,basis) -> dict[str,float]:
    """Max Frobenius Ward numerator / Frobenius tensor norm per group."""
    denom=float(np.linalg.norm(c))
    if denom<=1e-100:
        return {key:0. for key in ('SU2','U1','SU3')}
    out={}
    for name,gens in basis.generators().items():
        val=0.
        for tf,ts in gens:
            ward=(np.einsum('ki,kjab->ijab',tf,c,optimize=True)
                 +np.einsum('kj,ikab->ijab',tf,c,optimize=True)
                 +np.einsum('ca,ijcb->ijab',ts,c,optimize=True)
                 +np.einsum('cb,ijac->ijab',ts,c,optimize=True))
            val=max(val,float(np.linalg.norm(ward))/denom)
        out[name]=val
    return out


def project_a4_boundary(initial:np.ndarray, *, model:str='A_-4',
                        min_kept_fraction:float=1e-8) -> tuple[np.ndarray,dict]:
    """Orthogonal projection ONLY of the LLSS block of A(-4).

    Invariant template basis: six symmetric flavour combinations, each with
    unique SU(2) coupling for the A singlet-triplet scalar representations.
    Projection coefficients may depend on the unknown native basis; results
    are explicitly diagnostic and must never be counted as genuine SMEFT.
    """
    if model!='A_-4':
        raise ValueError('v7 projection authorised only for A_-4; other source maps unresolved')
    basis=build_basis(model)
    c=np.asarray(initial,np.complex128)
    shape=(basis.nf,basis.nf,basis.ns,basis.ns)
    if c.shape!=shape or not np.isfinite(c).all():
        raise ValueError('Bad canonical full Wilson boundary geometry/nonfinite input')
    if symmetry_residual(c)>max(1e-24,1e-9*float(np.max(np.abs(c)))):
        raise ValueError('Input lacks declared fermion/scalar exchange symmetry')
    ls=[i for pair in basis.leptons for i in pair]
    hs=list(range(4,basis.ns))
    ix=np.ix_(ls,ls,hs,hs)
    source=np.asarray(c[ix],complex)
    norm=float(np.linalg.norm(source))
    if norm<1e-100:
        raise ValueError('Zero LLSS source: cannot infer nontrivial Ward projection')
    # Use the established synthetic invariant *as an algebraic basis*, not as
    # a replacement source or as a physical coupling normalization.
    modes=[]
    for p in range(3):
        for q in range(p,3):
            flavour=np.zeros((3,3),complex)
            flavour[p,q]=flavour[q,p]=1.
            t=synthetic_llss(basis,flavour)
            modes.append((p,q,np.asarray(t[ix],complex)))
    replacement=np.zeros_like(source)
    components=[]
    for p,q,t in modes:
        tnorm=float(np.vdot(t,t).real)
        if tnorm<=1e-80: raise RuntimeError('Invariant template has zero norm')
        coef=np.vdot(t,source)/tnorm
        replacement+=coef*t
        components.append({'flavour_pair':[p+1,q+1],
                           'template_coefficient':{'real':float(coef.real),'imag':float(coef.imag)},
                           'contribution_norm':float(abs(coef)*np.sqrt(tnorm))})
    kept=float(np.linalg.norm(replacement))
    removed=float(np.linalg.norm(source-replacement))
    fraction=kept/norm
    if fraction<min_kept_fraction:
        raise ValueError(f'Ward invariant projection nearly zero: kept_fraction={fraction:.4g}')
    result=c.copy();result[ix]=replacement
    # Only the original LLSS block is replaced; all other fermion and
    # scalar sectors, including independently present hard C5, are unchanged.
    outside=c.copy();outside[ix]=0
    other_norm=float(np.linalg.norm(outside))
    restored_ward=_gauge_norm(result,basis)
    if max(restored_ward.values())>1e-9:
        raise ValueError('Projected boundary still violates a gauge Ward identity; '
                         f'other_blocks_norm={other_norm:.6g}; ward={restored_ward}')
    proj_check=0.
    for _,_,t in modes:
        proj_check=max(proj_check,float(abs(np.vdot(t,source-replacement))
                                     /max(np.linalg.norm(t)*norm,1e-100)))
    if proj_check>1e-10:raise RuntimeError('Projection is not orthogonal to template span')
    after=result.copy();after[ix]=0
    outside_unchanged=bool(np.array_equal(outside,after))
    return result,{
        'format':'t3-ward-projected-A4-boundary-diagnostic-v7',
        'status':STATUS,'model':model,
        'unverified_Matchete_RGBeta_native_scalar_basis_map':True,
        'synthetic_invariant_used_only_as_normalised_linear_basis':True,
        'original_full_tensor_norm':float(np.linalg.norm(c)),
        'original_LLSS_norm':norm,'projected_LLSS_norm':kept,
        'discarded_LLSS_norm':removed,
        'kept_LLSS_frobenius_fraction':fraction,
        'discarded_LLSS_frobenius_fraction':removed/norm,
        'original_Ward':_gauge_norm(c,basis),
        'projected_Ward':restored_ward,
        'orthogonal_projection_residual':proj_check,
        'non_LLSS_tensor_unchanged':outside_unchanged,
        'projection_components':components,
        'not_physical_matching':True,
        'all_fermions_integrated_at_MF':True,
        'both_scalars_remain_active_to_MS':True,
    }
