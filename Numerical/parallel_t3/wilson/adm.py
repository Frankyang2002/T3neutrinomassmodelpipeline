"""Matrix-free dimension-five psi^2 phi^2 candidate anomalous dimension.

First preference: execute the EXISTING repo engines GeneralEFT1SevenTermFlow,
GeneralEFT1CrossedWardCandidate or GeneralEFT1SixTermReferenceCandidate.
A standalone algebraically identical fallback is provided for installation
checks outside the repository. No choice is presently physically certified.
"""
from __future__ import annotations
import numpy as np
from .coordinates import symmetry_residual

SYMMETRIC_TERMS=("yukawa_wavefunction","scalar_pair","mixed_yukawa_gauge",
                 "crossed_yukawa","conjugate_coefficient_yukawa",
                 "scalar_anomalous_dimension","fermion_anomalous_dimension")
VARIANTS=('fonseca_eq22_six_term','ward_covariant_candidate','source_v1_unstarred')


def gauge_ward_residual(c:np.ndarray, t:np.ndarray, theta:np.ndarray) -> float:
    var=(np.einsum('ki,kjab->ijab',t,c,optimize=True)
         +np.einsum('kj,ikab->ijab',t,c,optimize=True)
         +np.einsum('ca,ijcb->ijab',theta,c,optimize=True)
         +np.einsum('cb,ijac->ijab',theta,c,optimize=True))
    return float(np.max(np.abs(var))/max(1e-300,np.max(np.abs(c))))


def _fallback(c,y,lam,gauge,variant):
    e=lambda sig,*x:np.einsum(sig,*x,optimize=True)
    z=c.conj(); terms={}
    terms['yukawa_wavefunction']=e('ijd,dab->ijab',y,e('kld,klab->dab',y.conj(),c))
    pair=e('abcd,ijcd->ijab',lam,c)
    for g,gens in gauge:
        for tf,ts in gens: pair-=2*g*g*e('ac,bd,ijcd->ijab',ts,ts,c)
    terms['scalar_pair']=(pair+pair.swapaxes(2,3))/2
    tmp=(2*e('jkc,klb,ilac->ijab',y,y.conj(),c)
         +e('jkb,klc,ilac->ijab',y,y.conj(),c))
    for g,gens in gauge:
        for tf,ts in gens:
            tmp+=4*g*g*e('lj,bc,ilac->ijab',tf,ts,c)
    terms['mixed_yukawa_gauge']=(tmp+tmp.swapaxes(0,1)+tmp.swapaxes(2,3)+tmp.transpose(1,0,3,2))
    if variant=='fonseca_eq22_six_term':
        cross=np.zeros_like(c)
    else:
        cross=e('jla,ikc,lkbc->ijab',y,y,z if variant=='ward_covariant_candidate' else c)
        cross=(cross+cross.swapaxes(0,1)+cross.swapaxes(2,3)+cross.transpose(1,0,3,2))
    terms['crossed_yukawa']=cross
    antitmp=e('ild,klab->ikabd',y,z)
    terms['conjugate_coefficient_yukawa']=(2*e('jkd,ikabd->ijab',y,antitmp)
                                          +e('ijd,kld,klab->ijab',y,y,z))
    gs=0.5*(e('ija,ijb->ab',y,y.conj())+e('ijb,ija->ab',y,y.conj()))
    gf=0.5*e('ika,jka->ij',y,y.conj())
    for g,gens in gauge:
        for tf,ts in gens:
            gs-=4*g*g*(ts@ts)
            cas=tf@tf
            gf-=3*g*g*(cas.T if variant!='source_v1_unstarred' else cas)
    x=e('bd,ijad->ijab',gs,c)
    terms['scalar_anomalous_dimension']=x+x.swapaxes(2,3)
    x=e('jk,ikab->ijab',gf,c)
    terms['fermion_anomalous_dimension']=x+x.swapaxes(0,1)
    return terms


def beta(c:np.ndarray, y:np.ndarray, lam:np.ndarray,
         gauge:tuple=(), *, variant:str='fonseca_eq22_six_term',
         prefer_repo:bool=True) -> tuple[dict[str,np.ndarray],str]:
    """Return (seven named channels, execution_backend), each 16*pi^2 dC/dlnmu.

    'source_v1_unstarred' is known to fail covariance in general bases.
    'ward_covariant_candidate' has a C* crossed channel, unverified loop factor.
    'fonseca_eq22_six_term' matches independent gaugeless Eq.(22), but full
    gauge completion has not been separately validated. NO variant is certified.
    """
    if variant not in VARIANTS: raise ValueError('Unsupported ADM candidate')
    c=np.asarray(c,np.complex128)
    y=np.asarray(y,np.complex128)
    lam=np.asarray(lam,float)
    nf,ns=y.shape[0],y.shape[2]
    if c.shape!=(nf,nf,ns,ns) or y.shape!=(nf,nf,ns) or lam.shape!=(ns,)*4:
        raise ValueError('ADM tensor shape mismatch')
    if not all(np.isfinite(a).all() for a in (c,y,lam)):
        raise ValueError('Nonfinite ADM inputs')
    if symmetry_residual(c)>max(1e-24,1e-10*np.max(np.abs(c))):
        raise ValueError('Wilson violates pair exchange symmetry')
    if not np.allclose(y,y.swapaxes(0,1),atol=1e-12):
        raise ValueError('Yukawa must be fermion-symmetric')
    if any(not np.allclose(lam,lam.swapaxes(0,k),atol=1e-12) for k in (1,2,3)):
        raise ValueError('Quartic must be fully symmetric')
    if prefer_repo:
        try:
            from RGE.running.intermediate.GeneralEFT1SevenTermFlow import Couplings,master_beta
            obj=Couplings(y,lam,tuple(gauge))
            if variant=='source_v1_unstarred':
                result=master_beta(c,obj)
            elif variant=='ward_covariant_candidate':
                from RGE.running.intermediate.GeneralEFT1CrossedWardCandidate import candidate_beta
                result=candidate_beta(c,obj)
            else:
                from RGE.running.intermediate.GeneralEFT1SixTermReferenceCandidate import six_term_beta
                result=six_term_beta(c,obj)
            tag='existing_repository_matrix_free_engine'
        except ImportError:
            result=_fallback(c,y,lam,gauge,variant)
            tag='portable_formula_mirror_no_repository_engine'
    else:
        result=_fallback(c,y,lam,gauge,variant)
        tag='portable_formula_mirror_explicit'
    if set(result)!=set(SYMMETRIC_TERMS): raise RuntimeError('Beta channel set changed')
    if not all(np.isfinite(a).all() for a in result.values()):
        raise FloatingPointError('Nonfinite Wilson beta')
    return result,tag


def total_beta(channels:dict) -> np.ndarray:
    return sum((channels[k] for k in SYMMETRIC_TERMS),
               np.zeros_like(channels[SYMMETRIC_TERMS[0]]))
