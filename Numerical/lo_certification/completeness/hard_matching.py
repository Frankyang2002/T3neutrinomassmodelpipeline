"""All-F-at-one-scale hard LLHH matching: conditional scalar-quartic triangle.

Convention is that of Numerical/joint_scalar_matching/core.py (checked on
GitHub main 2c29e5ba, 2026-11-02):
    delta_C5_soft = -project[sum_ab C_LLss_ab lambda_HHss_ab B0bar]/(2*16pi^2).

The hard dim-5 local term for a heavy *mass eigenstate* k, with its separately
normalized tree LLss tensor C_k, is
    C5_hard(muF) = project[sum_k,ab C_k lambda (1+ln(muF^2/Mk^2))]
                   /(2*16pi^2).
This follows by subtracting the EFT1 soft massless-heavy-propagator region from
an explicit one-loop scalar-scalar-heavy-fermion triangle. It is only one
operator topology, not complete matching for the most general UV Lagrangian.

The full zero-external-momentum UV triangle, without a heavy-mass expansion,
is also evaluated. Its difference from the dim-5 EFT result is ONE LOOP but
higher order in mS^2/MF^2 (dimension 7 and above), not a two-loop correction.
No physical C5 is certified without an independently normalized per-k tree
boundary, full UV operator/vertex coverage and complete one-loop subtraction.
"""
from __future__ import annotations

from math import isfinite, log, log1p, pi
from typing import Any
import numpy as np

LOOP = 16.0 * pi**2


def b0bar(a: float, b: float, mu: float) -> float:
    """MSbar scalar B0(0;a,b;mu^2), mass squared a,b positive."""
    a, b, mu = float(a), float(b), float(mu)
    if not all(map(isfinite, (a, b, mu))) or min(a, b, mu) <= 0:
        raise ValueError('Positive finite scalar mass squares and mu required')
    mean = (a+b)/2
    d = (a-b)/(a+b)
    if abs(d) < 1e-3:
        d2 = d*d
        return -log(mean / (mu*mu)) + d2/6 + d2*d2/20 + d2**3/42
    return 1.0-(a*log(a/(mu*mu))-b*log(b/(mu*mu)))/(a-b)


def _f(r: float) -> float:
    x = r-1.0
    if abs(x) < 1e-5:
        if x == 0: return -1.0
        return -(1+x)*log1p(x)/x
    return r*log(r)/(1-r)


def _df(r: float) -> float:
    x = r-1.0
    if abs(x) < 1e-4:
        return -0.5+x/3-x*x/4+x**3/5-x**4/6+x**5/7
    return (log(r)+1-r)/(1-r)**2


def triangle_ratio(M2: float, a: float, b: float) -> float:
    """Full UV triangle / (C_k lambda /2LOOP).

    K=(F(a/M2)-F(b/M2))/(a/M2-b/M2), F(r)=r ln r/(1-r).
    This is the EXACT zero-momentum three-propagator kernel for one massive
    fermion and two real scalar propagators, once the LLss-tree convention is
    fixed. Analytically K -> -1/2 when all three masses coincide, and
    K -> 1+ln(mS^2/M2) when both scalar masses are small and equal.
    """
    M2,a,b=map(float,(M2,a,b))
    if min(M2,a,b)<=0 or not all(map(isfinite,(M2,a,b))):
        raise ValueError('Mass squares must be positive finite')
    x,y=a/M2,b/M2
    if abs(x-y) > 1e-5*max(x,y,1e-100):
        return (_f(x)-_f(y))/(x-y)
    # An integral of the derivative is the exact divided difference. Avoid
    # subtracting nearly equal logarithms for coincident scalar eigenvalues.
    nodes,weights=np.polynomial.legendre.leggauss(32)
    vals=[_df(0.5*((1-t)*x+(1+t)*y)) for t in nodes]
    return float(0.5*np.dot(weights,vals))


def weinberg_tensor() -> np.ndarray:
    eps=np.array([[0.,1.],[-1.,0.]],complex)
    h=np.array([[1.,1j,0.,0.],[0.,0.,1.,1j]],complex)/(2**0.5)
    v=eps@h
    t=np.einsum('ia,jb->ijab',v,v)
    return (t+t.swapaxes(2,3))/2


def project(amplitude: np.ndarray) -> tuple[np.ndarray,float]:
    a=np.asarray(amplitude,dtype=complex)
    if a.shape!=(3,3,2,2,4,4): raise ValueError('Unexpected LLHH shape')
    t=weinberg_tensor()
    c=np.einsum('fgijAB,ijAB->fg',a,t.conj(),optimize=True)/np.vdot(t,t)
    rest=a-np.einsum('fg,ijAB->fgijAB',c,t,optimize=True)
    rel=float(np.linalg.norm(rest)/max(np.linalg.norm(a),1e-300))
    sym=float(np.linalg.norm(c-c.T)/max(np.linalg.norm(c),1e-300))
    if sym>1e-9: raise ValueError('C5 projection is not flavour symmetric')
    return (c+c.T)/2,rel


def _finite_array(x: Any, name: str, *, complex_ok: bool = True) -> np.ndarray:
    a=np.asarray(x,dtype=complex if complex_ok else float)
    if not np.isfinite(a).all(): raise ValueError(name+' non-finite')
    return a


def evaluate_hard_scalar_channel(
    C_by_fermion: np.ndarray, fermion_masses_gev: np.ndarray,
    lambda_hhss: np.ndarray, scalar_mass_squared_gev2: np.ndarray,
    mu_f_gev: float, mu_s_gev: float, *, projection_tolerance: float=1e-8,
    C_total_muF: np.ndarray | None = None,
    lambda_channels: np.ndarray | None = None,
    channel_names: tuple[str, ...] | None = None,
) -> dict:
    """Compute separate one-loop hard, scalar RGE log, joint soft and full UV triangle.

    Inputs MUST be physical-mass-eigenstate contributions to the LLSS tree
    boundary, not copies of the total tensor. Both heavy scalar multiplets
    appear in one real mass matrix; no scalar or fermion is thresholded alone.
    Fermions are all matched together at muF; Mk only labels terms in a sum.

    No claimed physical input normalization is supplied by this function.
    The returned ``dim5_sum`` excludes any other six RGE terms or different
    hard topologies, and ``higher_dimension_remainder`` is not silently
    promoted to a SMEFT dimension-five coefficient.
    """
    c=_finite_array(C_by_fermion,'C_by_fermion')
    ms=_finite_array(scalar_mass_squared_gev2,'scalar mass squared',complex_ok=False)
    q=_finite_array(lambda_hhss,'lambda_HHSS',complex_ok=False)
    masses=_finite_array(fermion_masses_gev,'fermion masses',complex_ok=False)
    mu_f,mu_s=float(mu_f_gev),float(mu_s_gev)
    if not all(map(isfinite,(mu_f,mu_s))) or not mu_f>=mu_s>0:
        raise ValueError('Require one muF >= one muS > 0')
    if c.ndim!=7 or c.shape[1:5]!=(3,3,2,2):
        raise ValueError('C_by_fermion must be (k,3,3,2,2,n,n)')
    k,n=c.shape[0],c.shape[-1]
    if c.shape[-2]!=n or ms.shape!=(n,n) or q.shape!=(4,4,n,n) or masses.shape!=(k,):
        raise ValueError('Inconsistent complete scalar/fermion tensor dimensions')
    if k<1 or n<2 or np.any(masses<=0): raise ValueError('Empty/invalid masses')
    if not np.allclose(c,c.transpose(0,2,1,4,3,5,6),rtol=1e-9,atol=1e-15):
        raise ValueError('LLSS lepton exchange symmetry failed')
    if not np.allclose(c,c.swapaxes(-1,-2),rtol=1e-9,atol=1e-15):
        raise ValueError('LLSS scalar exchange symmetry failed')
    if not np.allclose(q,q.swapaxes(0,1),rtol=1e-9,atol=1e-15) or not np.allclose(q,q.swapaxes(2,3),rtol=1e-9,atol=1e-15):
        raise ValueError('HHSS tensor symmetry failed')
    if not np.allclose(ms,ms.T,rtol=1e-12,atol=1e-10):
        raise ValueError('Scalar mass2 must be real symmetric')
    if lambda_channels is not None:
        qc=_finite_array(lambda_channels,'lambda_HHSS channels',complex_ok=False)
        if qc.ndim!=5 or qc.shape[1:]!=q.shape or not np.allclose(qc.sum(axis=0),q,rtol=1e-11,atol=1e-14):
            raise ValueError('Independent quartic channels must exactly sum to lambda_HHSS')
        if channel_names is None or len(channel_names)!=len(qc) or len(set(channel_names))!=len(qc):
            raise ValueError('Unique named quartic channel required for every tensor')
    elif channel_names:
        raise ValueError('Quartic channel names require channel tensors')
    if C_total_muF is not None:
        total=_finite_array(C_total_muF,'C_total_muF')
        if total.shape!=c.shape[1:] or not np.allclose(total,c.sum(axis=0),rtol=1e-10,atol=1e-15):
            raise ValueError('Per-fermion LLSS sources do not sum to source total')
    ms_e,U=np.linalg.eigh(ms)
    if np.min(ms_e)<=0: raise ValueError('Tachyonic heavy scalar eigenstate')
    cr=np.einsum('kfgijab,ac,bd->kfgijcd',c,U,U,optimize=True)
    qr=np.einsum('ABab,ac,bd->ABcd',q,U,U,optimize=True)
    hard=np.zeros((3,3,2,2,4,4),complex)
    soft=np.zeros_like(hard)
    full=np.zeros_like(hard)
    rg=np.zeros_like(hard)
    terms=[]
    for kk,m in enumerate(masses):
        for a in range(n):
            for b in range(a,n):
                mult=1 if a==b else 2
                w=mult*np.einsum('fgij,AB->fgijAB',cr[kk,...,a,b],qr[...,a,b],optimize=True)/(2*LOOP)
                kh=1+2*log(mu_f/m)
                ks=-b0bar(ms_e[a],ms_e[b],mu_s)
                kf=triangle_ratio(m*m,ms_e[a],ms_e[b])
                kr=2*log(mu_s/mu_f)
                hard+=w*kh
                soft+=w*ks
                full+=w*kf
                rg+=w*kr
                projected_w,projected_residual=project(w)
                terms.append({'projected_C5_base_not_independently_gauge_invariant':
                              {'real':projected_w.real.tolist(),'imag':projected_w.imag.tolist()},
                              'individual_projection_residual':projected_residual,
                              'kernels_apply_to_full_complex_base_without_conjugation':True,
                              'fermion':kk,'scalar_pair':[a,b],'multiplicity':mult,
                              'mass_f_gev':float(m),'mass_scalar_pair_gev':
                              [float(ms_e[a]**.5),float(ms_e[b]**.5)],
                              'hard_kernel':float(kh),'soft_kernel':float(ks),
                              'full_triangle_kernel':float(kf),
                              'rg_scalar_log_kernel':float(kr),
                              'hard_amplitude_frobenius_gev_inv':float(np.linalg.norm(w*kh))})
    arrays={name:project(v) for name,v in {'hard_muF':hard,'eft1_scalar_log':rg,
             'scalar_threshold_muS':soft,'dim5_sum':hard+rg+soft,
             'full_triangle':full,'higher_dimension_remainder':full-hard-rg-soft}.items()}
    max_res=max(r for _,r in arrays.values())
    # The higher-dimensional remainder need not be exactly pure Weinberg when
    # EWSB-split masses were supplied. Such residuals cannot be forced to zero.
    out={key:value for key,(value,_) in arrays.items()}
    channels={}
    if lambda_channels is not None:
        for label,tensor in zip(channel_names,lambda_channels):
            part=evaluate_hard_scalar_channel(c,masses,tensor,ms,mu_f,mu_s,
                projection_tolerance=projection_tolerance,C_total_muF=c.sum(axis=0))
            channels[label]={'matrices':part['matrices'],
                             'projection_residuals':part['residuals'],
                             'physical_C5_certified':False}
        for key in out:
            if not np.allclose(sum((v['matrices'][key] for v in channels.values()),
                                   np.zeros((3,3),complex)),out[key],rtol=1e-9,atol=1e-18):
                raise AssertionError('Quartic-source sum failed for '+key)
    return {'matrices':out,'quartic_channels':channels,'residuals':{key:r for key,(_,r) in arrays.items()},
            'projection_passed':bool(max_res<=projection_tolerance),
            'mass_hierarchy_max_scalar_over_min_fermion':float(np.sqrt(max(ms_e))/min(masses)),
            'dimension_five_truncation':True,'all_fermions_one_muF':True,
            'both_scalars_one_muS':True, 'per_pair_terms':terms,
            'physical_C5_certified':False,
            'complete_hard_one_loop_matching_certified':False,
            'other_six_EFT1_RGE_sources_included':False,
            'higher_dimension_remainder_loop_order':1,
            'higher_dimension_remainder_promoted_to_dim5':False,
            'normalization_requirement': 'independent signed GammaCC-to-Weyl factor and per-fermion tree LLSS tensors'}
