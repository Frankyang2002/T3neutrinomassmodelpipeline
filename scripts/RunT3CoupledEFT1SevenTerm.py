"""Exploratory A(-4) seven-term coupled dimension-5 EFT1 RGE.

Intended path: scripts/RunT3CoupledEFT1SevenTerm.py

No existing source files, authoritative datasets, or GitHub files are changed.
Reuses RunT3GeneralSequentialFitAblation's fitted UV/EFT1 trajectories,
fermion-tree LLSS tensor, HHSS quartic tensor and final SMEFT integration.

This runner genuinely integrates all seven terms of the repository's
RGE.running.intermediate.GeneralEFT1C5Mixing.calculate_master_beta
for a FULL complex symmetric master C[fermion,fermion,scalar,scalar] tensor.
The dynamical master coefficient includes induced components initially zero.
The evolving Wilson tensor is restricted to 9 lepton Weyl fields (L and
e^c); quark Yukawa contributions to the Higgs scalar anomalous dimension
are added explicitly using all three families and colours. This is a closed
subsystem when no lepton-quark Yukawas or BSM-scalar Yukawas remain in EFT1.

A(-4) has no extra candidate quartics in the general registry. Its HHHH,
SSSS and mixed S1-S2 quartic blocks are constructed with the explicitly
DOCUMENTED standard-potential ansatz below. HHSS is taken without alteration
from the exact same numerical conditional tensor as the corrected LLSS run.

IMPORTANT LIMITS (do not call this complete physical matching):
 * Renormalisable couplings/Yukawas/quartics are frozen to EFT1 muS endpoint
   while Wilson C is integrated. This isolates Wilson self-running and mixes
   all seven beta terms; no continuously changing coupling background.
 * The conversion of SM Yukawa matrices to the left-Weyl master y_ija is a
   stated candidate convention; heavy-source GammaCC sign/factor remains
   the inherited explicit spinor_factor.
 * The S2 triplet adjoint self quartic and S1/S2 scalar potential use the
   stated baseline potential ansatz, not a newly source-proven dictionary.
 * No hard one-loop matching or operators outside psi^2 phi^2 are supplied.
 * The existing finite scalar bubble is reused, NOT recomputed with the
   full seven-term Wilson tensor. A larger solution is not automatically
   a physically complete one-loop calculation.

PowerShell from repository root (prior quartic experiment installed):
 python -X utf8 -m scripts.RunT3CoupledEFT1SevenTerm --self-test
 python -X utf8 -m scripts.RunT3CoupledEFT1SevenTerm --source-run thesis16_merged_v1 --run-id seven_A_m4_v1 --models A_-4 --check-only
 python -X utf8 -m scripts.RunT3CoupledEFT1SevenTerm --source-run thesis16_merged_v1 --run-id seven_A_m4_v1 --models A_-4
"""
from __future__ import annotations

import math
import sys
from itertools import permutations
from pathlib import Path
from typing import Mapping

import numpy as np
from scipy.integrate import solve_ivp

LOOP = 1.0 / (16.0 * math.pi**2)
OUT_SCHEMA = 't3-a-m4-seven-term-full-master-dim5-ode-v2'
N_F = 9  # only L:0..5 and e^c:6..8 participate in evolving C
N_F_FULL = 45 # full SM Yukawas used to compute quark Higgs anomalous dimension
N_S = 12 # H:0..3; S1 singlet:4..5; S2 triplet:6..11
RTOL, ATOL = 1e-8, 1e-11 # scale state by norm(C_LLSS) before integration


def _finite_complex_matrix(value, name):
    """Decode the SAME serialized SM Yukawa format as the running pipeline.

    The canonical `_state_values` encodes each matrix entry separately as
    `{"real": ..., "imag": ...}`. The general LLSS `decode` helper
    alone accepts a whole-array real/imag dictionary but NOT the nested-list
    layout. Reuse the established SMEFT boundary decoder so both formats work.
    """
    from scripts.RunT3GeneralSequentialLeadingLog import _complex_matrix
    try:
        arr = _complex_matrix(value)
    except (TypeError, ValueError, KeyError) as exc:
        raise ValueError(f'{name}: cannot decode three-family complex Yukawa') from exc
    if arr.shape != (3,3) or not np.isfinite(arr).all():
        raise ValueError(f'{name}: expected complex finite 3x3')
    return arr


def _get_real(state: Mapping, name: str) -> float:
    if name not in state or state[name] is None:
        raise ValueError(f'Missing A(-4) EFT1 endpoint coupling: {name}')
    v=float(state[name])
    if not math.isfinite(v):
        raise ValueError(f'Nonfinite {name}')
    return v


def su2_complex(d):
    j=(d-1)*0.5
    levels=np.array([j-k for k in range(d)])
    plus=np.zeros((d,d),complex)
    for col,m in enumerate(levels):
        indices=np.where(np.isclose(levels,m+1))[0]
        if indices.size:
            plus[int(indices[0]),col]=math.sqrt(max(0.0,(j-m)*(j+m+1)))
    return ((plus+plus.conj().T)/2,
            (plus-plus.conj().T)/(2j),
            np.diag(levels).astype(complex))


def scalar_real_generator(t):
    """Exactly the interleaved (Re,Im) convention in GaugeGenerators.py."""
    dim=t.shape[0]
    theta=np.zeros((2*dim,2*dim),complex)
    for p in range(dim):
        for q in range(dim):
            z=t[p,q]
            theta[2*p,2*q]+=1j*z.imag
            theta[2*p,2*q+1]+=1j*z.real
            theta[2*p+1,2*q]+=-1j*z.real
            theta[2*p+1,2*q+1]+=1j*z.imag
    return theta


def scalar_bilinear_real(t):
    """Symmetric real M for (S^dag t S)=1/2 r^T M r."""
    dim=t.shape[0]
    m=np.zeros((2*dim,2*dim),float)
    for p in range(dim):
        for q in range(dim):
            z=t[p,q]
            m[2*p,2*q]+=z.real
            m[2*p,2*q+1]+=-z.imag
            m[2*p+1,2*q]+=z.imag
            m[2*p+1,2*q+1]+=z.real
    if not np.allclose(m,m.T,atol=1e-12):
        raise ValueError('Adjoint scalar bilinear is not symmetric')
    return m


def gauge_context(state: Mapping):
    """SU2 x U1 generators on nine lepton Weyl fields and all 12 scalars.

    Quarks cannot become external Wilson legs from the retained LLSS source
    with only SM Yukawas. They enter gamma_scalar via their Higgs Yukawas.
    SU3 hence has no direct action on the retained Wilson-tensor subspace.
    """
    g2=_get_real(state,'g2');gY=_get_real(state,'gY')
    scalar_blocks=((0,2,0.5),(4,1,-2.0),(6,3,-1.0))
    scalar_su2=[]
    for axis in range(3):
        out=np.zeros((N_S,N_S),complex)
        for start,d,_ in scalar_blocks:
            out[start:start+2*d,start:start+2*d]=scalar_real_generator(su2_complex(d)[axis])
        scalar_su2.append(out)
    scalar_u1=np.zeros((N_S,N_S),complex)
    for start,d,y in scalar_blocks:
        scalar_u1[start:start+2*d,start:start+2*d]=scalar_real_generator(y*np.eye(d))
    fermion_su2=[]
    for axis in range(3):
        f=np.zeros((N_F,N_F),complex)
        t=su2_complex(2)[axis]
        for gen in range(3):
            f[2*gen:2*gen+2,2*gen:2*gen+2]=t
        fermion_su2.append(f)
    fermion_u1=np.diag([-0.5]*6+[1.0]*3).astype(complex)
    return ((g2,tuple(zip(fermion_su2,scalar_su2))),
            (gY,((fermion_u1,scalar_u1),)))


def yukawa_context(state: Mapping):
    """Candidate left-Weyl SM y_ija; no heavy F Yukawa remains in EFT1.

    L H^dag e^c and Q H^dag d^c use H*; Q.epsilon.H u^c uses H.
    Complex H=(Re+i Im)/sqrt(2), i,j symmetric as y_ija psi_i psi_j phi_a/2.
    Existing fitted SM ye,yu,yd flavour matrices are used directly (no refit).
    """
    ye=_finite_complex_matrix(state['ye'],'ye')
    yu=_finite_complex_matrix(state['yu'],'yu')
    yd=_finite_complex_matrix(state['yd'],'yd')
    y=np.zeros((N_F_FULL,N_F_FULL,N_S),complex)
    def add(i,j,comp,z):
        y[i,j,comp]+=z; y[j,i,comp]+=z
    inv=1/math.sqrt(2)
    for g in range(3):
        for h in range(3):
            e=6+h
            add(2*g,e,0,ye[g,h]*inv)
            add(2*g,e,1,-1j*ye[g,h]*inv)
            add(2*g+1,e,2,ye[g,h]*inv)
            add(2*g+1,e,3,-1j*ye[g,h]*inv)
            for color in range(3):
                q=9+2*(3*g+color)
                u=27+3*h+color
                d=36+3*h+color
                # epsilon contraction Q_up H0 - Q_down H+.
                add(q,u,2,yu[g,h]*inv)
                add(q,u,3,1j*yu[g,h]*inv)
                add(q+1,u,0,-yu[g,h]*inv)
                add(q+1,u,1,-1j*yu[g,h]*inv)
                # Q_up H+* + Q_down H0*.
                add(q,d,0,yd[g,h]*inv)
                add(q,d,1,-1j*yd[g,h]*inv)
                add(q+1,d,2,yd[g,h]*inv)
                add(q+1,d,3,-1j*yd[g,h]*inv)
    return y


def quark_higgs_anomalous_dimension(full_y):
    """Exactly the quark Yukawa piece in master _collinear scalar gamma."""
    qy=full_y[9:,9:,:]
    gram=np.einsum("ija,ijb->ab",qy,qy.conj(),optimize=True)
    return 0.5*(gram+gram.T)


def add_pair_contraction(lam,p,q,k):
    """For V=(k/8)*(r.T@p@r)*(r.T@q@r) if p=q;
    for distinct p,q use V=(k/4)*(r.T@p@r)*(r.T@q@r).
    This routine receives k so that p=q contributes k*(p_ab p_cd + ...),
    and p!=q contributes k*(p_ab q_cd + q_ab p_cd + ...).
    """
    lam += k*(np.einsum('ab,cd->abcd',p,q)+
              np.einsum('ac,bd->abcd',p,q)+
              np.einsum('ad,bc->abcd',p,q))
    if p is not q:
        lam += k*(np.einsum('ab,cd->abcd',q,p)+
                  np.einsum('ac,bd->abcd',q,p)+
                  np.einsum('ad,bc->abcd',q,p))


def full_scalar_quartic_a4(state: Mapping, hhss: np.ndarray):
    """HHSS from existing source; baseline HHHH / SSSS from explicit ansatz.

    V4 (non-HHSS, candidate):
      + (lambdaH/2)(H†H)^2 + (lambdaS1/2)(S1†S1)^2
      + (lambdaS2/2)(S2†S2)^2
      + lambda12 (S1†S1)(S2†S2)
      + (lambdaS2Adj/2) sum_A(S2†T^A S2)^2.
    For phi=(r+i i)/sqrt(2), lambda_abcd=∂^4 V4/∂r_a... .
    HHSS coefficients use exactly scalar.quartic, avoiding a second map.
    A(-4) extras are empty in general_rgbeta_quartic_registry.json.
    """
    q=np.asarray(hhss,dtype=float)
    if q.shape!=(4,4,8,8):
        raise ValueError(f'A(-4) expected HHSS shape (4,4,8,8); got {q.shape}')
    lam=np.zeros((N_S,)*4,float)
    for A in range(4):
        for B in range(4):
            for a in range(8):
                for b in range(8):
                    value=q[A,B,a,b]
                    if value==0: continue
                    for tup in set(permutations((A,B,4+a,4+b))):
                        lam[tup]=value
    h=np.diag([1.]*4+[0.]*8)
    s1=np.diag([0.]*4+[1.]*2+[0.]*6)
    s2=np.diag([0.]*6+[1.]*6)
    add_pair_contraction(lam,h,h,_get_real(state,'lambdaH'))
    add_pair_contraction(lam,s1,s1,_get_real(state,'lambdaS1'))
    add_pair_contraction(lam,s2,s2,_get_real(state,'lambdaS2'))
    # Distinct-block conventional cross potential lambda12 |S1|² |S2|².
    # add_pair_contraction has 6 terms; lambda12 gives 1 for S1S1S2S2.
    add_pair_contraction(lam,s1,s2,_get_real(state,'lambda12'))
    c_adj=_get_real(state,'lambdaS2Adj')
    for t in su2_complex(3):
        adj=np.zeros((N_S,N_S))
        adj[6:12,6:12]=scalar_bilinear_real(t)
        add_pair_contraction(lam,adj,adj,c_adj)
    if not np.all(np.isfinite(lam)):
        raise ValueError('Nonfinite master quartic')
    for ax in range(1,4):
        if not np.allclose(lam,np.swapaxes(lam,0,ax),rtol=1e-10,atol=1e-12):
            raise ValueError('Master quartic tensor violates Bose symmetry')
    if not np.allclose(lam[:4,:4,4:,4:],q,rtol=1e-10,atol=1e-12):
        raise ValueError('HHSS baseline overwritten by other quartics')
    return lam


def embed_llss(source):
    a=np.asarray(source,dtype=complex)
    if a.shape!=(3,3,2,2,8,8):
        raise ValueError('Expected complete A(-4) LLSS (3,3,2,2,8,8)')
    c=np.zeros((N_F,N_F,N_S,N_S),complex)
    for p in range(3):
        for q in range(3):
            for i in range(2):
                for j in range(2):
                    c[2*p+i,2*q+j,4:,4:]=a[p,q,i,j]
    if not np.allclose(c,c.swapaxes(0,1),rtol=1e-11,atol=1e-20):
        raise ValueError('Flattened LLSS not symmetric in Weyl indices')
    if not np.allclose(c,c.swapaxes(2,3),rtol=1e-11,atol=1e-20):
        raise ValueError('Flattened LLSS not symmetric in scalar indices')
    return c


def project_leptonic_weinberg(master_tensor):
    from Numerical.joint_scalar_matching.core import project_weinberg
    z=np.asarray(master_tensor,complex)
    ll=z[np.ix_(np.arange(6),np.arange(6),np.arange(4),np.arange(4))]
    six=ll.reshape(3,2,3,2,4,4).transpose(0,2,1,3,4,5)
    return project_weinberg(six,strict=False)


def initial_context(state: Mapping, hhss: np.ndarray, source: np.ndarray,
                    mu_f: float,mu_s: float):
    from RGE.running.intermediate.GeneralEFT1C5Mixing import Context
    quartic=full_scalar_quartic_a4(state,hhss)
    full_y=yukawa_context(state)
    quark_gamma=quark_higgs_anomalous_dimension(full_y)
    context=Context(
        model='A_-4', nf=N_F, ns=N_S,
        leptons=((0,1),(2,3),(4,5)), higgs=(0,1,2,3),
        mu_f=mu_f, mu_s=mu_s,
        yukawa=full_y[:N_F,:N_F,:],quartic=quartic,
        gauge=gauge_context(state),
        operators=(('source_resolved_fermion_tree_LLSS',embed_llss(source)),),
        authority='synthetic',certificates={},
    )
    return context,quark_gamma


def solve_master(context,*,quark_gamma=None,rtol=RTOL,atol=ATOL):
    from RGE.running.intermediate.GeneralEFT1C5Mixing import (
        TERMS, calculate_master_beta,
    )
    c0=context.operators[0][1]
    norm=float(np.linalg.norm(c0))
    if not math.isfinite(norm) or norm<=0: raise ValueError('Zero master Wilson boundary')
    init=np.asarray(c0/norm,complex).ravel()
    shape=c0.shape
    def beta_terms(c):
        terms=calculate_master_beta(context,c)
        if quark_gamma is not None:
            gam=np.einsum("bd,ijad->ijab",quark_gamma,c,optimize=True)
            terms["scalar_anomalous_dimension"] += gam+gam.swapaxes(2,3)
        return terms
    initial_terms=beta_terms(c0)
    first_order=sum(initial_terms.values(),np.zeros_like(c0))
    defect=float(np.linalg.norm(first_order-first_order.swapaxes(0,1)) /
                 max(np.linalg.norm(first_order),1e-300))
    defect_scalars=float(np.linalg.norm(first_order-first_order.swapaxes(2,3)) /
                         max(np.linalg.norm(first_order),1e-300))
    if max(defect,defect_scalars)>1e-7:
        raise ValueError('Master beta violates coefficient symmetries: '
                         f'fermions={defect}, scalars={defect_scalars}')

    def rhs(t,v):
        c=v.reshape(shape)
        # beta of a complex tensor is real-linear: it can involve c.conj().
        terms=beta_terms(c)
        derivative=sum(terms.values(),np.zeros_like(c))
        return LOOP*derivative.ravel()
    logf,logs=math.log(context.mu_f),math.log(context.mu_s)
    sol=solve_ivp(rhs,(logf,logs),init,rtol=rtol,atol=atol,method='DOP853')
    if not sol.success: raise RuntimeError(f'Full master DOP853 failed: {sol.message}')
    end=norm*sol.y[:,-1].reshape(shape)
    if not np.isfinite(end).all(): raise ValueError('Master tensor integrated to nonfinite')
    logsum=LOOP*(logs-logf)*first_order
    return {
        'end':end,'first_order_tensor':logsum,'first_terms':initial_terms,
        'nfev':int(sol.nfev),'nsteps':len(sol.t),
        'tensor_relative_drift':float(np.linalg.norm(end-c0)/norm),
        'fermion_symmetry_residual':float(np.linalg.norm(end-end.swapaxes(0,1))/norm),
        'scalar_symmetry_residual':float(np.linalg.norm(end-end.swapaxes(2,3))/norm),
    }


def self_test():
    """Check Yukawa decoding, seven-term evaluator and coupled DOP853."""
    # Regression: `_state_values` emits a list of per-entry complex dicts.
    sample = np.array([[complex(0.03*i, -0.02*j) for j in range(3)]
                       for i in range(3)], dtype=complex)
    nested = [[{'real': float(z.real), 'imag': float(z.imag)}
               for z in row] for row in sample]
    dense = {'real': sample.real.tolist(), 'imag': sample.imag.tolist()}
    for raw in (nested, dense):
        decoded = _finite_complex_matrix(raw, 'test_yukawa')
        if not np.array_equal(decoded, sample):
            raise AssertionError('Complex Yukawa deserialization disagrees with input')
    print('PASS: fitted Yukawa nested-complex and dense-complex decoding')
    from RGE.running.intermediate.GeneralEFT1C5Mixing import Context
    q=np.zeros((4,4,2,2));q[2,2,0,0]=0.5
    scalar=np.zeros((6,6,6,6))
    for tup in set(permutations((2,2,4,4))):scalar[tup]=0.5
    C=np.zeros((6,6,6,6),complex);C[0,0,4,4]=1.e-13
    ctx=Context('A_-4',6,6,((0,1),(2,3),(4,5)),(0,1,2,3),
                1.e5,1.e3,np.zeros((6,6,6),complex),scalar,(),
                (('synthetic',C),),'synthetic',{})
    out=solve_master(ctx)
    k=LOOP*0.5*math.log(1.e3/1.e5)
    got=out['end'][0,0,2,2]/1.e-13
    remaining=out['end'][0,0,4,4]/1.e-13
    if not (np.isclose(got,math.sinh(k),rtol=1e-7,atol=1e-12) and
            np.isclose(remaining,math.cosh(k),rtol=1e-7,atol=1e-12)):
        raise AssertionError('Seven-term master ODE fails analytic quartic-only limit')
    if any(np.linalg.norm(out['first_terms'][term])>1.e-25
           for term in out['first_terms'] if term!='scalar_pair'):
        raise AssertionError('Yukawa/gauge terms nonzero with zero couplings')
    print('PASS: actual seven-term master evaluator + DOP853 matches sinh/cosh toy limit')
    print(f'  first-order expectation={k:.12e}; full result={got.real:.12e}')
    print(f'  nfev={out["nfev"]}')


def final_cases(model,c5_one,source,scalar,target,*,
                uv_on,eft_on,uv_nfev,eft_nfev,spinor_factor,vev):
    if model!='A_-4': raise ValueError('Seven-term experiment currently only A_-4')
    from Numerical.running.GeneralT3SequentialLeadingLog import encode,decode,response
    from scripts import RunT3GeneralSequentialFitAblation as original
    from scripts.RunT3CoupledEFT1Quartic import integrate_coupled_quartic
    from RGE.running.intermediate.GeneralEFT1C5Mixing import TERMS
    mu_f,mu_s=original.SCALES[1:3]
    c0=np.asarray(source.raw,complex)*complex(spinor_factor)
    quartic=np.asarray(scalar.quartic,float)
    state=original.scalar_snapshot_state(scalar)
    fixed=response(source,scalar,mu_f=mu_f,mu_s=mu_s,
                   spinor_factor=complex(spinor_factor),
                   check_scale_cancellation=True)
    base_log=decode(fixed['c5_eft1_leading_log'])
    finite=decode(fixed['c5_scalar_finite_bubble'])
    quartic_only=integrate_coupled_quartic(c0,quartic,mu_f,mu_s)
    # Reuse EXACT original project_weinberg convention with six-index tensor.
    from Numerical.joint_scalar_matching.core import project_weinberg
    quartic_c5,quartic_ward=project_weinberg(quartic_only['c_llhh_end'],strict=False)
    ctx,quark_gamma=initial_context(state,quartic,c0,mu_f,mu_s)
    sol=solve_master(ctx,quark_gamma=quark_gamma)
    c5_at_f,ward=project_leptonic_weinberg(sol['end'])
    c5_first,first_ward=project_leptonic_weinberg(sol['first_order_tensor'])
    c5_first_quartic,first_q_ward=project_leptonic_weinberg(
        LOOP*math.log(mu_s/mu_f)*sol['first_terms']['scalar_pair'])
    baseline_err=float(np.linalg.norm(c5_first_quartic-base_log)/max(np.linalg.norm(base_log),1e-300))
    if baseline_err>1e-8:
        raise ValueError('Master initial quartic beta does not reproduce previous LLSS leading log: '
                         f'relative difference={baseline_err:.5e}')
    term_projected={}
    for key in TERMS:
        projected,_=project_leptonic_weinberg(sol['first_terms'][key])
        term_projected[key]=encode(LOOP*math.log(mu_s/mu_f)*projected)
    relative=lambda a,b:float(np.linalg.norm(a-b)/max(np.linalg.norm(b),1e-300))
    stats={
        'full_7term_vs_fixed_log_rel':relative(c5_at_f,base_log),
        'full_7term_vs_quartic_ode_rel':relative(c5_at_f,quartic_c5),
        'full_7term_vs_full_initial_slope_rel':relative(c5_at_f,c5_first),
        'full_initial_slope_vs_quartic_fixed_rel':relative(c5_first,base_log),
        'quartic_only_ode_vs_fixed_log_rel':relative(quartic_c5,base_log),
        'master_quartic_initial_slope_parity_rel':baseline_err,
        'master_weinberg_residual':float(ward),
        'master_first_order_weinberg_residual':float(first_ward),
        'master_scalar_symmetry_residual':sol['scalar_symmetry_residual'],
        'master_fermion_symmetry_residual':sol['fermion_symmetry_residual'],
        'master_full_tensor_drift_rel':sol['tensor_relative_drift'],
        'master_ode_nfev':sol['nfev'],'master_ode_steps':sol['nsteps'],
    }
    # The original four-bit nomenclature is retained. W=1 is now the full
    # seven-term Wilson contribution, W=0 the same finite scalar bubble only.
    cases={}
    for w_on, c5 in ((False,finite),(True,finite+c5_at_f)):
        c5=original.c5_checked(c5)
        for sm_on in (False,True):
            key=f'{int(uv_on)}{int(eft_on)}{int(w_on)}{int(sm_on)}'
            out=original.low_energy(state,c5,on=sm_on,mu_s=mu_s,
                                    mu_low=original.SCALES[3],vev=vev,target=target)
            out.update({
                'case':key,'uv_on':uv_on,'eft1_on':eft_on,
                'wilson_log_on':w_on,'smeft_on':sm_on,
                'uv_nfev':uv_nfev,'eft1_nfev':eft_nfev,
                'c5_at_scalar_threshold':encode(c5),
                'c5_induced_wilson_log':encode(c5_at_f),
                'c5_fixed_order_llss_log':encode(base_log),
                'c5_coupled_quartic_only':encode(quartic_c5),
                'c5_initial_seven_term_slope':encode(c5_first),
                'c5_joint_scalar_bubble':encode(finite),
                'c5_matching_anchored_reference':encode(original.c5_checked(c5_one)),
                'c5_bookkeeping_residual_matching':encode(np.zeros((3,3),complex)),
                'c5_log_frobenius':float(np.linalg.norm(c5_at_f)),
                'c5_one_scale_frobenius':float(np.linalg.norm(c5_one)),
                'c5_log_over_one_scale_norm':float(np.linalg.norm(c5_at_f)/
                                               max(np.linalg.norm(finite+c5_at_f),1e-300)),
                'per_rge_term_initial_slope_log_contribution':term_projected,
                'per_general_interaction':fixed['per_channel'],
                'quartic_provenance':fixed['quartic_provenance'],
                'matching_scale_cancellation':fixed['matching_scale_cancellation_residual'],
                'full_master_rge_provisional':True,
                'renormalisable_couplings_frozen_to_muS':True,
                'wilson_fermion_subspace':'9 lepton Weyl fields; quark Higgs gamma included',
                'full_physical_one_loop_matching':False,
                'scalar_potential_convention':
                    'lambdaS1,S2/2 |S|^4 + lambda12 |S1|^2|S2|^2 + lambdaS2Adj/2 |S2^dag T S2|^2',
                'weyl_yukawa_embedding':
                    'candidate SM ye,yu,yd full 3-gen 3-colour left-chiral convention',
                'matching_approximation':
                    'exploratory_full_master_7terms_fixed_muS_background_plus_original_finite_bubble',
                **stats,
            })
            cases[key]=out
    print(f'{model} U={int(uv_on)} E={int(eft_on)} : '
          f'7term/fixed={stats["full_7term_vs_fixed_log_rel"]:.5e} '
          f'7term/quarticODE={stats["full_7term_vs_quartic_ode_rel"]:.5e} '
          f'quarticODE/fixed={stats["quartic_only_ode_vs_fixed_log_rel"]:.5e} '
          f'7term_initial_parity={baseline_err:.3e} '
          f'DOP853 nfev={sol["nfev"]}', flush=True)
    return cases


def main(argv=None):
    args=list(sys.argv[1:] if argv is None else argv)
    if args==['--self-test']:
        self_test();return 0
    if '--run-id' not in args:raise ValueError('Use a fresh --run-id')
    if '--resume' in args:raise ValueError('Resume forbidden: use fresh --run-id')
    if '--models' not in args or args[args.index('--models')+1]!='A_-4':
        raise ValueError('Seven-term exploratory run currently requires --models A_-4')
    from scripts import RunT3GeneralSequentialFitAblation as original
    original._final_cases=final_cases
    original.OUT=original.ROOT/'output/studies/general_scalar_extended/experimental_coupled_eft1_seven_term'
    old_save=original.save
    def save_with_metadata(path,data):
        if isinstance(data,dict):
            data=dict(data)
            if data.get('format')=='t3-full-general-16_U_E_W_S_ablation-v1':
                data['format']=OUT_SCHEMA
                data['seven_master_terms']=True
                data['running_background']='frozen_muS_EFT1_endpoint'
                data['complete_physical_C5']=False
            elif data.get('status')=='Computed' and 'cases' in data:
                data['full_sequential_general_C5_certified']=False
                data['matching_approximation']='experimental_seven_term_master_ode_A_m4'
        old_save(path,data)
    original.save=save_with_metadata
    return original.main(args)


if __name__=='__main__':
    raise SystemExit(main())
