"""Conditional but fully numerical, one-scale real-scalar bubble threshold.

Canonical master Wilson polynomial has no 1/2 or 1/4 factorial:
  L5 = C_ijab psi_i psi_j phi_a phi_b + h.c.
  L5^W = (1/2) C5_pq (L_p.eps.H)(L_q.eps.H) + h.c.
  C[nu_p,nu_q,H0R,H0R] = C5_pq/4.
V4 = lambda_abcd phi_a phi_b phi_c phi_d / 4!.

For the scalar bubble use the project's convention-conditional formula
  Delta C_ijAB = -(1/(2*16*pi^2)) sum_ab C_ijab lambda_ABab B0bar(0,ma2,mb2;MS)
with ALL heavy scalar eigenstate pairs included at the SINGLE MS threshold.
Other one-loop diagrams, field-normalisation/EOM effects and native-to-canonical
matching remain uncertified. No old direct LLSS->C5 leading log is added.
"""
from __future__ import annotations
import math
import numpy as np

LOOP=16.*math.pi**2


def finite(a: object, label: str, *, real: bool=False) -> np.ndarray:
    v=np.asarray(a)
    if not np.isfinite(v.real).all() or not np.isfinite(v.imag).all():
        raise ValueError(f'Nonfinite {label}')
    if real and np.max(abs(v.imag)) > 1e-12*max(1.,float(np.max(abs(v.real)))):
        raise ValueError(f'{label} must be real in real-field coordinates')
    return np.asarray(v.real if real else v,dtype=float if real else complex)


def symmetric3(x: np.ndarray, label: str) -> np.ndarray:
    a=finite(x,label)
    if a.shape != (3,3): raise ValueError(f'{label} must be complex symmetric 3x3')
    if np.linalg.norm(a-a.T)>1e-9*max(np.linalg.norm(a),1e-40):
        raise ValueError(f'{label} violates complex flavour symmetry (no Hermitian conjugation)')
    return (a+a.T)/2


def master_weinberg_template() -> np.ndarray:
    """Coefficient of C5 in full-sum master polynomial, with HH real legs.

    H=(h0+ih1,h2+ih3)/sqrt2, epsilon12=+1.
    This is HALF of joint_scalar_matching/core.py:weinberg_tensor():
    the latter uses a different external amplitude normalisation.
    """
    z=np.zeros((2,4),complex)
    z[0,2:]=[1/np.sqrt(2),1j/np.sqrt(2)]
    z[1,:2]=[-1/np.sqrt(2),-1j/np.sqrt(2)]
    return (np.einsum('iA,jB->ijAB',z,z)+
            np.einsum('iB,jA->ijAB',z,z))/4


def project_master(a: np.ndarray, *, tolerance: float=1e-8,
                   enforce_pure: bool=True) -> tuple[np.ndarray,float]:
    a=finite(a,'LLHH amplitude')
    if a.shape!=(3,3,2,2,4,4):
        raise ValueError('Expected [flavour,flavour,weak,weak,Hreal,Hreal]')
    T=master_weinberg_template()
    c5=np.einsum('fgijAB,ijAB->fg',a,T.conjugate(),optimize=True)/np.vdot(T,T)
    c5=symmetric3(c5,'projected Weinberg C5')
    err=a-np.einsum('fg,ijAB->fgijAB',c5,T,optimize=True)
    norm=float(np.linalg.norm(a))
    residual=float(np.linalg.norm(err)/max(norm,1e-100))
    if enforce_pure and residual>tolerance:
        raise ValueError(f'LLHH violates pure SU(2)-invariant Weinberg template: residual={residual:.4g}')
    return c5,residual


def b0bar(ma2: float, mb2: float, mu: float) -> float:
    """Finite MSbar zero-momentum two-scalar integral; squared masses GeV^2."""
    a,b,mu=float(ma2),float(mb2),float(mu)
    if not np.isfinite([a,b,mu]).all() or min(a,b,mu)<=0:
        raise ValueError('Scalar eigenmass^2 and matching scale must be finite positive')
    c=(a+b)/2; d=(a-b)/(a+b)
    if abs(d)<1e-3:
        t=d*d
        return -math.log(c/mu**2)+t/6+t*t/20+t**3/42+t**4/72
    return 1.-(a*math.log(a/mu**2)-b*math.log(b/mu**2))/(a-b)


def mass_eigenbasis(mass2: np.ndarray) -> tuple[np.ndarray,np.ndarray]:
    M=finite(mass2,'mass_squared_gev2',real=True)
    if M.ndim!=2 or M.shape[0]!=M.shape[1] or M.shape[0]<2:
        raise ValueError('Real scalar mass-squared matrix must be square, dimension >=2')
    if np.linalg.norm(M-M.T)>1e-10*max(np.linalg.norm(M),1.):
        raise ValueError('Scalar mass-squared Hessian must be symmetric')
    m,U=np.linalg.eigh((M+M.T)/2)
    if m[0]<=0:raise ValueError(f'Tachyonic/massless integrated scalar: m2={m[0]:g} GeV^2')
    if not np.allclose(U.T@M@U,np.diag(m),rtol=1e-10,atol=1e-9*max(1,float(np.max(m)))):
        raise ValueError('Mass eigenbasis diagonalisation failed')
    return m,U


def check_rotation(R: np.ndarray, n: int) -> np.ndarray:
    R=finite(R,'R_wilson_from_background',real=True)
    if R.shape!=(n,n) or not np.allclose(R.T@R,np.eye(n),atol=1e-10,rtol=1e-10):
        raise ValueError('A/B real-scalar basis transition must be explicitly supplied orthogonal nxn')
    return R


def check_quartic(q_full: np.ndarray, n: int) -> np.ndarray:
    full=finite(q_full,'full lambda_abcd',real=True)
    if full.shape!=(4+n,)*4:
        raise ValueError(f'Full scalar quartic must be real tensor {(4+n,)*4}')
    scale=max(float(np.max(abs(full))),1e-30)
    for permutation in ((1,0,2,3),(0,2,1,3),(0,1,3,2)):
        if np.max(abs(full-full.transpose(permutation)))>1e-9*scale:
            raise ValueError('V4=lam phi^4/4! needs a fully symmetric real quartic')
    return full



def scalar_gauge_ward_residuals(full_quartic: np.ndarray, mass2: np.ndarray,
                                generators_full: np.ndarray) -> dict:
    """Optional source-supplied *real-component* electroweak Ward diagnostic.

    G is the Hermitian purely imaginary generator of unitary gauge action on
    real fields; T=iG is a real antisymmetric infinitesimal rotation.
    The first four entries are the SM real Higgs, the remaining entries
    are the heavy background coordinates. C cannot generate missing G.
    """
    q=np.asarray(full_quartic,float)
    n=q.shape[0]-4
    m=np.asarray(mass2,float)
    gs=finite(generators_full,'full scalar gauge generators')
    if gs.ndim!=3 or gs.shape[1:]!=(n+4,n+4) or not len(gs):
        raise ValueError('Gauge generators must be [ng,4+n,4+n]')
    r_mass=0.;r_quartic=0.
    for g in gs:
        T=1j*g
        if np.max(abs(T.imag))>1e-10 or np.linalg.norm(T.real+T.real.T)>1e-9*max(1.,np.linalg.norm(T)):
            raise ValueError('Real-field gauge generator must be purely imaginary antisymmetric')
        T=T.real
        if np.max(abs(T[:4,4:]))>1e-10 or np.max(abs(T[4:,:4]))>1e-10:
            raise ValueError('Z2-even Higgs and Z2-odd heavy scalar generator blocks must not mix')
        t=T[4:,4:]
        r_mass=max(r_mass,float(np.linalg.norm(m@t-t@m)/max(np.linalg.norm(m),1.)))
        dQ=(np.einsum('ea,ebcd->abcd',T,q,optimize=True)+
            np.einsum('eb,aecd->abcd',T,q,optimize=True)+
            np.einsum('ec,abed->abcd',T,q,optimize=True)+
            np.einsum('ed,abce->abcd',T,q,optimize=True))
        r_quartic=max(r_quartic,float(np.linalg.norm(dQ)/max(np.linalg.norm(q),1.)))
    return {'mass_commutator_relative':r_mass,'quartic_invariance_relative':r_quartic,
            'generator_count':len(gs),'passed':max(r_mass,r_quartic)<1e-8}

def extract_lepton_sectors(full_C: np.ndarray, lepton_indices: np.ndarray,
                           higgs_indices: np.ndarray, heavy_indices: np.ndarray,
                           *, number_of_real_scalars: int) -> tuple[np.ndarray,np.ndarray]:
    """Extract LLHH and LLSS from FULL A evolved psi^2 phi^2 tensor.

    All other Wilson components remain present in the input artifact, but
    their potentially required *additional* matching is separately flagged.
    """
    C=finite(full_C,'full evolved Wilson C_ijab')
    if C.ndim!=4 or C.shape[0]!=C.shape[1] or C.shape[2:]!=(number_of_real_scalars,)*2:
        raise ValueError('Full Wilson tensor must be [nWeyl,nWeyl,nRealScalar,nRealScalar]')
    li=np.asarray(lepton_indices)
    h=np.asarray(higgs_indices)
    s=np.asarray(heavy_indices)
    if li.shape!=(3,2) or not np.issubdtype(li.dtype,np.integer):
        raise ValueError('Explicit 3x2 zero-based Weyl lepton indices required')
    if h.shape!=(4,) or s.ndim!=1 or h.dtype.kind not in 'iu' or s.dtype.kind not in 'iu':
        raise ValueError('Explicit zero-based Higgs[4] and heavy scalar indices required')
    if (len(set(li.flat))!=6 or len(set(h))!=4 or len(set(s))!=len(s) or
        set(h).intersection(s) or len(s)+4!=number_of_real_scalars):
        raise ValueError('Duplicate/missing scalar or lepton indices')
    if (min(li.flat)<0 or max(li.flat)>=C.shape[0] or
        min(np.r_[h,s])<0 or max(np.r_[h,s])>=number_of_real_scalars):
        raise ValueError('Index out of full Wilson dimension')
    l=li.reshape(-1)
    def part(indices):
        out=C[np.ix_(l,l,indices,indices)]
        return out.reshape(3,2,3,2,len(indices),len(indices)).transpose(0,2,1,3,4,5)
    hh,ss=part(h),part(s)
    for name,a in (('LLHH',hh),('LLSS',ss)):
        scale=max(float(np.max(abs(a))),1e-40)
        if np.max(abs(a-a.transpose(1,0,3,2,4,5)))>1e-8*scale:
            raise ValueError(f'{name} Weyl flavour/weak exchange is not symmetric')
        if np.max(abs(a-a.swapaxes(4,5)))>1e-8*scale:
            raise ValueError(f'{name} real scalar legs are not symmetric')
    return hh,ss


def rotate_llss(c: np.ndarray, R: np.ndarray) -> np.ndarray:
    return np.einsum('fgijab,ak,bl->fgijkl',c,R,R,optimize=True)


def bubble_amplitude(c: np.ndarray, q: np.ndarray, mass2: np.ndarray,
                     mu_s: float, *, detailed: bool=True) -> dict:
    """One simultaneous eigenbasis sum; off-diagonal eigenstate pairs counted twice."""
    c=finite(c,'LLSS')
    n=c.shape[-1]
    if c.shape!=(3,3,2,2,n,n):raise ValueError('LLSS must be [3,3,2,2,n,n]')
    q=finite(q,'HHSS quartic',real=True)
    if q.shape!=(4,4,n,n):raise ValueError('HHSS quartic must be [4,4,n,n]')
    if np.max(abs(q-q.swapaxes(0,1)))>1e-9*max(1e-40,np.max(abs(q))) or not np.allclose(q,q.swapaxes(2,3)):
        raise ValueError('Quartic HH/SS indices must be separately symmetric')
    m,U=mass_eigenbasis(mass2)
    cm=rotate_llss(c,U)
    qm=np.einsum('ABab,ak,bl->ABkl',q,U,U,optimize=True)
    B=np.array([[b0bar(m[a],m[b],mu_s) for b in range(n)] for a in range(n)])
    total=-np.einsum('fgijab,ABab,ab->fgijAB',cm,qm,B,optimize=True)/(2*LOOP)
    pair=[]
    if detailed:
        for a in range(n):
            for b in range(a,n):
                factor=1 if a==b else 2
                amp=-factor*np.einsum('fgij,AB->fgijAB',cm[...,a,b],qm[...,a,b],optimize=True)*B[a,b]/(2*LOOP)
                value,residual=project_master(amp,enforce_pure=False)
                pair.append({'a':a,'b':b,'multiplicity':factor,'m_a_gev':float(np.sqrt(m[a])),
                             'm_b_gev':float(np.sqrt(m[b])),'b0bar':float(B[a,b]),
                             'projected_c5':value,'projection_residual':residual})
    return {'amplitude':total,'masses2':m,'eigenvectors':U,'pair_data':pair,
            'all_eigenstates_removed_at_one_scale':True,'pair_count':n*(n+1)//2}


def match_scalar_threshold(*, full_wilson_ms: np.ndarray,
                            lepton_indices: np.ndarray,higgs_indices: np.ndarray,
                            heavy_indices: np.ndarray,
                            full_quartic_ms: np.ndarray,
                            scalar_mass2_background: np.ndarray,
                            R_wilson_from_background: np.ndarray,
                            mu_s: float, mode: str,
                            full_wilson_tree_mf: np.ndarray | None=None,
                            full_wilson_strict_ms: np.ndarray | None=None,
                            separately_supplied_c5_ms: np.ndarray | None=None,
                            independent_hard_c5: np.ndarray | None=None,
                            full_quartic_channels: np.ndarray | None=None,
                            channel_names: tuple[str,...]=(),
                            projection_tolerance: float=1e-8,
                            enforce_pure_weinberg: bool=True) -> dict:
    """Transform A Wilson into B real mass background, match all scalars at once.

    STRICT requires BOTH the MF order-zero full Wilson tensor (bubble)
    and the explicitly one-loop-truncated full Wilson at MS (C5).
    The RGE-improved mode uses the fully evolved MS tensor in both places.
    Running scalar couplings supplied at MS remain unexpanded background
    parameters, so one-loop strictness of their trajectories is not asserted.
    """
    if mode not in ('strict_one_loop','rge_improved'):
        raise ValueError('Unknown counting mode')
    ms=float(mu_s)
    if not np.isfinite(ms) or ms<=0:raise ValueError('Invalid common scalar matching scale')
    n=np.asarray(scalar_mass2_background).shape[0]
    R=check_rotation(R_wilson_from_background,n)
    mass2=finite(scalar_mass2_background,'mass2',real=True)
    if mass2.shape!=(n,n):raise ValueError('Mass2 must be nxn')
    eig,_=mass_eigenbasis(mass2)
    fullq=check_quartic(full_quartic_ms,n)
    q=fullq[:4,:4,4:,4:]
    hh,ss=extract_lepton_sectors(full_wilson_ms,lepton_indices,higgs_indices,heavy_indices,
                                 number_of_real_scalars=n+4)
    c5_improved,r0=project_master(hh,tolerance=projection_tolerance,
                                  enforce_pure=enforce_pure_weinberg)
    c5_eft1=c5_improved
    c5_fixed_order_status='rge_improved_resummed_endpoint'
    if separately_supplied_c5_ms is not None:
        supplied=symmetric3(separately_supplied_c5_ms,'separate C5_EFT1')
        err=np.linalg.norm(supplied-c5_improved)/max(np.linalg.norm(supplied),np.linalg.norm(c5_improved),1e-40)
        if err>projection_tolerance:
            raise ValueError(f'Separate evolved C5 disagrees with full Wilson at MS: {err:.4g}')
    if mode=='strict_one_loop':
        if full_wilson_tree_mf is None:
            raise ValueError('Strict one-loop mode requires full order-zero C(MF) from A')
        if full_wilson_strict_ms is None:
            raise ValueError('Strict one-loop requires A full one-loop-TRUNCATED EFT1 Wilson at MS; evolved resummed C5 cannot substitute')
        hh_strict,ss_strict=extract_lepton_sectors(
            full_wilson_strict_ms,lepton_indices,higgs_indices,heavy_indices,
            number_of_real_scalars=n+4)
        c5_eft1,r_strict=project_master(hh_strict,tolerance=projection_tolerance,
                                      enforce_pure=enforce_pure_weinberg)
        _,ss_use=extract_lepton_sectors(full_wilson_tree_mf,lepton_indices,higgs_indices,heavy_indices,
                                       number_of_real_scalars=n+4)
        kernel_source='tree_order0_at_MF'
        c5_fixed_order_status='A_full_one_loop_truncated_EFT1_at_MS'
    else:
        ss_use=ss
        kernel_source='evolved_full_at_MS_rge_improved'
    ss_b=rotate_llss(ss_use,R)
    bubble=bubble_amplitude(ss_b,q,mass2,ms)
    c5_bubble,r1=project_master(bubble['amplitude'],tolerance=projection_tolerance,
                                 enforce_pure=enforce_pure_weinberg)
    hard=(np.zeros((3,3),complex) if independent_hard_c5 is None else
          symmetric3(independent_hard_c5,'independent finite C5 hard term'))
    total=symmetric3(c5_eft1+c5_bubble+hard,'C5_SMEFT_at_MS')
    channels={}
    if full_quartic_channels is not None:
        arr=finite(full_quartic_channels,'quartic channels',real=True)
        if arr.ndim!=5 or arr.shape[1:]!=fullq.shape or len(arr)!=len(channel_names):
            raise ValueError('Full quartic channels must be [k,4+n,4+n,4+n,4+n] with names')
        if len(set(channel_names))!=len(channel_names):
            raise ValueError('Quartic channel names must be unique')
        if not np.allclose(arr.sum(axis=0),fullq,rtol=1e-9,atol=1e-12*max(1e-40,float(np.max(abs(fullq))))):
            raise ValueError('Quartic channels do not sum to full quartic')
        for name,Q in zip(channel_names,arr):
            amp=bubble_amplitude(ss_b,Q[:4,:4,4:,4:],mass2,ms,detailed=False)
            channels[name],_=project_master(amp['amplitude'],enforce_pure=False)
        if not np.allclose(sum(channels.values(),np.zeros((3,3),complex)),c5_bubble,
                           rtol=1e-8,atol=1e-26):
            raise ValueError('Per-channel scalar bubble does not sum to total')
    return {'c5_eft1':c5_eft1,'c5_eft1_rge_improved':c5_improved,
            'c5_eft1_power_counting_source':c5_fixed_order_status,
            'strict_c5_projection_residual':r_strict if mode=='strict_one_loop' else None,
            'c5_scalar_bubble':c5_bubble,'c5_independent_hard':hard,
            'c5_smeft_at_ms':total,'rge_input_C5_residual':r0,
            'scalar_bubble_projection_residual':r1,'bubble':bubble,
            'channel_c5':channels,'counting_mode':mode,
            'bubble_Wilson_order_source':kernel_source,
            'masses2_gev2':eig,'scalar_components_removed':n,
            'single_scalar_threshold_gev':ms,
            'no_direct_leading_log_reinserted':True,
            'scheme':'MSbar',
            'bubble_physics_status':'conditional_quartic_pair_only_not_full_finite_one_loop'}
