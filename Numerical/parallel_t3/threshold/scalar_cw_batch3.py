"""One-loop real-heavy-scalar Coleman--Weinberg Higgs-potential threshold.

The heavy scalars are integrated out TOGETHER at muS in the unbroken Z2 EFT1.
For real heavy fields S_a and four real Higgs fields h_i use exactly
 V2 = 1/2 S^T M2 S,
 V4 = lambda_ijab h_i h_j S_a S_b / 4 (lambda is fourth derivative),
 M2_S(h) = M2 + (1/2) lambda_ijab h_i h_j.

At one loop, in MSbar, the scalar part of the potential is
 V1 = Tr[M2_S(h)^2*(log(M2_S(h)/muS^2)-3/2)]/(64*pi^2).
The h^2 and h^4 Taylor coefficients are evaluated in the MASS EIGENBASIS;
all nondegenerate eigenstate pairs are retained. Not a C5 matching diagram.
All inputs are the saved strict Batch2 full real-scalar quartic and mass tensors.
"""
from __future__ import annotations

from itertools import permutations
from math import log, log1p, pi
import numpy as np

FORMAT='t3-scalar-cw-Higgs-potential-MSbar-common-MS-batch3-v1'


def _finite_real(a, name, *, shape=None):
    x=np.asarray(a)
    if shape is not None and x.shape!=shape:
        raise ValueError(f'{name}: shape {x.shape} expected {shape}')
    if not np.isfinite(x).all() or np.max(abs(x.imag))>1e-12*max(1.,float(np.max(abs(x.real)))):
        raise ValueError(f'{name}: nonfinite/non-real input')
    return np.asarray(x.real,dtype=float)


def fprime(m2:float,mu2:float)->float:
    if not np.isfinite([m2,mu2]).all() or m2<=0 or mu2<=0:
        raise ValueError('Positive finite squared masses/scales required')
    return 2*m2*(log(m2/mu2)-1.)


def divided_fprime(a:float,b:float,mu2:float)->float:
    """Stable divided difference of f'(x)=2x(log(x/mu2)-1)."""
    if not np.isfinite([a,b,mu2]).all() or min(a,b,mu2)<=0:
        raise ValueError('Positive finite squared masses/scales required')
    m=(a+b)/2
    d=(a-b)/(a+b)
    if abs(d)<1e-3:
        d2=d*d
        # ((1+d)ln(1+d)-(1-d)ln(1-d))/(2d)-1
        corr=-(d2/6+d2*d2/20+d2**3/42+d2**4/72)
    else:
        corr=((1+d)*log1p(d)-(1-d)*log1p(-d))/(2*d)-1.
    return 2*(log(m/mu2)+corr)


def _sym4(x:np.ndarray)->np.ndarray:
    """Symmetrise last four indices; channel and channel-pair axes unchanged."""
    return sum((x.transpose(tuple(range(x.ndim-4))+tuple(x.ndim-4+i for i in p))
                for p in permutations(range(4))),np.zeros_like(x))/24.


def higgs_quartic_reference()->np.ndarray:
    """Source real Higgs convention: V4=lambdaH*(H^dag H)^2/2.

    H=(h0+i h1,h2+i h3)/sqrt(2), so lambda4^H is the sum of
    Kronecker-delta pairings. Do not confuse with lambdaH*(H^dag H)^2.
    """
    eye=np.eye(4)
    return (np.einsum('ij,kl->ijkl',eye,eye)+
            np.einsum('ik,jl->ijkl',eye,eye)+
            np.einsum('il,jk->ijkl',eye,eye))


def one_loop_higgs_potential(mass_squared:np.ndarray, full_quartic:np.ndarray,
                             mu_s:float, *, channels:np.ndarray|None=None,
                             names:tuple[str,...]|None=None) -> dict:
    """Compute delta(mH^2)_ij and delta(lambda4)_ijkl, including all channels.

    Output is conditional on the FULL original EFT1 potential's source/Weyl
    frame identification. The source no-factorial fourth derivative is assumed.
    No SM gauge matching, C5 matching or external wavefunction correction
    is invented here. Those require distinct one-loop diagrams.
    """
    m=_finite_real(mass_squared,'heavy Hessian')
    if m.ndim!=2 or m.shape[0]!=m.shape[1] or m.shape[0]<2:
        raise ValueError('Require complete two-multiplet real heavy Hessian')
    n=len(m)
    q=_finite_real(full_quartic,'full real scalar quartic',shape=(n+4,)*4)
    if not np.allclose(m,m.T,atol=1e-7,rtol=1e-11):
        raise ValueError('Heavy Hessian is not symmetric')
    if not np.allclose(q,q.transpose(1,0,2,3),atol=1e-10,rtol=1e-10) or not np.allclose(q,q.transpose(0,2,1,3),atol=1e-10,rtol=1e-10) or not np.allclose(q,q.transpose(0,1,3,2),atol=1e-10,rtol=1e-10):
        raise ValueError('Quartic is not fully symmetric')
    if not (np.isfinite(mu_s) and mu_s>0):
        raise ValueError('Invalid common scalar scale')
    eig,U=np.linalg.eigh(m)
    if eig[0]<=0:
        raise ValueError('Tachyonic heavy scalar eigenstate')
    if channels is None:
        cq=q[np.newaxis]
        labels=('total_full_quartic',)
    else:
        cq=_finite_real(channels,'quartic channels')
        if cq.ndim!=5 or cq.shape[1:]!=q.shape:
            raise ValueError('Quartic channel tensor dimensions incorrect')
        labels=tuple(names or ())
        if len(labels)!=len(cq) or len(set(labels))!=len(labels) or not labels:
            raise ValueError('Unique quartic names required')
        if not np.allclose(cq.sum(axis=0),q,atol=1e-11,rtol=1e-10):
            raise ValueError('Quartic channels do not sum to the source quartic')
    # At fixed heavy background, only HHSS vertices contribute to h^2/h^4.
    hhss=cq[:,:4,:4,4:,4:]
    rotated=np.einsum('rijab,ac,bd->rijcd',hhss,U,U,optimize=True)
    if not np.allclose(rotated,rotated.swapaxes(-1,-2),atol=1e-9,rtol=1e-10):
        raise ValueError('Mass-eigenbasis HHSS coupling is not real symmetric')
    fp=np.array([fprime(float(x),mu_s*mu_s) for x in eig],float)
    dd=np.array([[divided_fprime(float(x),float(y),mu_s*mu_s) for y in eig] for x in eig],float)
    # Trace expansion: V_1^(2) = [sum_a f'(m2_a)*Q_ijaa] h_i h_j /(128pi²).
    # Hence delta mH²_ij = [sum_a f'(m2_a)*Q_ijaa] /(64pi²).
    md=np.einsum('a,rijaa->rij',fp,rotated,optimize=True)/(64*pi*pi)
    # V_1^(4) = sum_ab f'[ma2,mb2]*Q_ijab*Q_klba h_i h_j h_k h_l/(512pi²).
    # lambda4 is the FULL fourth derivative, i.e. 24 times its sym coefficient.
    raw=np.einsum('ab,rijab,sklab->rsijkl',dd,rotated,rotated,optimize=True)
    qd=_sym4(raw)*(24/(512*pi*pi))
    md_total=md.sum(axis=0)
    qd_total=qd.sum(axis=(0,1))
    ref=higgs_quartic_reference()
    dlam=float(np.vdot(ref,qd_total).real/np.vdot(ref,ref).real)
    lres=float(np.linalg.norm(qd_total-dlam*ref)/max(np.linalg.norm(qd_total),1e-100))
    mbar=float(np.trace(md_total)/4)
    mres=float(np.linalg.norm(md_total-mbar*np.eye(4))/max(np.linalg.norm(md_total),1e-100))
    channel_lambdas=np.einsum('ijkl,rsijkl->rs',ref,qd,optimize=True)/np.vdot(ref,ref)
    if not np.allclose(channel_lambdas.sum(),dlam,rtol=1e-10,atol=1e-14):
        raise AssertionError('Source-resolved one-loop quartic shift does not close')
    if not np.allclose(md_total,md.sum(axis=0),rtol=1e-10,atol=1e-14):
        raise AssertionError('Source-resolved mass threshold does not close')
    return {'format':FORMAT,'mass2_gev2':eig,'mass_eigenvectors':U,
            'higgs_mass2_channel_shifts_gev2':md,
            'higgs_mass2_shift_gev2':md_total,
            'higgs_quartic_channel_pairs':qd,
            'higgs_quartic_shift_tensor':qd_total,
            'delta_lambdaH_source_convention':dlam,
            'delta_mH_squared_source_convention_gev2':mbar,
            'quartic_projected_channel_pairs':channel_lambdas,
            'quartic_singlet_residual':lres,
            'quadratic_singlet_residual':mres,
            'channel_names':labels,
            'num_heavy_real_modes':n,
            'common_scalar_threshold_gev':float(mu_s),
            'scheme':'MSbar','unbroken_phase':True,
            'physical_source_normalization_certified':False,
            'C5_one_loop_matching_changed':False}
