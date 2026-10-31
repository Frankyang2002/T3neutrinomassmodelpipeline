"""Grouped (X,Y) scalar-pair bubble diagnostic in the *single* EFT1 -> SM step.

Only the subset with one insertion of C[LL phi phi] and one quartic
V=lambda_abcd phi_a phi_b phi_c phi_d/4! is calculated. Both massive
scalar eigenfields are integrated at a COMMON chosen matching scale.
No field-equation, wavefunction, other one-loop topology, or independent
UV/EFT counterterm calculation is claimed. NOT a physical C5 prediction.
"""
from __future__ import annotations
import math
from typing import Any
import numpy as np
from Numerical.experimental.T3BCScalarBubbleThreshold import b0bar_zero

LOOP=16.0*math.pi**2
SHAPE_C=(3,3,12,12,6)
SHAPE_Q=(12,12,12,12)


def require(ok: bool, reason: str) -> None:
    if not bool(ok):raise ValueError('Grouped B/C finite bubble: '+reason)


def group_bubble(C: np.ndarray, q: np.ndarray, masses_sq: tuple[float,float],
                 mu_gev: float) -> tuple[np.ndarray,dict[str,Any]]:
    """One-loop HH-block scalar-pair bubble involving >=1 heavy internal leg.

    For c,d in H,X,Y, internal masses squared are 0,mX^2,mY^2.
    The factor 1/(2*16pi^2) and B0bar are exactly those in the project's
    existing Y-only bubble diagnostic, generalized to group X and Y.
    """
    C=np.asarray(C,dtype=np.complex128)
    q=np.asarray(q,dtype=np.complex128)
    require(C.shape==SHAPE_C and q.shape==SHAPE_Q,'unexpected tensor shapes')
    require(np.isfinite(C).all() and np.isfinite(q).all(),'non-finite source tensors')
    require(np.linalg.norm(C-C.swapaxes(0,1))<1e-8*max(np.linalg.norm(C),1e-25),
            'lepton/gauge pair asymmetry')
    require(np.linalg.norm(C-C.swapaxes(2,3))<1e-8*max(np.linalg.norm(C),1e-25),
            'scalar pair asymmetry')
    require(np.linalg.norm(q-q.transpose(1,0,2,3))<1e-8*max(np.linalg.norm(q),1e-25),
            'quartic pair symmetry broken')
    require(mu_gev>0 and np.isfinite(mu_gev),'invalid grouped matching scale')
    m=np.asarray(masses_sq,float)
    require(m.shape==(2,) and np.isfinite(m).all() and np.all(m>0),'positive X/Y squared masses required')
    mfull=np.r_[np.zeros(4),np.full(4,m[0]),np.full(4,m[1])]
    delta=np.zeros((3,3,4,4,6),complex)
    sector={}
    weighted={}
    for first,name1 in ((slice(0,4),'H'),(slice(4,8),'X'),(slice(8,12),'Y')):
        for second,name2 in ((slice(0,4),'H'),(slice(4,8),'X'),(slice(8,12),'Y')):
            if name1==name2=='H':continue
            a=first.start
            b=second.start
            b0=b0bar_zero(mfull[a],mfull[b],mu_gev**2)
            term=np.einsum('abcd,ijcdf->ijabf',q[:4,:4,first,second],
                           C[:,:,first,second,:],optimize=True)
            delta-=b0*term/(2*LOOP)
            key=name1+name2
            sector[key]=float(np.linalg.norm(term))
            weighted[key]=float(b0)
    require(np.isfinite(delta).all(),'nonfinite grouped correction')
    require(np.linalg.norm(delta-delta.swapaxes(0,1)) <1e-7*max(np.linalg.norm(delta),1e-25),
            'grouped correction fermion symmetry broken')
    require(np.linalg.norm(delta-delta.swapaxes(2,3)) <1e-7*max(np.linalg.norm(delta),1e-25),
            'grouped correction scalar symmetry broken')
    return delta,{'internal_pairs':sector,'B0bar_by_ordered_pair':weighted,
                  'relative_operator_scope':'one insertion LLphi phi times quartic, all X/Y heavy internal pairs',
                  'finite_delta_HH_norm':float(np.linalg.norm(delta)),
                  'physical_finite_matching_certified':False}


def rge_log_check(C:np.ndarray,q:np.ndarray,masses_sq:tuple[float,float],
                  mu_gev:float,h:float=1e-4)->dict:
    """Check d(delta C)/d ln(mu)= -(heavy quartic beta source)/(16pi^2)."""
    high,_=group_bubble(C,q,masses_sq,mu_gev*math.exp(h))
    low,_=group_bubble(C,q,masses_sq,mu_gev*math.exp(-h))
    numerical=(high-low)/(2*h)
    analytic=np.zeros_like(numerical)
    for c in range(12):
        for d in range(12):
            if c<4 and d<4:continue
            analytic-=np.einsum('ab,ijf->ijabf',q[:4,:4,c,d],C[:,:,c,d,:],optimize=True)/LOOP
    residual=float(np.linalg.norm(numerical-analytic)/max(np.linalg.norm(analytic),1e-30))
    require(residual<5e-6,'scalar bubble UV logarithm failed beta cancellation')
    return {'relative_log_derivative_residual':residual,
            'log_derivative_norm':float(np.linalg.norm(analytic)),
            'counterterm_certification':False}


def project_delta(delta:np.ndarray,template:dict,ids:np.ndarray):
    """Embed the HH shift in the existing full tensor projector conventions."""
    from Numerical.experimental import T3BCCoupledWilsonTensor as tensor
    full=np.zeros(SHAPE_C,complex)
    full[:,:,:4,:4,:]=delta
    six,res=tensor.weinberg_projection(full,template,ids,threshold=5e-5)
    require(res<5e-5,'grouped bubble not compatible with Weinberg HH tensor')
    return tensor.flavor_to_matrix(six),float(res)
