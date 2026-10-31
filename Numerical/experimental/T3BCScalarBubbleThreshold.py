"""Experimental scalar-pair one-loop contribution to B/C Y/X matching.

This module computes only the single insertion (LL phi phi) x (phi^4)
zero-momentum scalar bubble. It is NOT the complete one-loop matching of the
full theory, especially when light/heavy scalar scales are nearby.

Conventions: quartic V = lambda_abcd phi_a phi_b phi_c phi_d/4!, and
16*pi^2 beta(C_ijab) contains +lambda_abcd C_ijcd (project tensor RGE).
Renormalised scalar bubble is B0bar(0;m_c^2,m_d^2;mu^2), with
   d B0bar / d ln(mu) = 2.
Thus the scalar-pair contribution to matching is
   delta C_ijab = -1/(2*16*pi^2) sum_cd lambda_abcd C_ijcd B0bar_cd.
The derivative of delta cancels the heavy-pair contribution to beta(C),
fixing the logarithm coefficient; its finite term still needs independent
native matching validation before physical certification.

Original real-field order (H[4],X[4],Y[4]); 6 physical flavour pairs.
"""
from __future__ import annotations
import math
import numpy as np

LOOP=16.0*math.pi**2


def require(ok, msg):
    if not ok:
        raise ValueError('T3 B/C scalar bubble: '+msg)


def b0bar_zero(m2_a, m2_b, mu2):
    """MSbar finite zero-momentum 2-propagator scalar integral (dimensionless).

    Includes stable equal-mass and massless limiting cases. No tadpoles,
    momentum corrections, scheme conversions or evanescent matching.
    """
    a=float(m2_a); b=float(m2_b); u=float(mu2)
    require(np.isfinite([a,b,u]).all() and a>=0 and b>=0 and a+b>0 and u>0,
            'invalid mass squares')
    if abs(a-b)<=1e-7*max(a,b):
        # Small-splitting Taylor to avoid catastrophic subtractive cancellation.
        center=(a+b)/2.0
        d=(a-b)/(a+b)
        return float(-math.log(center/u) + d*d/6.0 + d**4/20.0)
    def f(x):
        return 0.0 if x==0.0 else x*math.log(x/u)
    return float(1.0-(f(a)-f(b))/(a-b))


def partition_beta(Cmass, qmass):
    """All scalar-quartic heavy-leg contributions to the 8-real-field EFT2 beta.

    Retains full 3x3 Wilson fermion gauge indices and six flavour components.
    Y-only and Y-X are kept separate, so the finite-matching approximation
    (and the hard mass expansion) can be independently inspected.
    """
    C=np.asarray(Cmass,complex); q=np.asarray(qmass,float)
    require(C.shape==(3,3,12,12,6) and q.shape==(12,)*4,
            'incorrect mass basis dimensions')
    require(np.isfinite(C).all() and np.isfinite(q).all(), 'nonfinite inputs')
    require(np.linalg.norm(C-C.swapaxes(2,3))<=1e-8*max(1e-25,np.linalg.norm(C)),
            'Wilson scalar symmetry broken')
    require(np.linalg.norm(q-q.transpose(0,1,3,2))<=1e-8*max(1e-25,np.linalg.norm(q)),
            'quartic symmetry broken')
    def block(p, r):
        return np.einsum('abcd,ijcdf->ijabf',q[:8,:8,p,r],C[:,:,p,r,:],optimize=True)
    hh=block(slice(0,4),slice(8,12))+block(slice(8,12),slice(0,4))
    xx=block(slice(4,8),slice(8,12))+block(slice(8,12),slice(4,8))
    yy=block(slice(8,12),slice(8,12))
    by_direct=np.einsum('abcd,ijcdf->ijabf',q[:8,:8,:,8:],C[:,:,:,8:,:],optimize=True)
    by_direct+=np.einsum('abcd,ijcdf->ijabf',q[:8,:8,8:,:],C[:,:,8:,:,:],optimize=True)
    by_direct-=yy
    require(np.allclose(by_direct,hh+xx+yy,rtol=4e-12,atol=1e-24),
            'heavy-leg beta partition not exhaustive')
    return {'HY':hh,'XY':xx,'YY':yy}


def bubble_boundary(Cmass,qmass,mu_y,mx2,*,mode):
    """Return candidate finite one-loop Y-threshold tensor and diagnostics.

    mode='full-mass': keep the exact MX/MY dependence inside B0bar(0).
    mode='hard-massless-X': set m_X^2=0 inside mixed loop. The two are NOT
    formally the same Wilson matching scheme: the full-mass number retains
    unsuppressed power corrections/higher-dimensional effects for close masses.
    """
    require(mode in ('full-mass','hard-massless-X'),'unknown bubble scheme')
    my=float(mu_y); mx=float(mx2)
    require(my>0 and 0<=mx<my*my,'invalid Y/X hierarchy')
    sect=partition_beta(Cmass,qmass)
    x_arg=mx if mode=='full-mass' else 0.
    b_hy=b0bar_zero(my*my,0.,my*my)
    b_xy=b0bar_zero(my*my,x_arg,my*my)
    b_yy=b0bar_zero(my*my,my*my,my*my)
    delta=-(b_hy*sect['HY']+b_xy*sect['XY']+b_yy*sect['YY'])/(2*LOOP)
    require(delta.shape==(3,3,8,8,6) and np.isfinite(delta).all(),
            'nonfinite threshold correction')
    require(np.allclose(delta,delta.swapaxes(0,1),rtol=1e-8,atol=1e-22),
            'fermion symmetry broken by correction')
    require(np.allclose(delta,delta.swapaxes(2,3),rtol=1e-8,atol=1e-22),
            'scalar symmetry broken by correction')
    return delta,{
        'scheme':mode,'B0_YH':b_hy,'B0_YX':b_xy,'B0_YY':b_yy,
        'beta_quartic_YH_norm':float(np.linalg.norm(sect['HY'])),
        'beta_quartic_YX_norm':float(np.linalg.norm(sect['XY'])),
        'beta_quartic_YY_norm':float(np.linalg.norm(sect['YY'])),
        'delta_wilson_frobenius_norm':float(np.linalg.norm(delta)),
        'matching_sector':'one Wilson insertion times one scalar quartic, zero momentum',
        'complete_finite_matching_certified':False,
    }


def x_equal_mass_bubble_at_threshold(m_x):
    """ONLY the equal-mass XX quartic bubble in MSbar at mu=mX.

    Other contributions to the X threshold are not covered by this identity.
    """
    m=float(m_x)
    require(m>0 and math.isfinite(m),'invalid X mass')
    return b0bar_zero(m*m,m*m,m*m)
