"""Independent B/C scalar-mixing propagator and one-loop mass-kernel checks.

The tensor K_(ab),(cd) is defined by -MF^2 int_0^infty dt t/(t+MF^2)
 G_ac(t) G_bd(t), with G(t)=(t*I + M_scalar^2)^(-1).
It contains *only* mass dependence, not UV spinor/vertex normalisation,
external-state combinatorics, full Wilson matching, or EFT counterterms.
No physical C5 is computed. 16*pi^2 factors intentionally not applied.
"""
from __future__ import annotations

import math
from typing import Iterable

import mpmath as mp
import numpy as np
from scipy.integrate import quad_vec


def require(condition: bool, why: str) -> None:
    if not condition:
        raise ValueError('B/C mixed kernel: ' + why)


def _pos(x: float, label: str) -> float:
    v = float(x)
    require(math.isfinite(v) and v > 0, label + ' must be finite and positive')
    return v


def positive_mass_matrix(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    A = np.asarray(matrix, dtype=complex)
    require(A.shape == (2, 2) and np.isfinite(A).all(), 'scalar A must be finite 2x2')
    require(np.linalg.norm(A-A.conj().T) < 1e-12 * max(1.,np.linalg.norm(A)),
            'scalar mass matrix must be Hermitian')
    eig, U = np.linalg.eigh(A)
    require(np.all(eig > 0), 'scalar Hessian not positive definite')
    require(abs(eig[1]-eig[0]) > 1e-10*max(eig), 'degenerate mass basis')
    return eig, U


def signed_exact_master(mf_sq: float, ma_sq: float, mb_sq: float) -> float:
    """Full three-mass, zero-momentum Euclidean loop *mass* function.

    -MF^2 int_0^infty t dt / ((t+MF^2)(t+ma^2)(t+mb^2)).
    Uses a 55-digit second divided difference of f(x)=x log x.
    Equality limits are handled analytically, avoiding cancellation.
    """
    A = _pos(mf_sq, 'MF^2'); B = _pos(ma_sq, 'ma^2'); C = _pos(mb_sq, 'mb^2')
    with mp.workdps(55):
        a, b, c = map(lambda x: mp.mpf(str(x)), (A,B,C))
        scale = max(a,b,c)
        # Treat true equal masses only; near-equal masses handled by high precision.
        if a == b == c:
            d = 1/(2*a)
        elif b == c:
            d = (a*mp.log(a/b) - (a-b))/(a-b)**2 if a != b else 1/(2*a)
        elif a == c:
            d = (b*mp.log(b/a) - (b-a))/(b-a)**2 if a != b else 1/(2*a)
        elif a == b:
            d = (c*mp.log(c/a) - (c-a))/(c-a)**2 if a != c else 1/(2*a)
        else:
            # Common log reference cancels exactly, improving conditioning.
            vals = (a,b,c)
            d = sum(x*mp.log(x/scale) /
                    ((x-vals[(i+1)%3])*(x-vals[(i+2)%3]))
                    for i,x in enumerate(vals))
        ans = float(-a*d)
    require(math.isfinite(ans) and ans < 0, 'signed triangle master is not finite negative')
    return ans


def strict_heavy_mass_master(mf_sq:float, ma_sq:float, mb_sq:float) -> float:
    """Leading 1/MF hard+soft finite function from Workstream C v19."""
    A=_pos(mf_sq,'MF^2');B=_pos(ma_sq,'ma^2');C=_pos(mb_sq,'mb^2')
    with mp.workdps(55):
        a,b,c=map(lambda x:mp.mpf(str(x)),(A,B,C))
        if b==c:
            ans=1+mp.log(b/a)
        else:
            ans=(b*mp.log(b/a)-c*mp.log(c/a))/(b-c)
        return float(ans)


def hard_soft_master(mf_sq:float,ma_sq:float,mb_sq:float,scale_sq:float) -> dict:
    A=_pos(mf_sq,'MF^2'); B=_pos(ma_sq,'ma^2');C=_pos(mb_sq,'mb^2')
    Q=_pos(scale_sq,'renormalization scale squared')
    with mp.workdps(55):
        a,b,c,q=map(lambda x:mp.mpf(str(x)),(A,B,C,Q))
        h=1+mp.log(q/a)
        b0=1-((b*mp.log(b/q)-c*mp.log(c/q))/(b-c)
              if b!=c else (1+mp.log(b/q)))
        soft=-b0
        total=h+soft
        ref=strict_heavy_mass_master(A,B,C)
        residual=abs(float(total)-ref)/max(abs(ref),1.)
        require(residual<2e-12,'hard+soft subtraction identity failed')
        return {'hard':float(h),'soft':float(soft),'sum':float(total),
                'strict':ref,'relative_scale_cancellation_error':residual}


def spectral_kernel(A: np.ndarray, mf_sq: float) -> np.ndarray:
    """4x4 matrix in ordered species-pair coordinates (00,01,10,11)."""
    eig, U = positive_mass_matrix(A)
    mf=_pos(mf_sq,'MF^2')
    phi=np.array([signed_exact_master(mf,eig[i],eig[j])
                  for i in range(2) for j in range(2)],float)
    V=np.kron(U,U)
    K=(V*phi[np.newaxis,:])@V.conj().T
    require(np.linalg.norm(K-K.conj().T)<1e-11*max(np.linalg.norm(K),1.),
            'spectral kernel non-Hermitian')
    return K


def _resolvent(A: np.ndarray, t:float, order:int|None) -> np.ndarray:
    if order is None:
        return np.linalg.inv(A+t*np.eye(2))
    diag=np.diag(A).real
    D=np.diag(1/(t+diag))
    off=A-np.diag(np.diag(A))
    X=-D@off
    total=np.asarray(D,dtype=complex).copy()
    term=np.asarray(D,dtype=complex).copy()
    for _ in range(order):
        term=X@term
        total+=term
    return total


def euclidean_matrix_kernel(A: np.ndarray, mf_sq: float,
                            orders: Iterable[int]=(0,1,2,4,8),
                            *, epsrel:float=2e-10) -> dict:
    """Integrate exact 2x2 G⊗G, and selected Dyson insertion truncations.

    Numerical integration is independently over real log Euclidean momentum t;
    the reference spectral result uses 55-digit analytic masses/eigenvectors.
    """
    positive_mass_matrix(A)
    mass=_pos(mf_sq,'MF^2')
    orders=tuple(sorted(set(int(x) for x in orders)))
    require(orders and all(0<=n<=24 for n in orders), 'invalid mass-insertion orders')
    scale=min(float(np.linalg.eigvalsh(A)[0]),mass)
    lo=math.log(scale)-20.
    hi=math.log(max(float(np.linalg.eigvalsh(A)[1]),mass))+27.
    channels=(None,*orders)
    def integrand(x:float) -> np.ndarray:
        t=math.exp(x)
        pref=-mass*t*t/(t+mass)
        return np.concatenate([(pref*np.kron(_resolvent(A,t,n),_resolvent(A,t,n))).ravel()
                               for n in channels])
    value,err=quad_vec(integrand,lo,hi,epsabs=1e-11,epsrel=epsrel,limit=1000)
    kblocks=value.reshape(len(channels),4,4)
    analytic=spectral_kernel(A,mass)
    norm=max(float(np.linalg.norm(analytic)),1e-30)
    relative=float(np.linalg.norm(kblocks[0]-analytic)/norm)
    # 2e-8 covers very mild finite endpoint tail and roundoff.
    require(relative<3e-8,
            f'original-basis resolvent disagrees with exact mass-eigenstate spectral integral: {relative:.3g}')
    ins={str(n):float(np.linalg.norm(kblocks[i+1]-analytic)/norm)
         for i,n in enumerate(orders)}
    return {'mf_gev':math.sqrt(mass),'original_basis_vs_spectral_relative':relative,
            'quad_vec_error_estimate':float(err),
            'dyson_relative_errors':ins,
            'exact_kernel_4x4':{'real':analytic.real.tolist(),'imag':analytic.imag.tolist()},
            'mass_insertion_series_convergent_test':bool(ins[str(orders[-1])]<3e-8)}


def invariant_quartic_contraction(A:np.ndarray,mf_sq:float,
                                  *,seed:int=18) -> float:
    """Basis covariance of an arbitrary two-scalar-channel amplitude.

    This does NOT identify physical Yukawa/HH quartic tensors. It guarantees
    that any externally certified gauge-basis tensor can be contracted with
    K and rotated using U without introducing/omitting mass insertions.
    """
    _,U=positive_mass_matrix(A)
    K=spectral_kernel(A,mf_sq)
    rng=np.random.default_rng(seed)
    g=rng.normal(size=2)+1j*rng.normal(size=2)
    h=rng.normal(size=2)+1j*rng.normal(size=2)
    Q=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2))
    Q=0.5*(Q+Q.T)
    amp=np.vdot(np.kron(g,h),K@Q.reshape(4))
    pair=np.kron(U,U)
    gm=U.conj().T@g; hm=U.conj().T@h
    qm=U.conj().T@Q@U.conj()
    vals=np.array([signed_exact_master(mf_sq,positive_mass_matrix(A)[0][i],
                                       positive_mass_matrix(A)[0][j])
                   for i in range(2) for j in range(2)])
    other=np.vdot(np.kron(gm,hm),vals*qm.ravel())
    return float(abs(amp-other)/max(abs(amp),abs(other),1e-30))
