"""Exact REAL-scalar tensor transport for the conditional B_-1/C_-1 Y threshold.

Original real basis: (H_R1,H_I1,H_R2,H_I2,
                      S1_R1,S1_I1,S1_R2,S1_I2,
                      S2_R1,S2_I1,S2_R2,S2_I2).
New real basis:      (H_R1,H_I1,H_R2,H_I2,
                      X_R1,X_I1,X_R2,X_I2,
                      Y_R1,Y_I1,Y_R2,Y_I2).

For complex doublets z_a=(r_a+i s_a)/sqrt(2), the original fields obey
    S1 = U[0,0]*X + U[0,1]*Y
    S2 = -epsilon @ conj(U[1,0]*X + U[1,1]*Y).
The real matrix T maps old real fields to new: real_old=T @ real_new.
A Wilson tensor defined as C_ijabf psi_i psi_j real_old[a] real_old[b]
therefore transforms as C'_ijABf=T[a,A] T[b,B] C_ijabf.
No assumption is made about the fermion tensor/its charge conjugations.

Tree projection onto EFT2 (H,X) is just C'[:,:,:8,:8,:] AFTER rotation.
This is a field substitution, NOT finite one-loop Wilson matching.
"""
from __future__ import annotations

import numpy as np

EPS = np.array([[0.0, 1.0], [-1.0, 0.0]])
PAIRS = ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))


def require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError('T3 B/C Y projection: ' + message)


def _real_block(a: complex, conjugate: bool = False) -> np.ndarray:
    """Map (Re z, Im z) -> (Re a*z^(*) , Im a*z^(*)), without extra sqrt(2)."""
    a = complex(a)
    if conjugate:
        return np.array([[a.real, a.imag], [a.imag, -a.real]], dtype=float)
    return np.array([[a.real, -a.imag], [a.imag, a.real]], dtype=float)


def real_doublet(g: np.ndarray) -> np.ndarray:
    """Realification of a general 2x2 COMPLEX-linear action on a doublet."""
    g = np.asarray(g, dtype=complex)
    require(g.shape == (2, 2), 'doublet transformation must be 2x2')
    out = np.zeros((4, 4), dtype=float)
    for a in range(2):
        for b in range(2):
            out[2*a:2*a+2, 2*b:2*b+2] = _real_block(g[a,b])
    return out


def mixing_real_matrix(U: np.ndarray) -> np.ndarray:
    """12x12 orthogonal old <- new real-field transformation (H,S1,S2) <- (H,X,Y)."""
    U = np.asarray(U, dtype=complex)
    require(U.shape == (2, 2), 'U must be 2x2')
    require(np.isfinite(U).all(), 'nonfinite U')
    require(np.linalg.norm(U.conj().T @ U - np.eye(2)) < 5e-11,
            'U is not unitary')
    T = np.zeros((12, 12), dtype=float)
    T[:4,:4] = np.eye(4)
    for k, new_start in enumerate((4, 8)):
        for component in range(2):
            # S1_component = U[0,k] * eigen_component
            T[4+2*component:6+2*component,
              new_start+2*component:new_start+2*component+2] = _real_block(U[0,k])
            # S2_component = -epsilon_component,b * conj(U[1,k]*eigen_b)
            for b in range(2):
                coefficient = -EPS[component,b] * U[1,k].conjugate()
                T[8+2*component:10+2*component,
                  new_start+2*b:new_start+2*b+2] += _real_block(coefficient,conjugate=True)
    require(np.linalg.norm(T.T @ T - np.eye(12),ord=np.inf) < 1e-10,
            'real mixing map is not orthogonal')
    require(np.linalg.norm(T @ T.T - np.eye(12),ord=np.inf) < 1e-10,
            'real mixing map is not complete')
    return T


def real_field_vector(h: np.ndarray, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """z=(Re+i Im)/sqrt(2) gives interleaved real components in H,X,Y order."""
    vec = []
    for name, field in (('H',h),('X',x),('Y',y)):
        field = np.asarray(field,dtype=complex)
        require(field.shape == (2,), name + ' must be a complex SU2 doublet')
        for z in field:
            vec.extend([float(np.sqrt(2)*z.real), float(np.sqrt(2)*z.imag)])
    return np.asarray(vec,dtype=float)


def check_field_identity(U: np.ndarray, T: np.ndarray, *, seed: int = 7,
                         trials: int = 16) -> float:
    """Independent complex-field check; detects epsilon/conjugation/order errors."""
    rng = np.random.default_rng(seed)
    maxerr = 0.0
    for _ in range(trials):
        h,x,y = [rng.normal(size=2) + 1j*rng.normal(size=2) for _ in range(3)]
        s1 = U[0,0]*x + U[0,1]*y
        s2 = -EPS @ np.conj(U[1,0]*x + U[1,1]*y)
        old_expected = real_field_vector(h,s1,s2)
        old_observed = T @ real_field_vector(h,x,y)
        err = np.linalg.norm(old_expected-old_observed) / max(1.,np.linalg.norm(old_expected))
        maxerr = max(maxerr,float(err))
    require(maxerr < 3e-12,'complex epsilon/conjugation field identity failed')
    return maxerr


def _block_diag(blocks: tuple[np.ndarray,...]) -> np.ndarray:
    n=sum(b.shape[0] for b in blocks)
    out=np.zeros((n,n),dtype=float)
    offset=0
    for b in blocks:
        m=b.shape[0];out[offset:offset+m,offset:offset+m]=b;offset+=m
    return out


def check_gauge_covariance(T: np.ndarray, *, seed: int = 13,
                           trials: int = 16) -> float:
    """Check SU(2) and U(1) actions intertwine, and H/X projector is equivariant.

    H,S2 have hypercharge +1/2; S1,X,Y have -1/2 (Q=T3+Y).
    Fermion gauge invariance is NOT independently certified by this scalar check.
    """
    rng=np.random.default_rng(seed)
    P=np.eye(12,dtype=float)[:8,:]
    maxerr=0.0
    for _ in range(trials):
        q=rng.normal(size=4);q/=np.linalg.norm(q)
        a,b,c,d=q
        g=np.array([[a+1j*b,c+1j*d],[-c+1j*d,a-1j*b]])
        theta=float(rng.uniform(-np.pi,np.pi))
        plus=real_doublet(g*np.exp(0.5j*theta))
        minus=real_doublet(g*np.exp(-0.5j*theta))
        old=_block_diag((plus,minus,plus))
        new=_block_diag((plus,minus,minus))
        err=np.linalg.norm(old @ T - T @ new,ord=np.inf)
        err2=np.linalg.norm(P @ new - new[:8,:8] @ P,ord=np.inf)
        maxerr=max(maxerr,float(err),float(err2))
    require(maxerr<3e-12,'SU(2)xU(1) scalar covariance/projection failed')
    return maxerr


def symmetry_residual(C: np.ndarray) -> dict[str,float]:
    c=np.asarray(C,dtype=complex)
    require(c.ndim==5 and c.shape[:2]==(3,3) and c.shape[-1]==6 and
            c.shape[2]==c.shape[3] and c.shape[2] in (8,12),
            'Wilson shape must be (3,3,12|8,12|8,6)')
    require(np.isfinite(c).all(),'nonfinite Wilson tensor')
    norm=max(float(np.linalg.norm(c)),1e-300)
    return {'fermion_pair':float(np.linalg.norm(c-c.swapaxes(0,1))/norm),
            'scalar_pair':float(np.linalg.norm(c-c.swapaxes(2,3))/norm)}


def rotate_wilson(C: np.ndarray, T: np.ndarray) -> np.ndarray:
    """Quadratic field substitution on BOTH real-scalar indices (no complex conjugates)."""
    c=np.asarray(C,dtype=complex)
    require(c.shape==(3,3,12,12,6),'EFT1 Wilson shape changed')
    require(T.shape==(12,12) and np.isfinite(T).all(),'invalid real map')
    before=symmetry_residual(c)
    require(max(before.values())<1e-9,'input Wilson pair symmetries invalid')
    out=np.einsum('ap,ijabf,bq->ijpqf',T,c,T,optimize=True)
    after=symmetry_residual(out)
    require(max(after.values())<1e-9,'rotated Wilson pair symmetries invalid')
    return out


def light_projection(C_mass: np.ndarray) -> np.ndarray:
    require(np.asarray(C_mass).shape==(3,3,12,12,6),'expected full mass-basis Wilson tensor')
    return np.asarray(C_mass[:,:,:8,:8,:],dtype=complex).copy()


def check_tensor_identity(C: np.ndarray, C_mass: np.ndarray, T: np.ndarray,
                          *, seed:int=29, trials:int=10) -> float:
    """Numerical identity of original and rotated scalar polynomials for all flavor channels."""
    rng=np.random.default_rng(seed)
    maxerr=0.0
    for _ in range(trials):
        new=rng.normal(size=12)
        old=T @ new
        for i,j,f in ((0,0,0),(0,1,1),(1,2,4),(2,2,5)):
            lhs=old @ C[i,j,:,:,f] @ old
            rhs=new @ C_mass[i,j,:,:,f] @ new
            err=abs(lhs-rhs)/max(abs(lhs),abs(rhs),1e-25)
            maxerr=max(maxerr,float(err))
    require(maxerr < 5e-11,'rotated tensor polynomial identity failed')
    return maxerr


def check_light_identity(C_mass: np.ndarray, light: np.ndarray, *,seed:int=31,
                         trials:int=10) -> float:
    rng=np.random.default_rng(seed)
    maxerr=0.0
    for _ in range(trials):
        new8=rng.normal(size=8)
        new12=np.r_[new8,np.zeros(4)]
        lhs=np.einsum('ijabf,a,b->ijf',C_mass,new12,new12,optimize=True)
        rhs=np.einsum('ijabf,a,b->ijf',light,new8,new8,optimize=True)
        err=np.linalg.norm(lhs-rhs)/max(np.linalg.norm(lhs),np.linalg.norm(rhs),1e-25)
        maxerr=max(maxerr,float(err))
    require(maxerr<5e-12,'EFT2 light-field restriction failed')
    return maxerr


def classify_survivors(light: np.ndarray) -> dict[str,dict]:
    """Norms and occupied entries, without asserting which terms are gauge invariant."""
    require(light.shape==(3,3,8,8,6),'invalid EFT2 tensor')
    blocks={'HH':(slice(0,4),slice(0,4)),
            'HX':(slice(0,4),slice(4,8)),
            'XX':(slice(4,8),slice(4,8))}
    result={}
    total=max(float(np.max(abs(light))),1e-300)
    for name,(a,b) in blocks.items():
        data=light[:,:,a,b,:]
        result[name]={'frobenius_norm':float(np.linalg.norm(data)),
                      'max_abs':float(np.max(abs(data))),
                      'count_above_relative_1e_minus_12':int(np.count_nonzero(abs(data)>1e-12*total))}
    return result
