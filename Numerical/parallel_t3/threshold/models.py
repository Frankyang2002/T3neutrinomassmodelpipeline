"""Read-only-source-derived T3 model registry and scalar mass-basis constructors.

Conventions: real fields of each complex multiplet are interleaved
(Re S[0], Im S[0], ...), where S=(R+i I)/sqrt(2). All mass matrices
are real symmetric Hessians of V2 and have units GeV^2.  Source-native
self-bilinear and mixed-doublet field maps are EXPLICIT, not assumed to be
identical to the A/B tensor coordinate labels.

Reference source inspected: T3ModelCatalog.wl, special_mass_v6.py,
GeneralT3MixedScalarMassBasis.py (public main a944e003...).
"""
from __future__ import annotations
from dataclasses import dataclass
import re
import numpy as np

DIMS = {"A":(1,3,2),"B":(2,2,1),"C":(2,2,3),"D":(3,1,2),"E":(3,3,2)}
ALPHAS = {"A":(-4,-2,0),"B":(-3,-1,1),"C":(-3,-1,1),
          "D":(-2,0,2),"E":(-4,-2,0,2)}
ALL_MODELS = tuple(f"{k}_{a:+d}" if a >= 0 else f"{k}_{a}"
                   for k, aa in ALPHAS.items() for a in aa)
MIXED = frozenset(("B_-1","C_-1"))
SELF_S1 = frozenset(("A_+0","D_+0","E_+0"))
SELF_S2 = frozenset(("A_-2","D_-2","E_-2"))
EPS = np.array([[0., 1.],[-1.,0.]], dtype=complex)


def normalize_model(value: str) -> str:
    value = str(value).strip().replace('T3-', '')
    match = re.fullmatch(r'([A-E])_([mp]\d+|[+-]?\d+)', value)
    if not match:
        raise ValueError(f"Unsupported model label: {value}")
    part = match.group(2)
    alpha = (-int(part[1:]) if part.startswith('m') else
             int(part[1:]) if part.startswith('p') else int(part))
    model = f'{match.group(1)}_{alpha:+d}' if alpha >= 0 else f'{match.group(1)}_{alpha}'
    if model not in ALL_MODELS:
        raise ValueError(f'Not one of the production T3 models: {model}')
    return model


def dimensions(model: str) -> tuple[int,int,int]:
    return DIMS[normalize_model(model)[0]]


def n_real(model: str) -> int:
    d1,d2,_ = dimensions(model)
    return 2*(d1+d2)


def special_kind(model: str) -> str:
    model=normalize_model(model)
    if model in MIXED: return 'mixed_doublets'
    if model in SELF_S1: return 'self_s1'
    if model in SELF_S2: return 'self_s2'
    return 'ordinary'


def scalar_labels(model: str, basis: str='canonical_interleaved') -> list[str]:
    model=normalize_model(model)
    d1,d2,_=dimensions(model)
    if basis == 'canonical_interleaved':
        return [f'{name}[{i}].{part}' for name,d in [('S1',d1),('S2',d2)]
                for i in range(1,d+1) for part in ('R','I')]
    if basis != 'native_named':
        raise ValueError('basis must be canonical_interleaved or native_named')
    if model in MIXED:
        return [f'{s}[{i}].{p}' for s in ('X','Y') for i in (1,2) for p in ('R','I')]
    if model in SELF_S1:
        return ([f'S1Minus[{i}]' for i in range(1,d1+1)] +
                [f'S1Plus[{i}]' for i in range(1,d1+1)] +
                [f'S2[{i}].{p}' for i in range(1,d2+1) for p in ('R','I')])
    if model in SELF_S2:
        return ([f'S1[{i}].{p}' for i in range(1,d1+1) for p in ('R','I')] +
                [f'S2Minus[{i}]' for i in range(1,d2+1)] +
                [f'S2Plus[{i}]' for i in range(1,d2+1)])
    return scalar_labels(model,'canonical_interleaved')


def _finite(value: object, name: str) -> float:
    a=float(value)
    if not np.isfinite(a): raise ValueError(f'Nonfinite {name}')
    return a


def _interleave_from_named(dim: int) -> np.ndarray:
    """S[k]=(Minus[k]+i Plus[k])/sqrt2, no spherical W."""
    r=np.zeros((2*dim,2*dim))
    for j in range(dim):
        r[2*j,j]=1.
        r[2*j+1,dim+j]=1.
    return r


def _real_pack(z: np.ndarray) -> np.ndarray:
    x=np.asarray(z,complex)
    return np.stack((x.real,x.imag),axis=-1).reshape(-1)*np.sqrt(2.)


def _real_unpack(r: np.ndarray) -> np.ndarray:
    return (r[::2]+1j*r[1::2])/np.sqrt(2.)


def _mixed_rotation(unitary: np.ndarray) -> np.ndarray:
    """(S1,S2) canonical real = R @ (X,Y) canonical real."""
    out=np.zeros((8,8))
    for j in range(8):
        e=np.zeros(8);e[j]=1.
        x,y=_real_unpack(e[:4]),_real_unpack(e[4:])
        s1=unitary[0,0]*x+unitary[0,1]*y
        t2=unitary[1,0]*x+unitary[1,1]*y
        s2=-EPS @ t2.conjugate()
        out[:,j]=np.r_[_real_pack(s1),_real_pack(s2)]
    return out


@dataclass(frozen=True)
class ScalarMassAdapter:
    model: str
    kind: str
    mass2_canonical_gev2: np.ndarray
    R_canonical_from_native: np.ndarray
    mass2_native_gev2: np.ndarray
    native_labels: tuple[str,...]
    extra: dict


def from_parameters(model: str, m_s1_sq: float, m_s2_sq: float,
                    mu_sq: complex | None=None) -> ScalarMassAdapter:
    """Generate all 16 mass structures; explicit mu is REQUIRED for special cases.

    This is an unbroken-gauge quadratic-parameter construction, not an
    EWSB-resummed physical-mass construction. It does not generate quartics.
    """
    model=normalize_model(model)
    d1,d2,_=dimensions(model)
    m1,m2=_finite(m_s1_sq,'mS1Sq'),_finite(m_s2_sq,'mS2Sq')
    kind=special_kind(model)
    n=2*(d1+d2)
    R=np.eye(n)
    mass=np.diag(np.r_[np.full(2*d1,m1),np.full(2*d2,m2)])
    extra={'source':'T3 native ordinary interleaved real-component mass convention',
           'native_phase_alignment_certified':False}
    if kind != 'ordinary':
        if mu_sq is None:
            raise ValueError(f'{model}: complex quadratic parameter mu_sq must be explicitly supplied')
        mu=complex(mu_sq)
        if not np.isfinite([mu.real,mu.imag]).all():
            raise ValueError('Nonfinite complex mu_sq')
        extra.update(mu_sq_real=mu.real,mu_sq_imag=mu.imag)
        if kind == 'mixed_doublets':
            # (S1, epsilon S2*) = U (X,Y), with mu in A[0,1].
            A=np.array([[m1,mu],[mu.conjugate(),m2]],complex)
            ev,U=np.linalg.eigh(A)
            for j in range(2):
                k=int(np.argmax(np.abs(U[:,j])))
                U[:,j]*=np.conjugate(U[k,j]/abs(U[k,j]))
            R=_mixed_rotation(U)
            mass=R @ np.diag(np.repeat(ev,4)) @ R.T
            extra.update(complex_flavour_rotation_real=U.real.tolist(),
                         complex_flavour_rotation_imag=U.imag.tolist(),
                         eigen_doublet_masses2_gev2=ev.tolist(),
                         field_map='(S1, epsilon conjugate(S2)) = U (X,Y); epsilon_12=+1')
        else:
            dim=d1 if kind=='self_s1' else d2
            offset=0 if kind=='self_s1' else 2*d1
            m=m1 if kind=='self_s1' else m2
            # The native Matchete source uses +mu (Sdag)^2 in L_int,
            # so the real Hessian of V has -2Re(mu), -2Im(mu).
            local=np.array([[m-2*mu.real,-2*mu.imag],
                            [-2*mu.imag,m+2*mu.real]])
            for j in range(dim):
                ids=[offset+2*j,offset+2*j+1]
                mass[np.ix_(ids,ids)]=local
            R[offset:offset+2*dim,offset:offset+2*dim]=_interleave_from_named(dim)
            extra.update(split_masses2_gev2=[m-2*abs(mu),m+2*abs(mu)],
                         phase_handling='Complex mu retained as R/I offdiagonal; not set real',
                         field_map='S[k]=(Minus[k]+i Plus[k])/sqrt2, unchanged Cartesian SU2 index')
    if not np.allclose(R.T@R,np.eye(n),atol=2e-12):
        raise RuntimeError('Nonorthogonal scalar native map')
    mass=.5*(mass+mass.T)
    eig=np.linalg.eigvalsh(mass)
    if eig[0]<=0:raise ValueError(f'{model}: non-positive mass-squared eigenstate {eig[0]:g} GeV^2')
    native=R.T@mass@R
    if not np.allclose(np.linalg.eigvalsh(native),eig,rtol=1e-11,atol=1e-7):
        raise RuntimeError('Scalar field map does not preserve mass eigenvalues')
    return ScalarMassAdapter(model,kind,mass,R,native,
                             tuple(scalar_labels(model,'native_named')),extra)
