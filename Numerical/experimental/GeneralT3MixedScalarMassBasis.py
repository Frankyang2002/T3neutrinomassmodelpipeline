"""Experimental exact linear scalar mass-basis transport for ordinary T3 B_-1/C_-1.

Declared SU(2) convention: epsilon = [[0,1],[-1,0]],
  Phi=(S1, epsilon @ S2.conjugate()), Phi[a] = U[a,b] Chi[b],
  A=[[mS1Sq,mu2],[mu2*,mS2Sq]], U^dag A U=diag(mLowSq,mHighSq).

This transports *arbitrary supplied component tensors* and their original
index contractions (including a second independent contraction of identical
field content). It does not invent Matchete-CG or RGBeta component tensors.
No matching or field-registry mutation occurs.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from typing import Any, Mapping, Sequence

import numpy as np

EPS = np.array([[0., 1.], [-1., 0.]], dtype=np.complex128)
I2 = np.eye(2, dtype=np.complex128)
MIXED_MODELS = ('B_-1', 'C_-1')
SPECIES = ('H', 'Hdag', 'S1', 'S1dag', 'S2', 'S2dag')
OLD_HYPERCHARGE2 = {'H':1, 'Hdag':-1, 'S1':-1, 'S1dag':1,
                     'S2':1, 'S2dag':-1}
NEW_HYPERCHARGE2 = {'H':1, 'Hdag':-1, 'X':-1, 'Xdag':1,
                     'Y':-1, 'Ydag':1}


def _finite_real(value: float, name: str) -> float:
    v = float(value)
    if not np.isfinite(v):
        raise ValueError(f'{name} must be finite')
    return v


def _finite_complex(value: complex, name: str) -> complex:
    v = complex(value)
    if not (np.isfinite(v.real) and np.isfinite(v.imag)):
        raise ValueError(f'{name} must be finite')
    return v


def cjson(z: complex) -> dict[str, float]:
    z = _finite_complex(z, 'serialization scalar')
    return {'real': float(z.real), 'imag': float(z.imag)}


def matrix_json(m: np.ndarray) -> dict[str, list[list[float]]]:
    a = np.asarray(m, dtype=np.complex128)
    if a.ndim != 2 or not (np.isfinite(a.real).all() and np.isfinite(a.imag).all()):
        raise ValueError('Nonfinite/nonmatrix complex array')
    return {'real': a.real.tolist(), 'imag': a.imag.tolist()}


@dataclass(frozen=True)
class MixedDoubletMassBasis:
    model: str
    m_s1_sq: float
    m_s2_sq: float
    mu_sq: complex
    eigenvalues_sq: np.ndarray
    U: np.ndarray

    @classmethod
    def build(cls, model: str, m_s1_sq: float, m_s2_sq: float,
              mu_sq: complex) -> 'MixedDoubletMassBasis':
        if model not in MIXED_MODELS:
            raise ValueError('Only T3 B_-1/C_-1 have this mixed-doublet basis')
        m1 = _finite_real(m_s1_sq, 'mS1Sq')
        m2 = _finite_real(m_s2_sq, 'mS2Sq')
        mu = _finite_complex(mu_sq, 'mu_sq')
        A = np.array([[m1, mu], [mu.conjugate(), m2]], dtype=np.complex128)
        values, U = np.linalg.eigh(A)
        # A deterministic phase convention for reproducible logs/output.
        for a in range(2):
            pivot = int(np.argmax(abs(U[:, a])))
            if abs(U[pivot, a]) > 0:
                U[:, a] *= np.conjugate(U[pivot, a] / abs(U[pivot, a]))
        obj = cls(model, m1, m2, mu, values, U)
        obj.verify()
        return obj

    @property
    def A(self) -> np.ndarray:
        return np.array([[self.m_s1_sq, self.mu_sq],
                         [self.mu_sq.conjugate(), self.m_s2_sq]], dtype=complex)

    @property
    def tachyonic(self) -> bool:
        return bool(np.any(self.eigenvalues_sq < 0.))

    @property
    def degenerate(self) -> bool:
        return bool(np.isclose(self.eigenvalues_sq[0], self.eigenvalues_sq[1],
                               rtol=1e-12, atol=1e-12))

    @property
    def mixing_sin_2theta(self) -> float:
        delta = self.m_s1_sq - self.m_s2_sq
        length = float(np.hypot(delta, 2*abs(self.mu_sq)))
        return 0.0 if length == 0.0 else 2*abs(self.mu_sq)/length

    def verify(self, atol: float = 1e-8) -> None:
        scale = max(abs(self.m_s1_sq), abs(self.m_s2_sq),
                    abs(self.mu_sq), 1.0)
        if not np.allclose(self.U.conj().T @ self.U, I2,
                           rtol=1e-12, atol=1e-12):
            raise RuntimeError('Mass-basis transformation not unitary')
        if not np.allclose(self.U.conj().T @ self.A @ self.U,
                           np.diag(self.eigenvalues_sq), rtol=1e-11,
                           atol=max(atol, scale*1e-12)):
            raise RuntimeError('Mixed-doublet mass matrix not diagonalized')
        avg = (self.m_s1_sq+self.m_s2_sq)/2.0
        gap = np.hypot((self.m_s1_sq-self.m_s2_sq)/2.0, abs(self.mu_sq))
        if not np.allclose(self.eigenvalues_sq, [avg-gap,avg+gap],
                           rtol=1e-11, atol=max(atol,scale*1e-12)):
            raise RuntimeError('Analytic scalar mass eigenvalue mismatch')

    def original_fields(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Return (S1,S2) from two hypercharge -1/2 eigen-doublets."""
        x = np.asarray(x, dtype=complex)
        y = np.asarray(y, dtype=complex)
        if x.shape != (2,) or y.shape != (2,):
            raise ValueError('X and Y must have two SU(2) components')
        s1 = self.U[0,0]*x + self.U[0,1]*y
        tilde_s2 = self.U[1,0]*x + self.U[1,1]*y
        s2 = -EPS @ tilde_s2.conjugate()
        return s1, s2

    def eigen_fields(self, s1: np.ndarray, s2: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        s1 = np.asarray(s1, dtype=complex)
        s2 = np.asarray(s2, dtype=complex)
        if s1.shape != (2,) or s2.shape != (2,):
            raise ValueError('S1 and S2 must have two SU(2) components')
        vec = self.U.conj().T @ np.stack([s1, EPS @ s2.conjugate()])
        return vec[0], vec[1]

    def old_potential(self, s1: np.ndarray, s2: np.ndarray) -> float:
        value = (self.m_s1_sq*np.vdot(s1,s1) + self.m_s2_sq*np.vdot(s2,s2)
                 + 2*np.real(self.mu_sq*np.vdot(s1,EPS @ s2.conjugate())))
        return float(np.real(value))

    def new_potential(self, x: np.ndarray, y: np.ndarray) -> float:
        return float(self.eigenvalues_sq[0]*np.vdot(x,x).real
                     + self.eigenvalues_sq[1]*np.vdot(y,y).real)

    def info(self) -> dict:
        return {'model':self.model,'basis':'(S1, epsilon*conjugate(S2)) = U @ (X,Y)',
                'epsilon':EPS.real.tolist(),
                'eigenstate_hypercharge':-0.5,
                'mass_squared_flavour_matrix_gev2':matrix_json(self.A),
                'unitary_flavour_rotation_U':matrix_json(self.U),
                'masses_squared_gev2':self.eigenvalues_sq.tolist(),
                'masses_gev':[float(np.sqrt(z)) if z>=0 else None
                              for z in self.eigenvalues_sq],
                'tachyonic':self.tachyonic,
                'degenerate_basis_nonunique':self.degenerate,
                'sin_2theta':self.mixing_sin_2theta,
                'unitarity_verified':True,'diagonal_mass_verified':True,
                'matchete_field_redefinition_installed':False,
                'rgbeta_matchete_component_map_certified':False,
                'threshold_matching_performed':False}


def field_map(field: str, U: np.ndarray) -> tuple[tuple[str, np.ndarray], ...]:
    """OLD component field f_i = sum_out M[i,j] NEW component field_j.

    For S2 and S2dag the -epsilon factor is indispensable. Omitting it is
    wrong even if the flavour-space mass matrix diagonalizes numerically.
    """
    if field == 'H': return (('H',I2.copy()),)
    if field == 'Hdag': return (('Hdag',I2.copy()),)
    if field == 'S1':
        return tuple((n, U[0,k]*I2) for k,n in enumerate(('X','Y')))
    if field == 'S1dag':
        return tuple((n, U[0,k].conjugate()*I2)
                     for k,n in enumerate(('Xdag','Ydag')))
    if field == 'S2':
        return tuple((n, -U[1,k].conjugate()*EPS)
                     for k,n in enumerate(('Xdag','Ydag')))
    if field == 'S2dag':
        return tuple((n, -U[1,k]*EPS)
                     for k,n in enumerate(('X','Y')))
    raise ValueError('Unknown species: '+str(field))


def transport_operator(fields: Sequence[str], U: np.ndarray) -> list[dict[str, Any]]:
    """Ordered tensor pushforward for a symbolic interaction T_{i...} fields_i.

    The returned term's index maps do NOT contain the original invariant tensor.
    The same T_{i...} must be contracted against the term maps. This retains
    each distinct source CG tensor, including the B/C 2-dimensional sector.
    """
    old_fields = tuple(fields)
    if not old_fields or len(old_fields)>4 or any(f not in SPECIES for f in old_fields):
        raise ValueError('Need 1-4 supported ordered scalar field factors')
    if not isinstance(U,np.ndarray) or U.shape!=(2,2):
        raise ValueError('U must be 2x2 numpy matrix')
    terms = []
    for choice in product(*(field_map(name,U) for name in old_fields)):
        names = [item[0] for item in choice]
        if sum(OLD_HYPERCHARGE2[f] for f in old_fields) != sum(NEW_HYPERCHARGE2[f] for f in names):
            raise ValueError('Field replacement violates U(1) hypercharge')
        terms.append({'fields':names,
                      'index_maps':[matrix_json(item[1]) for item in choice]})
    return terms


def _tensor_contract(tensor: np.ndarray, vectors: Sequence[np.ndarray]) -> complex:
    out = np.asarray(tensor,dtype=complex)
    for v in vectors:
        out = np.tensordot(out,v,axes=(0,0))
    return complex(out)


def evaluate_transported_tensor(tensor: np.ndarray, fields: Sequence[str],
                                U: np.ndarray, x: np.ndarray,
                                y: np.ndarray, h: np.ndarray) -> tuple[complex,complex]:
    """Compare arbitrary original SU(2) tensor vs its full rotated expansion."""
    fields = tuple(fields)
    n=len(fields)
    t=np.asarray(tensor,dtype=complex)
    if t.shape!=(2,)*n:
        raise ValueError('One doublet component index required per field')
    # For this check, the user-provided U defines the basis; don't use basis.U.
    s1=(U[0,0]*x+U[0,1]*y)
    s2=-EPS @ (U[1,0]*x+U[1,1]*y).conjugate()
    old={'H':h,'Hdag':h.conjugate(),'S1':s1,'S1dag':s1.conjugate(),
         'S2':s2,'S2dag':s2.conjugate()}
    new={'H':h,'Hdag':h.conjugate(),'X':x,'Xdag':x.conjugate(),
         'Y':y,'Ydag':y.conjugate()}
    before=_tensor_contract(t,[old[f] for f in fields])
    after=0j
    for choice in product(*(field_map(f,U) for f in fields)):
        new_vectors=[M @ new[name] for name,M in choice]
        after += _tensor_contract(t,new_vectors)
    return before,after


def verify_tensor_transport(fields: Sequence[str], U: np.ndarray,
                            *, seed: int = 4, trials: int = 5) -> float:
    """An identity test for arbitrary complex coefficient tensors, not gauge cert."""
    rng=np.random.default_rng(seed)
    residual=0.0
    for _ in range(trials):
        z=lambda *size: rng.normal(size=size)+1j*rng.normal(size=size)
        old,new=evaluate_transported_tensor(z(*(2,)*len(fields)),fields,U,
                                             z(2),z(2),z(2))
        err=abs(old-new)/max(1.,abs(old),abs(new))
        residual=max(residual,float(err))
    if residual>3e-12:
        raise RuntimeError(f'Operator tensor rotation failed: {fields} residual={residual}')
    return residual


def quadratic_mass_tensor_transport(basis: MixedDoubletMassBasis) -> dict[str, Any]:
    """Transport EVERY scalar mass term, including holomorphic mixing + h.c.

    This is an independent index-level check of the epsilon convention. The
    resulting coefficients of Chi_a^dag Chi_b must equal the diagonal masses.
    """
    sources=(
        (basis.m_s1_sq,('S1dag','S1'),I2),
        (basis.m_s2_sq,('S2dag','S2'),I2),
        (basis.mu_sq,('S1dag','S2dag'),EPS),
        (basis.mu_sq.conjugate(),('S1','S2'),EPS),
    )
    entries=np.zeros((2,2),dtype=np.complex128)
    names={'X':0,'Y':1,'Xdag':0,'Ydag':1}
    for coefficient,fields,tensor in sources:
        for name0,m0 in field_map(fields[0],basis.U):
            for name1,m1 in field_map(fields[1],basis.U):
                component_tensor=coefficient*(m0.T@tensor@m1)
                coeff=complex(np.trace(component_tensor)/2.)
                if not np.allclose(component_tensor,coeff*I2,
                                   rtol=1e-12,atol=1e-7):
                    raise RuntimeError('Quadratic transformed tensor not an SU2 singlet')
                if name0.endswith('dag') and not name1.endswith('dag'):
                    entries[names[name0],names[name1]]+=coeff
                elif name1.endswith('dag') and not name0.endswith('dag'):
                    # Bosonic fields commute; swapping their order also swaps
                    # the component indices, which are proportional to I here.
                    entries[names[name1],names[name0]]+=coeff
                else:
                    raise RuntimeError('Rotated mixed mass term violates charge')
    residual=float(np.max(np.abs(entries-np.diag(basis.eigenvalues_sq))))
    scale=max(1.,abs(basis.m_s1_sq),abs(basis.m_s2_sq),abs(basis.mu_sq))
    if residual>2e-10*scale:
        raise RuntimeError('Full bilinear tensor transport did not diagonalize mass')
    return {'rotated_bilinear_coefficient_matrix':matrix_json(entries),
            'absolute_offdiagonal_max_gev2':float(np.max(np.abs(entries-np.diag(np.diag(entries))))),
            'mass_tensor_diagonalization_residual_gev2':residual,
            'holomorphic_term_and_hermitian_conjugate_included':True}


def verify_su2_covariance(basis: MixedDoubletMassBasis, trials:int=16,
                          *,seed:int=0) -> float:
    """Check epsilon*g.conj() == g*epsilon and inverse/forward gauge action."""
    rng=np.random.default_rng(seed)
    max_error=0.
    for _ in range(trials):
        q=rng.normal(size=4);q/=np.linalg.norm(q)
        a,b,c,d=q
        g=np.array([[a+1j*b,c+1j*d],[-c+1j*d,a-1j*b]])
        x=rng.normal(size=2)+1j*rng.normal(size=2)
        y=rng.normal(size=2)+1j*rng.normal(size=2)
        s1,s2=basis.original_fields(x,y)
        t1,t2=basis.original_fields(g@x,g@y)
        err=max(np.max(abs(t1-g@s1)),np.max(abs(t2-g@s2)),
                np.max(abs(EPS@g.conjugate()-g@EPS)))
        max_error=max(max_error,float(err))
    if max_error>3e-12:
        raise RuntimeError('SU2-doublet covariance check failed')
    return max_error
