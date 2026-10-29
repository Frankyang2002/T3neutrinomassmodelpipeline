"""Exact spherical-to-real mass basis for the six neutral T3 self-bilinears.

For d=1 or 3 and the declared spherical invariant C, the potential is
  V2=m^2 S^dag S + mu2 (S^dag)^T C (S^dag) + h.c. (NO extra 1/2).

Let C=W W^T, W^T C W=I; S=exp(i arg(mu2)/2) W (R+i I)/sqrt(2),
where R and I are independent real d-component scalar fields.  Their
mass squared values are m^2+2|mu2| and m^2-2|mu2| respectively.

The NORMAL S-S^dag propagator is the average of the two eigenpropagators.
The ANOMALOUS S-S propagator is nonzero and is NOT automatically included
by using the ordinary topology's C5 with averaged scalar propagators.

This is an explicit declared spherical convention, NOT a certified native
Matchete/Cartesian Clebsch-Gordan coefficient map.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from typing import Any, Mapping, Sequence
import numpy as np


def invariant_metric(dim: int) -> np.ndarray:
    if dim == 1:
        return np.eye(1)
    if dim == 3:
        return np.array([[0., 0., 1.], [0., -1., 0.], [1., 0., 0.]])
    raise ValueError('Only integer-isospin scalar dimensions 1 and 3 are supported')


def real_intertwiner(dim: int) -> np.ndarray:
    """Deterministic unitary W: W^T C W = I, C W* = W."""
    C=invariant_metric(dim)
    eigen,V=np.linalg.eigh(C)
    phases=np.where(eigen > 0,1.0+0j,1.0j)
    W=V.astype(complex)@np.diag(phases)
    I=np.eye(dim)
    if not (np.allclose(W.conj().T@W,I,atol=1e-13) and
            np.allclose(W.T@C@W,I,atol=1e-13) and
            np.allclose(C@W.conj(),W,atol=1e-13)):
        raise RuntimeError('SU(2) spherical-to-real unitary intertwiner failed')
    return W


def _finite_complex(z: complex, name: str) -> complex:
    out=complex(z)
    if not np.isfinite([out.real,out.imag]).all():
        raise ValueError(f'{name} must be finite')
    return out


@dataclass(frozen=True)
class SelfBilinearMassBasis:
    species: str
    dim: int
    diagonal_mass_squared: float
    mu_squared: complex
    eigen_masses_squared: tuple[float,float]
    intertwiner: np.ndarray
    phase: complex

    @classmethod
    def build(cls, species: str, dim: int, mass_squared: float,
              mu_squared: complex) -> 'SelfBilinearMassBasis':
        if species not in ('S1','S2'):
            raise ValueError('Self-bilinear must belong to S1 or S2')
        m=float(mass_squared)
        mu=_finite_complex(mu_squared,'mu_squared')
        if not np.isfinite(m):
            raise ValueError('Diagonal scalar squared mass is nonfinite')
        W=real_intertwiner(dim)
        phase=np.exp(0.5j*np.angle(mu)) if mu != 0 else (1+0j)
        masses=(float(m+2*abs(mu)),float(m-2*abs(mu)))
        basis=cls(species,dim,m,mu,masses,W,complex(phase))
        basis.verify()
        return basis

    @property
    def tachyonic(self) -> bool:
        return min(self.eigen_masses_squared) <= 0

    def verify(self) -> None:
        m=self.diagonal_mass_squared
        mu=self.mu_squared
        C=invariant_metric(self.dim)
        B=2*mu.conjugate()*C
        H=np.block([[m*np.eye(self.dim)+B.real,-B.imag],
                    [-B.imag,m*np.eye(self.dim)-B.real]])
        expected=np.sort(np.repeat(self.eigen_masses_squared,self.dim))
        if not np.allclose(np.linalg.eigvalsh(H),expected,rtol=1e-12,atol=1e-8):
            raise RuntimeError('Self-bilinear Hessian and mass eigenvalues disagree')
        if not np.allclose(self.phase*self.phase*mu.conjugate(),abs(mu),rtol=1e-12,atol=1e-9):
            raise RuntimeError('Self-bilinear field phase does not diagonalize mu2')

    def fields(self, real: Sequence[float], imaginary: Sequence[float]) -> np.ndarray:
        R=np.asarray(real,dtype=float);I=np.asarray(imaginary,dtype=float)
        if R.shape!=(self.dim,) or I.shape!=(self.dim,):
            raise ValueError('Real multiplets have incorrect SU(2) dimensions')
        return self.phase*(self.intertwiner@(R+1j*I))/np.sqrt(2)

    def potential_old(self, scalar: np.ndarray) -> float:
        s=np.asarray(scalar,dtype=complex)
        if s.shape!=(self.dim,):raise ValueError('Incorrect scalar component count')
        # mu2 (Sdag)^T C Sdag + h.c. = 2 Re(mu2 S*^T C S*)
        return float(self.diagonal_mass_squared*np.vdot(s,s).real+
                     2*np.real(self.mu_squared*(s.conjugate().T@
                                               invariant_metric(self.dim)@s.conjugate())))

    def potential_mass_basis(self, real: Sequence[float], imaginary: Sequence[float]) -> float:
        R=np.asarray(real,dtype=float);I=np.asarray(imaginary,dtype=float)
        if R.shape!=(self.dim,) or I.shape!=(self.dim,):
            raise ValueError('Real multiplet dimensions differ')
        a,b=self.eigen_masses_squared
        return float(0.5*(a*(R@R)+b*(I@I)))

    def propagators(self, euclidean_p_squared: float) -> tuple[np.ndarray,np.ndarray]:
        """Return <S_i S_j*> and <S_i S_j> at Euclidean momentum p.

        Normal G_N=I/2 (1/(p²+m_+²)+1/(p²+m_-²));
        anomalous G_A=e^{i arg(mu)} C/2 (1/(p²+m_+²)-1/(p²+m_-²)).
        """
        p=float(euclidean_p_squared)
        if not np.isfinite(p) or p<0:
            raise ValueError('Euclidean momentum squared must be nonnegative and finite')
        a,b=self.eigen_masses_squared
        if p+a<=0 or p+b<=0:
            raise ValueError('Non-positive scalar Euclidean propagator pole denominator')
        ga,gb=1/(p+a),1/(p+b)
        normal=0.5*(ga+gb)*np.eye(self.dim,dtype=complex)
        anomalous=0.5*(ga-gb)*self.phase**2*invariant_metric(self.dim)
        return normal,anomalous

    def field_replacements(self, dagger: bool=False) -> dict[str,np.ndarray]:
        """Source spherical scalar = map_R R + map_I I, indexed by gauge slot."""
        common=(self.phase.conjugate()*self.intertwiner.conjugate() if dagger
                else self.phase*self.intertwiner)/np.sqrt(2)
        return {f'{self.species}_R':common,
                f'{self.species}_I':(-1j if dagger else 1j)*common}


def rotate_component_tensor(tensor: np.ndarray, fields: Sequence[str],
                            basis: SelfBilinearMassBasis) -> dict[tuple[str,...],np.ndarray]:
    """Push forward a supplied source CG/vertex tensor, retaining every index.

    Original tensors must be explicitly supplied. This function does not invent
    native Matchete or RGBeta CG coefficients. All input axes are of dimension
    1 or 3 for the selected self-bilinear field; other axes remain untouched.
    """
    fields=tuple(fields)
    t=np.asarray(tensor,dtype=complex)
    if t.ndim!=len(fields) or not fields or len(fields)>4:
        raise ValueError('Source tensor axes must match 1-4 supplied fields')
    opts=[]
    found=False
    for i,f in enumerate(fields):
        if f in (basis.species,basis.species+'dag'):
            if t.shape[i]!=basis.dim:
                raise ValueError('Selected scalar tensor dimension mismatch')
            opts.append(list(basis.field_replacements(f.endswith('dag')).items()))
            found=True
        else:
            opts.append([(f,np.eye(t.shape[i],dtype=complex))])
    if not found:raise ValueError('No selected self-bilinear field in tensor')
    out={}
    for combo in product(*opts):
        names=tuple(name for name,_ in combo)
        value=t
        for axis,(_name,M) in enumerate(combo):
            value=np.tensordot(value,M,axes=([axis],[0]))
            value=np.moveaxis(value,-1,axis)
        out[names]=out.get(names,0)+value
    return out


def evaluate_tensor(tensor: np.ndarray, vectors: Sequence[np.ndarray]) -> complex:
    value=np.asarray(tensor,dtype=complex)
    for vector in vectors:
        value=np.tensordot(value,np.asarray(vector,dtype=complex),axes=(0,0))
    return complex(value)


def parse_matchete_constant(raw: str) -> complex:
    """Read only explicit numeric Mathematica InputForm constants from Stage22.

    This intentionally is NOT a general Mathematica expression interpreter.
    """
    import re
    import sympy as sp
    if not isinstance(raw,str) or not re.fullmatch(r'[0-9ISqrt\s\+\-\*\/\^\(\)\[\]]+',raw):
        raise ValueError('Unsupported Matchete tensor component: '+repr(raw))
    expression=raw.replace('Sqrt[','sqrt(').replace('] ', ') ').replace(']',')').replace('^','**')
    try:
        expr=sp.sympify(expression,locals={'I':sp.I,'sqrt':sp.sqrt})
    except (sp.SympifyError,ValueError,TypeError,SyntaxError) as exc:
        raise ValueError('Invalid numeric Matchete tensor component: '+repr(raw)) from exc
    if expr.free_symbols or expr.is_finite is not True:
        raise ValueError('Nonnumeric Matchete component: '+repr(raw))
    value=complex(sp.N(expr,17))
    return _finite_complex(value,'Matchete tensor component')


def native_tensors_for_model(path: Any, model: str,
                             basis: SelfBilinearMassBasis) -> list[dict[str,Any]]:
    """Read saved Stage22 *extra* quartic source tensors, not ordinary topology CG.

    Source component-order conventions are kept as-is. Identification of their
    complex native Matchete basis with the spherical C convention used here
    is a DECLARED HYPOTHESIS, not an established cross-package map.
    """
    import json
    from pathlib import Path
    obj=json.loads(Path(path).read_text(encoding='utf-8-sig'))
    if (obj.get('format')!='stage22-exact-tensors-v1' or
            obj.get('source')!='MatcheteInvariantTensors' or
            obj.get('model_count')!=16 or obj.get('quartic_slot_count')!=50):
        raise ValueError('Unexpected saved Stage22 Matchete tensor catalogue')
    record=obj.get('models',{}).get(model)
    if not isinstance(record,dict):
        raise ValueError('Missing self-bilinear model from Matchete catalogue')
    data=record.get('tensors')
    expected=6 if model.startswith('E_') else 5
    if not isinstance(data,list) or len(data)!=expected:
        raise ValueError('Wrong number of extra quartics for self-bilinear model')
    observed=set()
    result=[]
    for item in data:
        symbol=item['symbol']
        if not isinstance(symbol,str) or not symbol.startswith('lambdaT3ExtraQ') or symbol in observed:
            raise ValueError('Duplicate/unrecognised extra Matchete quartic')
        observed.add(symbol)
        species=item['species']
        dims=item['dimensions']
        if (not isinstance(species,list) or len(species)!=4 or
            not isinstance(dims,list) or len(dims)!=4 or
            item.get('potential_prefactor')!='1'):
            raise ValueError('Malformed native quartic tensor dimensions/prefactor')
        for s,d in zip(species,dims):
            if not isinstance(d,int) or d not in (1,2,3):
                raise ValueError('Malformed Matchete SU(2) component dimension')
            if s in (basis.species,basis.species+'dag') and d!=basis.dim:
                raise ValueError('Matchete selected field has inconsistent SU(2) dimension')
        raw=item['components']
        if not isinstance(raw,list) or len(raw)!=int(np.prod(dims)):
            raise ValueError('Matchete tensor flattened component count differs')
        tensor=np.array([parse_matchete_constant(x) for x in raw],dtype=complex).reshape(dims)
        result.append({'symbol':symbol,'species':tuple(species),'tensor':tensor,
                       'rotated':rotate_component_tensor(tensor,species,basis)})
    return result
