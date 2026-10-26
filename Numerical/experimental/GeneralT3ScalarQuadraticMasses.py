"""Experimental general T3 scalar quadratic mass-mixing Hessians.

Models use a specified *spherical SU(2) basis*, NOT a certified Matchete
component convention. Generated Hessians are a useful numerical cross-check
of which previously omitted quadratic interactions affect the scalar masses.
They are NOT automatically inserted into RGBeta or the matching thresholds.

Quadratic convention:
  V2 = mS1^2 S1^dag S1 + mS2^2 S2^dag S2
       + [mu2_extra (S_a^dag)^T C (S_b^dag) + h.c.].

The self-field term has no implicit factor of 1/2, matching the independent
PlusHc[mu2_extra * product] structure in GeneralScalarOperators.wl.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Mapping, Any
import numpy as np

SPECIES=('H','Hdag','S1','S1dag','S2','S2dag')


def _metric(dim:int) -> np.ndarray:
    if dim==1:
        return np.eye(1,dtype=float)
    if dim==2:
        return np.array([[0,1],[-1,0]],dtype=float)
    if dim==3:
        # j=1 singlet in {|m=+1>,|m=0>,|m=-1>} spherical basis.
        return np.array([[0,0,1],[0,-1,0],[1,0,0]],dtype=float)
    raise ValueError('Only singlet, doublet and triplet currently supported.')


def quadratic_structure(manifest:Mapping[str,Any],model:str) -> list[dict]:
    record=manifest['models'].get(model)
    if record is None:raise ValueError(f'Unknown model {model}')
    return list(record['quadratic_structures_not_enabled'])


def quadratic_mass_hessian(manifest:Mapping[str,Any],model:str,
                           m_s1_sq:float,m_s2_sq:float,extra_mu_sq:complex=0j) -> dict:
    """Real 2(dS1+dS2)x2(dS1+dS2) Hessian and mass-squared eigenvalues."""
    from scripts.GenerateGeneralT3RGBetaBasis import REP
    point=manifest['models'].get(model)
    if point is None:raise ValueError(f'Unknown model {model}')
    family=model.split('_')[0]
    d1,d2,_=REP[family]
    masses=(float(m_s1_sq),float(m_s2_sq))
    if not all(np.isfinite(m) for m in masses):raise ValueError('Nonfinite input mass-squared')
    mu=complex(extra_mu_sq)
    if not (np.isfinite(mu.real) and np.isfinite(mu.imag)):
        raise ValueError('Nonfinite quadratic mixing coefficient')
    n=d1+d2
    M=np.diag([masses[0]]*d1+[masses[1]]*d2).astype(complex)
    B=np.zeros((n,n),dtype=complex)
    structures=quadratic_structure(manifest,model)
    if len(structures)>1:
        raise NotImplementedError('More than one extra quadratic invariant needs a multivariate mass Hessian.')
    if structures:
        counts=tuple(structures[0]['counts'])
        names=[s for s,c in zip(SPECIES,counts) for _ in range(c)]
        if len(names)!=2 or any(n not in ('S1dag','S2dag') for n in names):
            raise NotImplementedError('Unsupported non-holomorphic quadratic structure.')
        if names==['S1dag','S1dag']:
            B[:d1,:d1]=2*mu.conjugate()*_metric(d1)
        elif names==['S2dag','S2dag']:
            B[d1:,d1:]=2*mu.conjugate()*_metric(d2)
        elif names==['S1dag','S2dag']:
            if d1!=d2:raise ValueError('Off-diagonal SU(2) invariant requires equal multiplet sizes.')
            B[:d1,d1:]=mu.conjugate()*_metric(d1)
            B[d1:,:d1]=B[:d1,d1:].T
        else:
            raise NotImplementedError(f'Quadratic structure {names} not implemented.')
    elif abs(mu)>0:
        raise ValueError('Nonzero extra quadratic supplied to model with no allowed extra quadratic.')
    if not np.allclose(B,B.T,rtol=0,atol=1e-12):
        raise RuntimeError('Quadratic scalar B matrix not symmetric.')
    # q=(Re S1,Re S2,Im S1,Im S2), S=(x+i y)/sqrt(2).
    hessian=np.block([[M.real+B.real, -B.imag],[-B.imag, M.real-B.real]])
    if not np.allclose(hessian,hessian.T,rtol=0,atol=1e-12):
        raise RuntimeError('Real scalar mass Hessian not symmetric.')
    values=np.linalg.eigvalsh(hessian)
    return {'status':'SphericalBasisQuadraticPrototype',
            'model':model,'additional_quadratic_structures':len(structures),
            'dimensions':[d1,d2], 'extra_mass_coefficient_gev2':{'real':mu.real,'imag':mu.imag},
            'mass_squared_eigenvalues_gev2':values.tolist(),
            'tachyonic':bool(np.any(values<0)),
            'real_scalar_hessian_gev2':hessian.tolist(),
            'matchete_mass_basis_verified':False,
            'included_in_threshold_matching':False}


def load_manifest(path:Path|str):
    return json.loads(Path(path).read_text(encoding='utf-8'))
