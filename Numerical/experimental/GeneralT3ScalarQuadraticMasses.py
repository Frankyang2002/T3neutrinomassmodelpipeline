"""Candidate quadratic scalar mass Hessians for the 16 ordinary T3 models.

This is a fixed spherical SU(2) component convention, NOT a certified
Matchete<->RGBeta Cartesian field identification. It extends the existing
prototype without changing the callable signature or production thresholds.

V2 = mS1Sq S1†S1 + mS2Sq S2†S2
     + [mu2_extra (Sa†)^T C (Sb†) + h.c.].

There is NO hidden 1/2 in the holomorphic contribution. Accordingly,
self-pair mass splittings are 2*abs(mu2_extra), while mixed-pair mass
off-diagonal terms are abs(mu2_extra). For SU(2) spin j=1, the invariant
metric is the spherical-basis matrix C_(m,n)=(-1)^(1-m) delta_(m,-n).
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Mapping, Any

import numpy as np

SPECIES = ('H', 'Hdag', 'S1', 'S1dag', 'S2', 'S2dag')
REP = {'A': (1,3,2), 'B': (2,2,1), 'C': (2,2,3),
       'D': (3,1,2), 'E': (3,3,2)}


def _metric(dim: int) -> np.ndarray:
    """Invariant bilinear in the declared spherical isospin basis."""
    if dim == 1:
        return np.eye(1, dtype=float)
    if dim == 2:
        return np.array([[0., 1.], [-1., 0.]])
    if dim == 3:
        return np.array([[0., 0., 1.], [0., -1., 0.], [1., 0., 0.]])
    raise ValueError('Only SU(2) dimensions 1, 2 and 3 are supported')


def quadratic_structure(manifest: Mapping[str, Any], model: str) -> list[dict]:
    records = manifest.get('models')
    if not isinstance(records, dict) or model not in records:
        raise ValueError(f'Unknown quadratic T3 model: {model}')
    structures = records[model].get('quadratic_structures_not_enabled')
    if not isinstance(structures, list) or len(structures) > 1:
        raise ValueError(f'Unexpected quadratic multiplicity for {model}')
    if records[model].get('quadratic_complex_couplings_not_enabled') != len(structures):
        raise ValueError(f'Quadratic registry mismatch for {model}')
    return structures


def quadratic_kind(manifest: Mapping[str, Any], model: str) -> str:
    structures = quadratic_structure(manifest, model)
    if not structures:
        return 'none'
    counts = structures[0].get('counts')
    if not isinstance(counts, list) or len(counts) != 6 or sum(counts) != 2:
        raise ValueError(f'Malformed quadratic counts: {model}')
    names = [name for name, count in zip(SPECIES, counts) for _ in range(count)]
    kinds = {('S1dag', 'S1dag'): 'self_s1',
             ('S2dag', 'S2dag'): 'self_s2',
             ('S1dag', 'S2dag'): 'mixed'}
    kind = kinds.get(tuple(names))
    if kind is None or structures[0].get('invariants') != 1:
        raise ValueError(f'Unrecognised quadratic invariant: {model}/{names}')
    family, alpha_s = model.split('_')
    if family not in REP:
        raise ValueError(f'Unrecognised T3 representation: {model}')
    d1, d2, _ = REP[family]
    alpha = int(alpha_s)
    if kind == 'self_s1' and not (alpha == 0 and d1 % 2 == 1):
        raise ValueError('S1 holomorphic invariant requires a neutral integer-isospin scalar')
    if kind == 'self_s2' and not (alpha == -2 and d2 % 2 == 1):
        raise ValueError('S2 holomorphic invariant requires a neutral integer-isospin scalar')
    if kind == 'mixed' and not (alpha == -1 and d1 == d2 == 2):
        raise ValueError('Mixed bilinear quantum-number inconsistency')
    return kind


def analytic_mass_squared(manifest: Mapping[str, Any], model: str,
                          m_s1_sq: float, m_s2_sq: float,
                          extra_mu_sq: complex = 0j) -> np.ndarray:
    """Independent analytic spectrum of the gauge-unbroken quadratic potential."""
    family = model.split('_')[0]
    d1, d2, _ = REP[family]
    kind = quadratic_kind(manifest, model)
    m1, m2 = float(m_s1_sq), float(m_s2_sq)
    mu = complex(extra_mu_sq)
    if not all(np.isfinite(x) for x in (m1, m2, mu.real, mu.imag)):
        raise ValueError('Nonfinite quadratic inputs')
    if kind == 'none':
        if mu != 0:
            raise ValueError(f'{model} has no independent extra quadratic')
        values = [m1] * (2*d1) + [m2] * (2*d2)
    elif kind == 'self_s1':
        values = ([m1-2*abs(mu)]*d1 + [m1+2*abs(mu)]*d1
                  + [m2]*(2*d2))
    elif kind == 'self_s2':
        values = ([m1]*(2*d1) + [m2-2*abs(mu)]*d2
                  + [m2+2*abs(mu)]*d2)
    else:
        middle = 0.5*(m1+m2)
        split = np.hypot(0.5*(m1-m2), abs(mu))
        values = ([middle-split]*(2*d1) + [middle+split]*(2*d1))
    return np.sort(np.asarray(values, dtype=float))


def quadratic_mass_hessian(manifest: Mapping[str, Any], model: str,
                           m_s1_sq: float, m_s2_sq: float,
                           extra_mu_sq: complex = 0j) -> dict[str, Any]:
    """Construct 2(d1+d2) Cartesian-real Hessian; test independent spectrum."""
    if model.split('_')[0] not in REP:
        raise ValueError('Unknown ordinary T3 model')
    d1, d2, _ = REP[model.split('_')[0]]
    kind = quadratic_kind(manifest, model)
    m1, m2 = float(m_s1_sq), float(m_s2_sq)
    mu = complex(extra_mu_sq)
    if not all(np.isfinite(x) for x in (m1, m2, mu.real, mu.imag)):
        raise ValueError('Nonfinite quadratic mass input')
    if kind == 'none' and mu != 0:
        raise ValueError(f'{model}: nonzero mu2 has no allowed operator')
    n = d1+d2
    A = np.diag([m1]*d1 + [m2]*d2).astype(complex)
    B = np.zeros((n,n),dtype=complex)
    if kind == 'self_s1':
        B[:d1,:d1] = 2*mu.conjugate()*_metric(d1)
    elif kind == 'self_s2':
        B[d1:,d1:] = 2*mu.conjugate()*_metric(d2)
    elif kind == 'mixed':
        if d1 != d2:
            raise ValueError('Mixed bilinear requires equal SU(2) dimensions')
        B[:d1,d1:] = mu.conjugate()*_metric(d1)
        B[d1:,:d1] = B[:d1,d1:].T
    if not np.allclose(B, B.T, atol=1e-13, rtol=0):
        raise RuntimeError('Holomorphic mass tensor is not symmetric')
    H = np.block([[A.real+B.real, -B.imag], [-B.imag, A.real-B.real]])
    if not np.allclose(H, H.T, atol=1e-13, rtol=0):
        raise RuntimeError('Cartesian mass Hessian is not symmetric')
    eigen = np.linalg.eigvalsh(H)
    predicted = analytic_mass_squared(manifest, model, m1, m2, mu)
    # Check the entire spectrum, not just the trace or average eigenvalue.
    if not np.allclose(eigen, predicted, rtol=2e-12, atol=2e-9):
        raise RuntimeError(f'{model}: numerical Hessian does not match analytic spectrum')
    return {
        'status': 'SphericalBasisQuadraticPrototype',
        'model': model,
        'quadratic_kind': kind,
        'additional_quadratic_structures': int(kind != 'none'),
        'dimensions': [d1, d2],
        'extra_mass_coefficient_gev2': {'real':mu.real, 'imag':mu.imag},
        'mass_squared_eigenvalues_gev2': eigen.tolist(),
        'analytic_mass_squared_eigenvalues_gev2': predicted.tolist(),
        'analytic_spectrum_verified': True,
        'minimum_mass_squared_gev2': float(eigen[0]),
        'tachyonic': bool(np.any(eigen < 0)),
        'real_scalar_hessian_gev2': H.tolist(),
        'su2_component_basis': 'chosen spherical m=j,...,-j',
        'quadratic_potential_prefactor': 'mu2 * invariant + h.c.; no 1/2',
        'matchete_mass_basis_verified': False,
        'rgbeta_mass_normalisation_verified': False,
        'included_in_threshold_matching': False,
    }


def load_manifest(path: Path | str) -> dict[str, Any]:
    result = json.loads(Path(path).read_text(encoding='utf-8'))
    if not isinstance(result, dict) or result.get('model_count') != 16:
        raise ValueError('Incomplete quadratic manifest')
    if set(result.get('models',{})) != {
            'A_-4','A_-2','A_+0','B_-3','B_-1','B_+1',
            'C_-3','C_-1','C_+1','D_-2','D_+0','D_+2',
            'E_-4','E_-2','E_+0','E_+2'}:
        raise ValueError('Wrong ordinary T3 manifest model set')
    if sum(len(quadratic_structure(result,name)) for name in result['models']) != 8:
        raise ValueError('Expected exactly eight catalogued independent quadratic slots')
    return result
