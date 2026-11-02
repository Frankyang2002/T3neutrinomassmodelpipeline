"""Conditional real-linear numerical Wilson evolution in an explicit basis.

The general one-loop T3 master RGE includes terms linear in conjugate(C):
    16*pi**2 dC/dln(mu) = A(mu) C + B(mu) conjugate(C).
A complex-linear Gamma@C alone cannot in general represent that equation.

This module does NOT derive A/B or certify a physical Wilson normalization.
Both complete, separately supplied sparse A and B matrices are mandatory,
even when one is identically zero. Zero-boundary closure is supported.
"""
from __future__ import annotations

import math
from typing import Any, Mapping

import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import coo_matrix, csr_matrix

LOOP = 1 / (16 * math.pi ** 2)
FORMAT = 't3-explicit-ordered-wilson-real-linear-gamma-v1'
EQUATION = '16pi2*dC/dlnmu=A(mu)@C+B(mu)@conjugate(C)'
BASIS = 'ordered_Matchete_real_components_with_flavour'


def parse_sparse_matrix(payload: Mapping[str, Any], n: int, *, max_nnz: int = 2_000_000) -> csr_matrix:
    """Read a sparse matrix in COO coordinates. No unlisted-row inference."""
    if not isinstance(payload, Mapping) or payload.get('format') != 'complex-coo-v1':
        raise ValueError('Expected explicit complex-coo-v1 matrix')
    if payload.get('shape') != [n, n]:
        raise ValueError('Matrix shape differs from full coordinate basis')
    rows = payload.get('row')
    cols = payload.get('col')
    re = payload.get('real')
    im = payload.get('imag')
    if not all(isinstance(v, list) for v in (rows, cols, re, im)):
        raise ValueError('COO requires row, col, real, imag arrays')
    size = len(rows)
    if size != len(cols) or size != len(re) or size != len(im) or size > max_nnz:
        raise ValueError('Inconsistent or excessive sparse matrix entry count')
    if not all(isinstance(v, int) and not isinstance(v, bool) and 0 <= v < n for v in rows + cols):
        raise ValueError('Sparse indices outside basis')
    if len(set(zip(rows, cols))) != size:
        raise ValueError('Duplicate sparse matrix coordinate')
    real, imag = np.asarray(re, dtype=float), np.asarray(im, dtype=float)
    if not np.isfinite(real).all() or not np.isfinite(imag).all():
        raise ValueError('Nonfinite matrix coefficient')
    return coo_matrix((real + 1j * imag, (rows, cols)), shape=(n, n)).tocsr()


def encode_sparse_matrix(matrix: Any) -> dict:
    """Canonical small/test/sample matrix serializer (not a gamma generator)."""
    coo = csr_matrix(matrix, dtype=np.complex128).tocoo()
    return {'format':'complex-coo-v1', 'shape':list(coo.shape),
            'row':coo.row.tolist(), 'col':coo.col.tolist(),
            'real':coo.data.real.tolist(), 'imag':coo.data.imag.tolist()}


def validate_real_linear_source(gamma: Mapping[str, Any], source_coordinates: list[str],
                                tensor_sha256: str, model: str, mu_f: float, *,
                                max_dim: int = 20000) -> tuple[list[str], np.ndarray, list[csr_matrix], list[csr_matrix]]:
    if gamma.get('format') != FORMAT or gamma.get('equation') != EQUATION:
        raise ValueError('Unsupported real-linear Wilson equation/format')
    if gamma.get('basis') != BASIS or gamma.get('model') != model:
        raise ValueError('Wilson basis or model mismatch')
    if gamma.get('ordered_tensor_sha256') != tensor_sha256:
        raise ValueError('Gamma coefficients not tied to exact boundary SHA-256')
    if gamma.get('closure_includes_zero_boundary_components') is not True:
        raise ValueError('Gamma source must declare its complete closure')
    if gamma.get('conjugate_term_included') is not True:
        raise ValueError('Gamma source must explicitly account for C* contributions')
    if 'mu_F_gev' in gamma and float(gamma['mu_F_gev']) != mu_f:
        raise ValueError('Matching threshold differs from gamma source')
    labels = gamma.get('coordinates')
    if not isinstance(labels, list) or not all(isinstance(x,str) and x for x in labels):
        raise ValueError('Invalid full gamma coordinate inventory')
    if not 1 <= len(labels) <= max_dim or len(labels) != len(set(labels)):
        raise ValueError('Duplicate, empty or oversized gamma coordinate basis')
    if len(set(source_coordinates)) != len(source_coordinates) or not set(source_coordinates).issubset(labels):
        raise ValueError('Gamma basis omits/duplicates nonzero boundary coordinate')
    samples = gamma.get('samples')
    if not isinstance(samples, list) or len(samples) < 2:
        raise ValueError('At least two A/B samples required')
    n = len(labels)
    mus, A, B = [], [], []
    for sample in samples:
        if not isinstance(sample, Mapping) or set(('mu_gev','A','B')) - set(sample):
            raise ValueError('Each sample requires mu_gev, A, B (even if B=0)')
        mu = float(sample['mu_gev'])
        if not math.isfinite(mu) or mu <= 0:
            raise ValueError('Invalid gamma sampling scale')
        mus.append(mu)
        A.append(parse_sparse_matrix(sample['A'], n))
        B.append(parse_sparse_matrix(sample['B'], n))
    scales = np.asarray(mus, dtype=float)
    if not np.all(np.diff(scales) > 0):
        raise ValueError('Gamma scales must be strictly increasing')
    mu_s = float(gamma.get('mu_S_gev', float('nan')))
    if not math.isfinite(mu_s) or not 0 < mu_s < mu_f:
        raise ValueError('Require 0 < mu_S < mu_F')
    if not (scales[0] <= mu_s and mu_f <= scales[-1]):
        raise ValueError('A/B samples do not cover both thresholds')
    return labels, np.log(scales), A, B


def evolve_real_linear(boundary_coords: list[str], boundary_values: np.ndarray,
                       gamma: Mapping[str, Any], *, model: str,
                       tensor_sha256: str, mu_f: float,
                       points: int = 41, max_dim: int = 20000,
                       rtol: float = 1e-9, atol: float = 1e-12) -> dict:
    """Evolve a *conditional* RGE with explicitly supplied A, B, closure basis."""
    if not isinstance(points, int) or not 2 <= points <= 501:
        raise ValueError('Require 2..501 samples')
    boundary = np.asarray(boundary_values, dtype=np.complex128)
    if boundary.shape != (len(boundary_coords),) or not np.isfinite(boundary).all():
        raise ValueError('Invalid ordered boundary vector')
    labels, log_mu, A, B = validate_real_linear_source(
        gamma, boundary_coords, tensor_sha256, model, mu_f, max_dim=max_dim)
    position = {x:i for i,x in enumerate(labels)}
    initial = np.zeros(len(labels), dtype=np.complex128)
    for label, val in zip(boundary_coords, boundary):
        initial[position[label]] = val

    def derivative(t: float, z: np.ndarray) -> np.ndarray:
        high = min(max(int(np.searchsorted(log_mu, t, side='right')), 1), len(log_mu)-1)
        low = high-1
        frac = (t-log_mu[low]) / (log_mu[high]-log_mu[low])
        linear = ((1-frac)*(A[low] @ z) + frac*(A[high] @ z))
        antilinear = ((1-frac)*(B[low] @ z.conjugate()) + frac*(B[high] @ z.conjugate()))
        return LOOP * (linear + antilinear)

    t0 = math.log(mu_f)
    t1 = math.log(float(gamma['mu_S_gev']))
    sol = solve_ivp(derivative, (t0, t1), initial,
                    t_eval=np.linspace(t0,t1,points), rtol=rtol, atol=atol,
                    method='DOP853')
    if not sol.success or sol.y.shape != (len(labels),points) or not np.isfinite(sol.y).all():
        raise RuntimeError('Real-linear Wilson ODE failed: '+sol.message)
    end = sol.y[:, -1]
    return {'format':'t3-conditional-real-linear-ordered-wilson-trajectory-v1',
            'status':'ConditionalRealLinearGammaTransport', 'model':model,
            'equation':EQUATION, 'ordered_tensor_sha256':tensor_sha256,
            'mu_F_gev':mu_f, 'mu_S_gev':float(gamma['mu_S_gev']),
            'coordinates':labels, 'source_coordinate_count':len(boundary_coords),
            'zero_initial_closure_coordinates':len(labels)-len(boundary_coords),
            'scales_gev':np.exp(sol.t).tolist(),
            'trajectory_real':sol.y.T.real.tolist(),
            'trajectory_imag':sol.y.T.imag.tolist(),
            'initial_norm_gev_inv':float(np.linalg.norm(initial)),
            'end_norm_gev_inv':float(np.linalg.norm(end)),
            'relative_change':float(np.linalg.norm(end-initial)/max(1e-300,np.linalg.norm(initial))),
            'nfev':int(sol.nfev),
            'anomalous_dimension_derived_by_this_module':False,
            'canonical_C_ijab_normalization_certified':False,
            'physical_eft1_wilson_rge_certified':False,
            'scalar_threshold_matching_performed':False}
