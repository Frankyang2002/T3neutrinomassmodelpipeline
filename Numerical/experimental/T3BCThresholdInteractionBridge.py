"""B/C running mixed-doublet mass-basis diagnostic, using original source code.

The source snapshot is the *candidate RGBeta* UV/EFT1 full-mass trajectory.
A numerical identity between scalar quadratic forms is checked; it does NOT
identify a cross-package physical quartic/Yukawa spinor normalization and
never produces a Weinberg matching coefficient.
"""
from __future__ import annotations

import math
from typing import Any, Mapping

import numpy as np


def extract_complex(x: Mapping[str, float]) -> complex:
    if not isinstance(x, dict) or set(x) != {'real', 'imag'}:
        raise ValueError('Expected native real/imag scalar')
    z = complex(x['real'], x['imag'])
    if not (math.isfinite(z.real) and math.isfinite(z.imag)):
        raise ValueError('Nonfinite mass-mixing coupling')
    return z


def bridge_snapshot(model: str, row: Mapping[str, Any], *, trials: int = 8) -> dict:
    from Numerical.experimental.GeneralT3MixedScalarMassBasis import (
        MixedDoubletMassBasis, verify_tensor_transport, quadratic_mass_tensor_transport,
        verify_su2_covariance,
    )
    if model not in ('B_-1', 'C_-1'):
        raise ValueError('Only B_-1/C_-1 are mixed doublets')
    # The source snapshot stores mS1Sq/mS2Sq at its top level. The nested
    # 'state' is a separately serialized T3UVState/IntermediateScalarState.
    mu = extract_complex(row['muT3GeneralS1S2_gev2'])
    basis = MixedDoubletMassBasis.build(model, float(row['mS1Sq_gev2']),
                                        float(row['mS2Sq_gev2']), mu)
    if basis.tachyonic or basis.degenerate:
        raise ValueError('Tachyonic/degenerate scalar eigenmass at running snapshot')
    rng = np.random.default_rng(20261101)
    max_res = 0.0
    for _ in range(trials):
        x = rng.normal(size=2) + 1j*rng.normal(size=2)
        y = rng.normal(size=2) + 1j*rng.normal(size=2)
        s1, s2 = basis.original_fields(x, y)
        rx, ry = basis.eigen_fields(s1, s2)
        if not (np.allclose(rx, x, rtol=1e-12, atol=1e-12) and
                np.allclose(ry, y, rtol=1e-12, atol=1e-12)):
            raise ValueError('Old/new mixed scalar field map not invertible')
        a = basis.old_potential(s1, s2)
        b = basis.new_potential(x, y)
        res = abs(a-b)/max(abs(a), abs(b), 1)
        max_res = max(max_res, res)
        if res > 1e-10:
            raise ValueError('Original/rotated quadratic scalar action inconsistent')
    # Recheck the actual source's full holomorphic bilinear + Hermitian
    # conjugate transport, not just eigenvalues of the 2x2 matrix.
    quadratic_tensor = quadratic_mass_tensor_transport(basis)
    su2_covariance = float(verify_su2_covariance(basis, trials=8))
    tensor_fields = (
        ('H','H','S1','S2dag'),
        ('Hdag','Hdag','S2','S2'),
        ('Hdag','Hdag','S1dag','S1dag'),
        ('S1dag','S2','S2dag','S2dag'),
        ('S1dag','S2dag'),
    )
    tensor_residuals = {'/'.join(fields):float(verify_tensor_transport(fields,basis.U,seed=20261101,trials=4))
                        for fields in tensor_fields}
    result = basis.info()
    mass_sq = np.asarray(result['masses_squared_gev2'], dtype=float)
    saved = row['scalar_spectrum']
    if 'mass_squared_eigenvalues_gev2' in saved:
        full = np.asarray(saved['mass_squared_eigenvalues_gev2'], dtype=float)
        if full.shape != (8,) or not np.allclose(
            np.sort(full), np.sort(np.repeat(mass_sq, 4)), rtol=1e-10, atol=1e-7):
            raise ValueError('Running scalar Hessian spectrum disagrees with mixed doublet')
    return {'model': model, 'mu_gev': float(row['mu_gev']),
            'mass_basis': result,
            'quadratic_action_relative_residual_max': max_res,
            'source_bilinear_tensor_transport': quadratic_tensor,
            'su2_covariance_residual': su2_covariance,
            'generic_tensor_rotation_residuals': tensor_residuals,
            'valid_running_scalar_eigenbasis': True,
            'rotated_full_scalar_quartic_yukawa_matching_certified': False,
            'matchete_rgbeta_component_map_certified': False,
            'physical_matching_performed': False}


def bridge_trajectory(model: str, report: Mapping[str, Any], *, stride: int = 4) -> dict:
    if report.get('model') != model or report.get('solver_success') is not True:
        raise ValueError('Invalid B/C coupled UV/EFT1 running report')
    output = {}
    for name in ('uv_trajectory', 'eft1_trajectory'):
        rows = report[name]
        if not isinstance(rows, list) or not rows:
            raise ValueError('Missing running scalar trajectory')
        indices = sorted(set([0, len(rows)-1] + list(range(0, len(rows), max(1, stride)))))
        output[name] = [bridge_snapshot(model, rows[i]) for i in indices]
    uend = output['uv_trajectory'][-1]['mass_basis']['masses_squared_gev2']
    estart = output['eft1_trajectory'][0]['mass_basis']['masses_squared_gev2']
    if not np.allclose(uend, estart, rtol=5e-11, atol=1e-7):
        raise ValueError('Unexplained scalar mass discontinuity in tree-level F handoff')
    return {'model': model, 'status': 'ConditionalRGBetaMixedDoubletBasisChecked',
            'uv_trajectory': output['uv_trajectory'],
            'eft1_trajectory': output['eft1_trajectory'],
            'scalar_mass_continuity_at_grouped_F': True,
            'physical_sequential_matching_certified': False}
