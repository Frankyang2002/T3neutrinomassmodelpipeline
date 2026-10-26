"""Evaluate source-derived B/C(-1) extra C5 with RGBeta couplings mapped by -1.

Only a conditional *UV-endpoint/frozen-input* sensitivity estimate. The two
extra Matchete couplings are -lambda_RGBeta (not complex-conjugated), as
derived from the explicit scalar-Lagrangian convention for common field basis.
"""
from __future__ import annotations

from typing import Any, Mapping

import numpy as np

from Numerical.experimental.GeneralT3C5Pilot import ResponseKernel, evaluate_response
from Numerical.experimental.GeneralT3ConditionalUVBridge import (
    RGBETA_BY_MAT,
    _encode,
    _induced_rgbeta_couplings,
    _validate_reference,
    _validate_tensor_mapping,
)
from Numerical.experimental.GeneralT3LagrangianConvention import MAPPING_FIELDS


def source_mapped_uv_to_c5(
    model: str,
    parameters: Mapping[str, Any],
    stage13: Mapping[str, Any],
    tensor_map: Mapping[str, Any],
    convention: Mapping[str, Any],
    kernels: tuple[ResponseKernel, ...],
    reference: Mapping[str, Any],
    *,
    stage14: Mapping[str, Any] | None = None,
    vev_gev: float = 246.22,
) -> dict[str, Any]:
    if model not in {'B_m1', 'C_m1'}:
        raise ValueError('Only B/C(-1) models are supported')
    _validate_tensor_mapping(tensor_map, model)
    if convention.get('status') != 'SourceLagrangianConventionDerived' or convention.get('source_checks_passed') is not True:
        raise ValueError('Source-derived convention is not certified')
    if convention.get('full_threshold_matching_verified') is not False:
        raise ValueError('Cannot treat this as validated threshold matching')
    if convention.get('same_field_basis_assumed') is not True:
        raise ValueError('Expected same-field-basis convention')
    actual = [r for r in convention.get('results', []) if r.get('model') == model]
    if len(actual) != 2:
        raise ValueError('Expected two source-derived convention entries')
    for entry in actual:
        mat = entry.get('matchete_coupling')
        if mat not in MAPPING_FIELDS or entry.get('rgbeta_coupling') != MAPPING_FIELDS[mat][0]:
            raise ValueError('Unrecognised convention coupling pair')
        if entry.get('matchete_in_terms_of_rgbeta') != 'negative' or entry.get('relative_sign') != -1:
            raise ValueError('Unexpected convention sign')
        if entry.get('requires_complex_conjugation') is not False:
            raise ValueError('Unexpected complex conjugation in source mapping')
    if set(row['matchete_coupling'] for row in actual) != set(RGBETA_BY_MAT):
        raise ValueError('Source convention incomplete')
    induced = _induced_rgbeta_couplings(stage13, model, parameters)
    reference_norm = _validate_reference(reference, model)
    vev = float(vev_gev)
    if not np.isfinite(vev) or vev <= 0:
        raise ValueError('Invalid Higgs vev')
    if set(item.coupling for item in kernels) != set(RGBETA_BY_MAT) or len(kernels) != 2:
        raise ValueError('Exactly two formal Matchete response kernels required')
    params = dict(parameters)
    params['matchete_extra_couplings'] = {
        mat: _encode(-induced[rgb]) for mat, rgb in RGBETA_BY_MAT.items()
    }
    evaluated = evaluate_response(kernels, params)
    c5 = np.asarray(evaluated['delta_c5']['real'], dtype=float) + 1j * np.asarray(evaluated['delta_c5']['imag'], dtype=float)
    if c5.shape != (3, 3) or not np.isfinite(c5).all() or not np.allclose(c5, c5.T, atol=1e-22, rtol=1e-10):
        raise ValueError('Evaluated extra C5 must be finite and symmetric')
    mnu = -(vev ** 2 / 2) * 1e9 * c5
    mass_norm = float(np.linalg.norm(mnu, 'fro'))
    singular = sorted(float(s) for s in np.linalg.svd(mnu, compute_uv=False))
    if stage14 is not None:
        if stage14.get('model') != model or stage14.get('status') != 'ConditionalUVInducedC5SensitivityOnly':
            raise ValueError('Stage 14 model/status mismatch')
        cases = [row for row in stage14.get('mapping_cases', [])
                 if all(row.get('mapping_choices', {}).get(mat) == 'negative' for mat in RGBETA_BY_MAT)]
        if len(cases) != 1:
            raise ValueError('Stage 14 did not contain exactly one all-negative scenario')
        reference_case = cases[0]
        if not np.isclose(mass_norm, float(reference_case['extra_only_mnu_norm_ev']), rtol=1e-9, atol=1e-27):
            raise ValueError('Source-derived result differs from matching Stage 14 scenario')
        if not np.allclose(singular, reference_case['extra_only_takagi_singular_values_ev'], rtol=1e-8, atol=1e-24):
            raise ValueError('Stage 14 singular values inconsistent with source-derived result')
    return {
        'status': 'SourceMappedConditionalUVInducedC5',
        'model': model,
        'source_convention': 'lambda_matchete = - lambda_rgbeta, same field word, no conjugation',
        'mapped_matchete_extra_couplings': params['matchete_extra_couplings'],
        'induced_rgbeta_quartics': {k: _encode(v) for k, v in induced.items()},
        'delta_c5': evaluated['delta_c5'],
        'delta_mnu_ev': {
            'real': mnu.real.tolist(),
            'imag': mnu.imag.tolist(),
            'frobenius_norm': mass_norm,
            'takagi_singular_values_ev': singular,
        },
        'reference_stage11_c5_norm_gev_inverse': reference_norm,
        'relative_to_stage11_extra_only_c5_norm': float(evaluated['delta_c5']['frobenius_norm']) / reference_norm,
        'stage14_negative_case_reproduced': bool(stage14 is not None),
        'uv_only_endpoint_gev': float(stage13['endpoint_gev']),
        'scale_crosses_heavy_threshold': False,
        'same_field_basis_assumed': True,
        'threshold_matching': False,
        'full_c5_running': False,
        'physical_neutrino_spectrum': False,
        'disclaimer': (
            'Source-convention-mapped UV-endpoint sensitivity only. The quartics are radiatively generated '
            'at one loop, so one-loop C5 insertion is formally O(hbar^2). The toy Yukawas and masses '
            'are frozen; heavy/scalar threshold matching and the original C5 are not included.'
        ),
    }
