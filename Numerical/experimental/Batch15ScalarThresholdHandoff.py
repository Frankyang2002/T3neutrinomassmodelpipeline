"""Batch 15: lossless renormalisable EFT1 -> scalar-threshold pre-matching handoff.

Project path: Numerical/experimental/Batch15ScalarThresholdHandoff.py

This is an audit of the Batch13/14 RGBeta-native candidate pipeline. It is NOT
an EFT1 Wilson boundary, a two-scalar one-loop match, or a physical C5 result.
The quadratic terms can split or mix real scalar components; diagonal masses
are only diagnostic proxies and are not used as physical eigenmasses.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from Numerical.experimental.Batch13GeneralUVRunning import validate_uv_export
from Numerical.experimental.Batch14EFT1Running import (
    load_saved_uv_endpoint, validate_eft1_export,
)
from Numerical.running.FermionThresholdBoundary import (
    build_intermediate_scalar_boundary, fermion_threshold_masses,
)
from Numerical.running.IntermediateScalarStateVector import (
    pack_intermediate_scalar_state, unpack_intermediate_scalar_state,
)
from Numerical.running.IntermediateScalarRGBetaEvaluator import (
    intermediate_scalar_state_environment,
)
from Numerical.running.ScalarThresholdBoundary import build_final_sm_boundary

GENERAL_MASS_KEYS = frozenset((
    'muT3GeneralS1Sq', 'muT3GeneralS2Sq', 'muT3GeneralS1S2',
))
FERMION_KEYS = frozenset(('MF', 'y1', 'y2', 'y1Conj', 'y2Conj'))


def read_json(path: Path | str) -> dict[str, Any]:
    obj = json.loads(Path(path).read_text(encoding='utf-8-sig'))
    if not isinstance(obj, dict):
        raise ValueError(f'{path}: expected a JSON object')
    return obj


def sha256_file(path: Path | str) -> str:
    h = hashlib.sha256()
    with Path(path).open('rb') as file:
        for data in iter(lambda: file.read(1024 * 1024), b''):
            h.update(data)
    return h.hexdigest()


def _isclose_complex(a: Any, b: Any, *, rtol: float = 5e-10,
                     atol: float = 1e-14) -> bool:
    return bool(np.allclose(np.asarray(a, dtype=complex),
                            np.asarray(b, dtype=complex),
                            rtol=rtol, atol=atol))


def _decode_value(raw: Any, block) -> np.ndarray | complex | float:
    """Decode exactly the complex-array format used in Batches 13 and 14."""
    if not isinstance(raw, Mapping) or set(raw) != {'real', 'imag'}:
        raise ValueError(f'{block.name}: expected real/imag encoded value')
    real, imag = np.asarray(raw['real'], float), np.asarray(raw['imag'], float)
    if real.shape != block.shape or imag.shape != block.shape:
        raise ValueError(f'{block.name}: incorrect value shape; expected {block.shape}')
    if not (np.isfinite(real).all() and np.isfinite(imag).all()):
        raise ValueError(f'{block.name}: nonfinite encoded value')
    if block.kind == 'real_scalar':
        if float(imag) != 0:
            raise ValueError(f'{block.name}: complex value in real state slot')
        return float(real)
    if block.kind == 'complex_scalar':
        return complex(float(real), float(imag))
    if block.kind == 'complex_matrix':
        return real + 1j * imag
    raise ValueError(f'{block.name}: unrecognised canonical block type')


def decode_eft1_endpoint(raw: Mapping[str, Any], template, *, mu_gev: float):
    """Reconstruct an EFT1 state, including all additional quartic/quadratic slots."""
    if not isinstance(raw, Mapping):
        raise ValueError('EFT1 endpoint needs a coupling object')
    _, layout = pack_intermediate_scalar_state(template)
    expected = set(layout.names)
    if set(raw) != expected:
        raise ValueError('EFT1 endpoint coupling mismatch: '
                         f'missing={sorted(expected-set(raw))}, '
                         f'extra={sorted(set(raw)-expected)}')
    vec = np.zeros(layout.size, dtype=float)
    for block in layout.blocks:
        value = _decode_value(raw[block.name], block)
        if block.kind == 'real_scalar':
            vec[block.start] = float(value)
        elif block.kind == 'complex_scalar':
            z = complex(value)
            vec[block.start:block.stop] = (z.real, z.imag)
        else:
            m = np.asarray(value, dtype=complex)
            n = m.size
            vec[block.start:block.start+n] = m.real.ravel()
            vec[block.start+n:block.stop] = m.imag.ravel()
    return unpack_intermediate_scalar_state(vec, layout, mu_gev=float(mu_gev))


def compare_states(a, b, *, label: str) -> None:
    ea, eb = intermediate_scalar_state_environment(a), intermediate_scalar_state_environment(b)
    if set(ea) != set(eb):
        raise ValueError(f'{label}: coupling inventory differs')
    bad = [k for k in sorted(ea) if not _isclose_complex(ea[k], eb[k])]
    if bad:
        raise ValueError(f'{label}: differing couplings: {bad}')


def scalar_threshold_audit(state, *, mu_f_gev: float, mu_s_gev: float) -> dict[str, Any]:
    """Report bare diagonal masses; never identify them as eigenmasses with mu2."""
    m1sq = float(state.mS1Sq)
    m2sq = float(state.mS2Sq)
    if not all(math.isfinite(x) for x in (m1sq, m2sq)):
        raise ValueError('Nonfinite diagonal scalar mass-square coefficient')
    if min(m1sq, m2sq) <= 0:
        raise ValueError('Nonpositive diagonal scalar mass-square coefficient')
    proxy = (math.sqrt(m1sq), math.sqrt(m2sq))
    if not (0 < mu_s_gev < mu_f_gev):
        raise ValueError('Invalid scale ordering muF > muS > 0')
    quadratic = {k: complex(v) for k, v in state.general_couplings.items()
                 if k in GENERAL_MASS_KEYS}
    if len(quadratic) > 1:
        raise ValueError('Unexpected multiple independent scalar bilinears')
    for k, v in quadratic.items():
        if not (math.isfinite(v.real) and math.isfinite(v.imag)):
            raise ValueError(f'Nonfinite {k}')
    no_mass_mixing = not quadratic
    # Even if a candidate extra bilinear happens to be zero at this scale,
    # its normalization / mass-matrix mapping is not source certified.
    hierarchy_proxy = bool(min(proxy) > mu_s_gev and max(proxy) < mu_f_gev)
    return {
        'diagonal_mass_squared_gev2': {'mS1Sq': m1sq, 'mS2Sq': m2sq},
        'diagonal_mass_proxies_gev': list(proxy),
        'diagonal_proxy_threshold_hierarchy_passed': hierarchy_proxy,
        'additional_quadratic_symbols': sorted(quadratic),
        'additional_quadratic_values': {
            k: {'real': v.real, 'imag': v.imag} for k, v in quadratic.items()},
        'physical_scalar_mass_eigenvalues_computed': False,
        'physical_scalar_mass_basis_certified': False,
        'mass_mixing_requires_resolution': not no_mass_mixing,
        'simultaneous_two_scalar_matching_performed': False,
    }


def handoff(model: str, uv_payload: Mapping[str,Any], eft1_payload: Mapping[str,Any],
            uv_endpoint_path: Path, eft1_endpoint: Mapping[str,Any]) -> dict[str,Any]:
    """Rebuild both threshold states and certify lossless renormalisable handoff.

    Explicitly decline to construct the C5/Wilson boundary: that requires the
    appropriate physical matching calculation, not only an RGBeta trajectory.
    """
    validate_uv_export(uv_payload, model)
    names = validate_eft1_export(eft1_payload, uv_payload, model)
    if set(names) & FERMION_KEYS:
        raise ValueError('Heavy fermion parameter leaked into EFT1')
    if eft1_endpoint.get('model') != model or eft1_endpoint.get('status') != 'CandidateFThresholdAndEFT1Integrated':
        raise ValueError(f'{model}: invalid Batch14 EFT1 endpoint status/model')
    for flag in ('fermion_threshold_wilson_matching_performed',
                 'general_scalar_mass_matrix_diagonalised',
                 'matchete_rgbeta_normalisation_certified',
                 'physical_neutrino_mass_prediction'):
        if eft1_endpoint.get(flag) is not False:
            raise ValueError(f'{model}: unexpected physical certification flag {flag}')
    if eft1_endpoint.get('both_scalars_active_in_eft1') is not True:
        raise ValueError(f'{model}: both scalars must remain active in EFT1')
    mu_f, mu_s = float(eft1_endpoint['mu_f_gev']), float(eft1_endpoint['mu_s_gev'])
    if not all(math.isfinite(x) for x in (mu_f, mu_s)) or not mu_f > mu_s > 0:
        raise ValueError(f'{model}: invalid scale ordering')
    uv = load_saved_uv_endpoint(uv_endpoint_path, uv_payload, model)
    if not math.isclose(uv.mu_gev, mu_f, rel_tol=1e-12):
        raise ValueError(f'{model}: UV endpoint does not lie at muF')
    f_masses = fermion_threshold_masses(uv)
    if f_masses.minimum_mass_gev <= mu_f:
        raise ValueError(f'{model}: grouped fermion-threshold hierarchy violated')
    projected = build_intermediate_scalar_boundary(uv, matching_scale_gev=mu_f)
    initial = decode_eft1_endpoint(eft1_endpoint['initial'], projected, mu_gev=mu_f)
    compare_states(projected, initial, label=f'{model} muF tree scalar handoff')
    final = decode_eft1_endpoint(eft1_endpoint['final'], projected, mu_gev=mu_s)
    if set(intermediate_scalar_state_environment(final)) != set(names):
        raise ValueError(f'{model}: final EFT1/JSON beta inventory mismatch')
    audit = scalar_threshold_audit(final, mu_f_gev=mu_f, mu_s_gev=mu_s)
    # The scalar-only SM boundary discards both scalar fields; it does not
    # compute a Wilson coefficient nor assert matching-scale independence.
    sm = build_final_sm_boundary(final, matching_scale_gev=mu_s)
    return {
        'format': 't3-batch15-renormalisable-scalar-threshold-handoff-v1',
        'model': model,
        'status': 'CandidateScalarThresholdHandoffAudited',
        'mu_f_gev': mu_f, 'mu_s_gev': mu_s,
        'fermion_singular_masses_at_muF_gev': f_masses.masses_gev.tolist(),
        'n_eft1_couplings': len(names),
        'eft1_scalar_couplings_retained': sorted(set(names)-set(('gY','g2','g3','yu','yd','ye','lambdaH'))),
        'scalar_threshold': audit,
        'sm_boundary_at_muS': {
            'gY': sm.sm.gY, 'g2': sm.sm.g2, 'g3': sm.sm.g3,
            'lambdaH': sm.sm.lambdaH,
            'yukawa_matrices_present': ['yu','yd','ye'],
        },
        'both_scalars_integrated_out_in_renormalisable_sm_projection': True,
        'fermion_threshold_wilson_matching_performed': False,
        'eft1_seven_term_wilson_rge_performed': False,
        'scalar_threshold_c5_matching_performed': False,
        'general_quartic_matchete_rgbeta_mapping_certified': False,
        'physical_neutrino_mass_prediction': False,
        'blocking_requirements_for_physical_c5': [
            'Matchete-to-RGBeta scalar-coupling normalisation map',
            'F-threshold LLSS master Wilson tensor with certified spinor/operator factors',
            'EFT1 Wilson RGE in the same scalar and fermion component basis',
            'Complete simultaneous S1/S2 scalar matching including finite terms',
        ] + (['Physical real-component scalar mass matrix for additional quadratic term']
             if audit['mass_mixing_requires_resolution'] else []),
    }
