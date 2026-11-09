"""Read-only A/B readiness checks for one common scalar threshold.

This module validates the *actual supplied endpoints* before executing C.
It never derives or fabricates missing A/B Wilson, quartic, phase, or mass data.
Its `ready_conditional` status indicates a mathematically usable candidate,
not a physical one-loop matching certificate.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any
import math
import numpy as np

from .models import ALL_MODELS, n_real, normalize_model
from .matching import (
    check_quartic, check_rotation, extract_lepton_sectors, mass_eigenbasis,
    match_scalar_threshold, project_master, scalar_gauge_ward_residuals,
    symmetric3,
)
from .pipeline import (
    B_FORMAT, W_FORMAT, _contract, _load_pair, _load_sm,
)

PREFLIGHT_FORMAT = 't3-parallel-c-real-endpoint-preflight-v2'


def _maxabs(x: np.ndarray) -> float:
    return float(np.max(np.abs(x))) if np.size(x) else 0.0


def _check_full_wilson(c: np.ndarray, size: int, *, label: str) -> dict:
    """Assert full A tensor is finite and (ij),(ab) exchange symmetric.

    Checks use max-norm, scaled to the input's *own* Wilson magnitude.
    Individual complex entries are not conjugated in either exchange.
    """
    c = np.asarray(c)
    if c.ndim != 4 or c.shape[0] != c.shape[1] or c.shape[2:] != (size, size):
        raise ValueError(f'{label} must have shape (nW,nW,{size},{size})')
    if c.shape[0] < 6 or not np.isfinite(c).all():
        raise ValueError(f'{label} has too few Weyl entries or nonfinite entries')
    scale = max(_maxabs(c), 1e-100)
    ij = _maxabs(c - c.swapaxes(0, 1)) / scale
    ab = _maxabs(c - c.swapaxes(2, 3)) / scale
    if max(ij, ab) > 1e-8:
        raise ValueError(f'{label}: full fermion/scalar pair symmetry failed: {ij:g}, {ab:g}')
    return {'n_weyl': int(c.shape[0]), 'n_real_scalars_with_higgs': size,
            'flavour_pair_relative_residual': ij,
            'scalar_pair_relative_residual': ab}


def _verify_channels(q: np.ndarray, arrays: dict, meta: dict) -> dict:
    names = tuple(meta.get('quartic_channel_names', ()))
    if not all(isinstance(s, str) and s.strip() for s in names) or len(set(names)) != len(names):
        raise ValueError('Quartic channel names must be unique nonempty strings')
    channels = arrays.get('quartic_real_channels_ms')
    if channels is None:
        if names:
            raise ValueError('Quartic channel names supplied without component tensors')
        return {'channel_count': 0, 'status': 'unresolved_independent_invariant_decomposition'}
    a = np.asarray(channels)
    if a.ndim != 5 or a.shape[1:] != q.shape or len(a) != len(names):
        raise ValueError('Independent quartic-channel tensor shapes/names mismatch')
    if not len(names):
        raise ValueError('Quartic tensor channels require explicit names')
    for name, component in zip(names, a):
        if not np.isfinite(component).all():
            raise ValueError(f'Nonfinite quartic channel {name}')
        if not np.allclose(component, component.real, rtol=0, atol=1e-12):
            raise ValueError(f'Complex real-field quartic channel {name}')
        # Symmetry validation is a component-by-component invariant check.
        scale = max(_maxabs(component), 1e-40)
        for p in ((1,0,2,3), (0,2,1,3), (0,1,3,2)):
            if _maxabs(component-component.transpose(p)) > 1e-8*scale:
                raise ValueError(f'Quartic channel {name} not fully symmetric')
    relative = _maxabs(np.sum(a, axis=0)-q) / max(_maxabs(q), 1e-100)
    if relative > 1e-9:
        raise ValueError(f'Quartic channels omit/double-count tensors: relative={relative:g}')
    return {'channel_count': len(names), 'names': list(names),
            'sum_relative_residual': relative, 'status': 'numerically_consistent_supplied_components'}


def preflight_pair(wilson_manifest: Path, background_manifest: Path, *,
                   mode: str = 'rge_improved', mu_low: float = 173.,
                   allow_projection_residual: bool = False,
                   expected_model: str | None = None,
                   permit_synthetic: bool = False) -> dict[str, Any]:
    """Return checked, serializable readiness, with no file writes or RGE calls.

    Checks the precise tensor conventions, provenance, full Wilson and
    quartic inventory, the common eigenmass threshold, and Wilson projection.
    A/B lack of numerical endpoints is recorded as a blocker, never silently
    filled using the synthetic C fixtures or stale F-first results.
    """
    wp, bp = Path(wilson_manifest), Path(background_manifest)
    expected = normalize_model(expected_model) if expected_model is not None else None
    report: dict[str, Any] = {
        'format': PREFLIGHT_FORMAT,
        'model': expected,
        'wilson_manifest': str(wp),
        'background_manifest': str(bp),
        'requested_mode': mode,
        'mu_low_gev': float(mu_low),
        'status': 'blocked',
        'blockers': [], 'checks': {}, 'warnings': [],
        'no_synthetic_substitution': True,
        'real_physical_prediction_certified': False,
    }

    def stage(name: str, callback):
        try:
            result = callback()
            # Avoid serialising the full NPZ content into a diagnostic report.
            if name.endswith('input_manifest_and_checksum'):
                meta, arrays, json_sha, npz_sha = result
                shown = {'format': meta.get('format'),
                         'model': meta.get('model'),
                         'npz_keys': sorted(arrays),
                         'manifest_sha256': json_sha, 'npz_sha256': npz_sha}
            elif name == 'model_seed_convention_scale_contract':
                shown = {'model': result[0], 'mu_f_gev': result[1],
                         'mu_s_gev': result[2]}
            else:
                shown = result
            report['checks'][name] = {'status': 'passed', 'details': shown}
            return result
        except (ValueError, KeyError, TypeError, FileNotFoundError, OSError, RuntimeError) as error:
            msg = f'{type(error).__name__}: {error}'
            report['checks'][name] = {'status': 'blocked', 'reason': msg}
            report['blockers'].append(f'{name}: {msg}')
            return None

    wa_load = stage('A_input_manifest_and_checksum', lambda: _load_pair(wp, W_FORMAT))
    ba_load = stage('B_input_manifest_and_checksum', lambda: _load_pair(bp, B_FORMAT))
    if wa_load is None or ba_load is None:
        return report
    w, wa, w_sha, w_npz_sha = wa_load
    b, ba, b_sha, b_npz_sha = ba_load
    report['input_sha256'] = {'A_manifest': w_sha, 'A_npz': w_npz_sha,
                              'B_manifest': b_sha, 'B_npz': b_npz_sha}
    # Keep manifest content out of readiness JSON: source files may contain
    # large metadata and arbitrary private paths.
    info = stage('model_seed_convention_scale_contract', lambda: _contract(w, b))
    if info is None:
        return report
    model, mf, ms = info
    report['model'] = model
    report['mu_f_gev'], report['mu_s_gev'] = mf, ms
    if expected is not None and model != expected:
        report['blockers'].append(f'model: expected {expected}, received {model}')
        return report
    n = n_real(model)
    if (w.get('synthetic_fixture', False) or b.get('synthetic_fixture', False)) and not permit_synthetic:
        report['blockers'].append('input_kind: synthetic fixture forbidden for real integration')
        return report
    report['input_kind'] = ('synthetic_contract_probe' if w.get('synthetic_fixture', False)
                            else 'real_A_B_candidate')
    # An immutable synthetic model-source identifier can never certify a real endpoint.
    if report['input_kind'] == 'real_A_B_candidate' and any(
            'synthetic' in str(w.get(k,'')).lower() for k in
            ('lagrangian_source_id', 'matching_context_id', 'matched_seed_id')):
        report['blockers'].append('provenance: synthetic source ID on supposedly real data')
        return report
    if not (math.isfinite(mu_low) and 0 < float(mu_low) <= ms):
        report['blockers'].append(f'scales: require 0 < mu_low <= MS={ms}')
        return report

    def scalars():
        matrix = np.asarray(ba['mass_squared_gev2_ms'])
        if matrix.shape != (n,n): raise ValueError(f'B mass Hessian expected {(n,n)}')
        eigen, _ = mass_eigenbasis(matrix)
        R = check_rotation(ba['R_wilson_from_background'], n)
        labels_a = w.get('scalar_source_real_labels')
        labels_b = b.get('scalar_background_real_labels')
        if (not isinstance(labels_a, list) or len(labels_a) != n or len(set(labels_a)) != n or
            not isinstance(labels_b, list) or len(labels_b) != n or len(set(labels_b)) != n):
            raise ValueError('Complete unique A/B scalar component labels required')
        return {'n_real_heavy': n, 'positive_eigenvalues_gev2': eigen.tolist(),
                'mass_basis_orthogonality_relative': _maxabs(R.T@R-np.eye(n)),
                'scalar_labels_present': True}
    scalar_info = stage('all_heavy_scalar_mass_eigenbasis', scalars)

    def quartics():
        q = check_quartic(ba['quartic_real_full_ms'], n)
        return {'shape': list(q.shape), 'full_tensor_component_count': int(q.size),
                **_verify_channels(q, ba, b)}
    quartic_info = stage('full_scalar_quartic_and_independent_channels', quartics)

    def wilson():
        c = np.asarray(wa['C_psi2phi2_ms'])
        summary = _check_full_wilson(c, n+4, label='A full evolved Wilson')
        hh, _ = extract_lepton_sectors(c, w['lepton_indices'],
                                        w['higgs_source_indices'], w['heavy_source_indices'],
                                        number_of_real_scalars=n+4)
        c5, residual = project_master(hh, enforce_pure=not allow_projection_residual)
        if 'C5_EFT1_ms' in wa:
            separate = symmetric3(wa['C5_EFT1_ms'], 'A separately supplied C5')
            mismatch = _maxabs(separate - c5)/max(_maxabs(separate), _maxabs(c5), 1e-100)
            if mismatch > 1e-8:
                raise ValueError(f'Independent C5 projection mismatch: {mismatch:g}')
        else:
            mismatch = None
        tree = wa.get('C_psi2phi2_tree_mf')
        if tree is not None:
            _check_full_wilson(tree, n+4, label='A zero-loop Wilson at MF')
        elif mode == 'strict_one_loop':
            raise ValueError('Strict one-loop mode requires A zero-order full Wilson at MF')
        strict_endpoint=wa.get('C_psi2phi2_ms_strict_1loop')
        if mode=='strict_one_loop':
            if (strict_endpoint is None or w.get('strict_1loop_endpoint_at_ms') is not True
                    or not isinstance(w.get('strict_1loop_origin'),str)
                    or not w['strict_1loop_origin'].strip()):
                raise ValueError('Strict mode needs A one-loop-truncated FULL Wilson at MS and explicit origin/order metadata')
            _check_full_wilson(strict_endpoint,n+4,label='A strict one-loop Wilson at MS')
            strict_hh,_=extract_lepton_sectors(strict_endpoint,w['lepton_indices'],
                          w['higgs_source_indices'],w['heavy_source_indices'],
                          number_of_real_scalars=n+4)
            project_master(strict_hh,enforce_pure=not allow_projection_residual)
        if wa.get('C5_independent_finite_at_MS') is not None:
            symmetric3(wa['C5_independent_finite_at_MS'], 'independent finite C5')
            if w.get('finite_hard_ms_kind') != 'independent_finite_only':
                raise ValueError('Finite MS term is not independently certified as finite-only')
        return {**summary, 'weinberg_projection_relative_residual': residual,
                'C5_EFT1_norm_gev_inverse': float(np.linalg.norm(c5)),
                'C5_separate_projection_relative_mismatch': mismatch,
                'order_zero_full_wilson_available': tree is not None,
                'other_operator_classes_matching_certified': False}
    wilson_info = stage('full_evolved_wilson_and_C5_projection', wilson)

    def sm():
        constants, y = _load_sm({**ba, **b})
        return {'gauge_and_quartic_parameters': constants,
                'yukawa_shapes': {name: list(val.shape) for name, val in y.items()}}
    sm_info = stage('SM_couplings_for_existing_one_loop_runner', sm)

    def ward():
        g = ba.get('scalar_gauge_generators_full')
        if g is None:
            if b.get('gauge_ward_certified') is True:
                raise ValueError('B claims Ward certification without generators')
            return {'status': 'not_supplied_not_certified'}
        ans = scalar_gauge_ward_residuals(ba['quartic_real_full_ms'],
                                          ba['mass_squared_gev2_ms'], g)
        if not ans['passed']:raise ValueError(f'Gauge Ward residuals fail: {ans}')
        return ans
    if scalar_info and quartic_info:
        stage('optional_B_gauge_ward_test', ward)

    if not report['blockers']:
        def threshold():
            x = match_scalar_threshold(
                full_wilson_ms=wa['C_psi2phi2_ms'],
                full_wilson_tree_mf=wa.get('C_psi2phi2_tree_mf'),
                full_wilson_strict_ms=wa.get('C_psi2phi2_ms_strict_1loop'),
                lepton_indices=np.asarray(w['lepton_indices']),
                higgs_indices=np.asarray(w['higgs_source_indices']),
                heavy_indices=np.asarray(w['heavy_source_indices']),
                full_quartic_ms=ba['quartic_real_full_ms'],
                scalar_mass2_background=ba['mass_squared_gev2_ms'],
                R_wilson_from_background=ba['R_wilson_from_background'],
                mu_s=ms, mode=mode,
                separately_supplied_c5_ms=wa.get('C5_EFT1_ms'),
                independent_hard_c5=wa.get('C5_independent_finite_at_MS'),
                full_quartic_channels=ba.get('quartic_real_channels_ms'),
                channel_names=tuple(b.get('quartic_channel_names',())),
                enforce_pure_weinberg=not allow_projection_residual)
            return {
                'eigenstates_integrated_out_at_single_MS': x['scalar_components_removed'],
                'unordered_eigenstate_pairs': x['bubble']['pair_count'],
                'already_evolved_C5_kept_once': x['no_direct_leading_log_reinserted'],
                'bubble_source': x['bubble_Wilson_order_source'],
                'C5_EFT1_order_source':x['c5_eft1_power_counting_source'],
                'scalar_bubble_norm_gev_inverse': float(np.linalg.norm(x['c5_scalar_bubble'])),
                'C5_SMEFT_at_MS_norm_gev_inverse': float(np.linalg.norm(x['c5_smeft_at_ms'])),
                'scalar_bubble_weinberg_residual': x['scalar_bubble_projection_residual'],
                'complete_finite_matching_certified': False,
            }
        stage('single_MS_finite_bubble_match_candidate', threshold)

    if not report['blockers']:
        report['status'] = 'ready_conditional'
        for key, owner, reason in (
            ('wilson_adm_certified', 'A', 'full EFT1 dim-5 anomalous dimension'),
            ('wilson_dim5_closure_certified', 'A', 'operator closure'),
            ('signed_GammaCC_to_Weyl_certified', 'A', 'spinor normalisation'),
            ('scalar_operator_inventory_certified', 'B', 'full allowed quartic inventory'),
            ('native_matchete_rgbeta_normalization_certified', 'B', 'Matchete/RGBeta scalar tensor mapping'),
        ):
            value = (w if owner == 'A' else b).get(key, False)
            if value is not True:
                report['warnings'].append(f'{owner}: {reason} not certified')
        report['warnings'].append('C: independent additional MS one-loop diagrams, wavefunction/EOM matching uncomputed')
        report['warnings'].append('Even with certified A/B inputs, result is a conditional one-loop scalar-bubble candidate')
        if allow_projection_residual:
            report['warnings'].append('SU(2) Weinberg projection allowed nonzero residual: not a pure SMEFT dimension-five match')
    return report
