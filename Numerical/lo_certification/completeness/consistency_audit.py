"""Read-only all-16 conditional hard-portal numerical and evidence cross-audit.

This is a consistency check WITHIN one declared scalar-portal topology, not
an independent full UV calculation or a physical normalization certificate.
No existing output, upstream worker, or GitHub repository is mutated.
"""
from __future__ import annotations

from hashlib import sha256
from pathlib import Path
from typing import Any
import json
import math
import numpy as np

from .hard_matching import evaluate_hard_scalar_channel
from .inventory import ALL_MODELS


class ConsistencyError(ValueError):
    pass


def require(ok: bool, message: str) -> None:
    if not bool(ok):
        raise ConsistencyError(message)


def sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def json_dict(path: Path) -> dict:
    obj = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(obj, dict), f'Expected JSON object: {path}')
    return obj


def slug(model: str) -> str:
    require(model in ALL_MODELS, f'Unknown T3 model: {model}')
    return model.replace('-', 'm').replace('+', 'p')


def decode_matrix(obj: Any) -> np.ndarray:
    require(isinstance(obj, dict) and set(obj) == {'real', 'imag'},
            'Complex matrix needs exactly real/imag')
    a = np.asarray(obj['real'], dtype=float) + 1j * np.asarray(obj['imag'], dtype=float)
    require(a.shape == (3, 3) and np.isfinite(a).all(), 'Expected finite complex 3x3 matrix')
    return a


def relative_norm_error(actual: np.ndarray, expected: np.ndarray) -> float:
    return float(np.linalg.norm(actual - expected) /
                 max(float(np.linalg.norm(expected)), 1e-250))


def _run_kernel(z: dict, mf: float, ms: float, *, channels: bool = False,
                mass_scale: float = 1.0) -> dict:
    return evaluate_hard_scalar_channel(
        z['C_LLSS_by_fermion'], z['fermion_masses_gev'], z['lambda_HHSS'],
        z['scalar_mass_squared_gev2'] * mass_scale, mf, ms,
        C_total_muF=z['C_LLSS_total_muF'],
        lambda_channels=z['lambda_HHSS_channels'] if channels else None,
        channel_names=tuple(z['channel_names']) if channels else None,
    )


def audit_model(model: str, run_dir: Path) -> dict:
    """Recompute and independently stress-test saved conditional tensors.

    Both muF/muS variations keep Wilson coefficients, masses, and quartics
    FIXED. This verifies only algebraic partition/scale cancellation, not RGE
    transport of running couplings. The scalar mass2 scan compares the exact
    triangle with its associated leading-power expansion, not an external UV
    amplitude.
    """
    name = slug(model)
    meta_path = run_dir / 'conditional_inputs' / (name + '.json')
    calc_path = run_dir / 'conditional_hard' / (name + '.json')
    meta = json_dict(meta_path)
    stored = json_dict(calc_path)
    require(meta.get('model') == model and stored.get('model') == model,
            f'{model}: source/diagnostic model mismatch')
    require(meta.get('tensor_convention') == 'T3_REAL_LLSS_HHSS_WEINBERG_v1' and
            meta.get('scheme') == 'MSbar', f'{model}: mismatched tensor/scheme')
    require(meta.get('source_normalization_certified') is False and
            meta.get('physical_C5_certified') is False and
            stored.get('physical_C5_certified') is False,
            f'{model}: unexpected physical or source-normalization claim')
    require(stored.get('projection_passed') is True and
            stored.get('all_fermions_one_muF') is True and
            stored.get('both_scalars_one_muS') is True,
            f'{model}: missing projection or common threshold')
    require(stored.get('other_six_EFT1_RGE_sources_included') is False and
            stored.get('complete_hard_one_loop_matching_certified') is False,
            f'{model}: this script is for the incomplete portal topology only')
    provenance = meta.get('provenance')
    require(isinstance(provenance, dict), f'{model}: no source provenance')
    back_meta = run_dir / 'scalar_backgrounds' / (name + '.json')
    back_npz = run_dir / 'scalar_backgrounds' / (name + '.npz')
    require(sha(back_meta) == meta.get('source_scalar_background_metadata_sha256') and
            sha(back_npz) == meta.get('source_scalar_background_npz_sha256') and
            sha(back_npz) == provenance.get('scalar_bundle_sha256'),
            f'{model}: scalar background changed since hard calculation')
    npz = meta.get('tensor_npz')
    require(npz == name + '.npz', f'{model}: unexpected hard input NPZ name')
    npz_path = meta_path.parent / npz
    with np.load(npz_path, allow_pickle=False) as f:
        keys = {'C_LLSS_by_fermion', 'C_LLSS_total_muF', 'fermion_masses_gev',
                'lambda_HHSS', 'lambda_HHSS_channels', 'scalar_mass_squared_gev2'}
        require(keys <= set(f.files), f'{model}: missing input tensors')
        z = {k: np.asarray(f[k]) for k in keys}
    z['channel_names'] = meta.get('quartic_channel_names')
    mf, ms = float(meta['mu_f_gev']), float(meta['mu_s_gev'])
    require(math.isfinite(mf) and math.isfinite(ms) and mf >= ms > 0,
            f'{model}: invalid common thresholds')
    recomputed = _run_kernel(z, mf, ms, channels=True)
    require(recomputed['projection_passed'], f'{model}: recomputed Ward projection failed')
    matrix_keys = ('hard_muF', 'eft1_scalar_log', 'scalar_threshold_muS',
                   'dim5_sum', 'full_triangle', 'higher_dimension_remainder')
    errors = {}
    stored_matrices = stored.get('matrices', {})
    for key in matrix_keys:
        a = decode_matrix(stored_matrices[key])
        b = recomputed['matrices'][key]
        errors[key] = relative_norm_error(a, b)
        require(errors[key] < 1e-9 or np.linalg.norm(a-b) < 1e-26,
                f'{model}: saved {key} disagrees with fresh calculation')
    for key in matrix_keys:
        r = float(recomputed['residuals'][key])
        require(r <= 1e-8, f'{model}: {key} Ward projection failed ({r})')
    # Each independent native quartic channel must sum to each output.
    qsum = np.sum(z['lambda_HHSS_channels'], axis=0)
    require(np.allclose(qsum, z['lambda_HHSS'], rtol=1e-11, atol=1e-14),
            f'{model}: quartic channels do not reconstruct total')

    baseline = recomputed['matrices']['dim5_sum']
    f_change = _run_kernel(z, mf * 1.4, ms)['matrices']['dim5_sum']
    s_change = _run_kernel(z, mf, ms * 0.73)['matrices']['dim5_sum']
    scale_f = relative_norm_error(f_change, baseline)
    scale_s = relative_norm_error(s_change, baseline)
    require(max(scale_f, scale_s) < 1e-9,
            f'{model}: hard/RG/scalar matching-scale cancellation failed')

    def remainder_at_scale(s: float) -> float:
        r = recomputed if s == 1.0 else _run_kernel(z, mf, ms, mass_scale=s)
        return float(np.linalg.norm(r['matrices']['higher_dimension_remainder']))

    rem_1, rem_quarter, rem_sixteenth = (remainder_at_scale(s) for s in (1., .25, .0625))
    tiny = max(float(np.linalg.norm(baseline)), 1e-250) * 1e-12
    if rem_1 > tiny and rem_quarter > tiny and rem_sixteenth > tiny:
        exponent_first = math.log(rem_1 / rem_quarter) / math.log(4.)
        exponent_second = math.log(rem_quarter / rem_sixteenth) / math.log(4.)
        power_status = ('PASS_APPROX_DIM7_POWER' if
                        0.5 < exponent_first < 1.5 and 0.5 < exponent_second < 1.5
                        else 'FAIL_NON_DIM7_SCALING')
    else:
        exponent_first = exponent_second = None
        power_status = 'UNDETERMINED_CANCELLATION_OR_ZERO'
    require(power_status != 'FAIL_NON_DIM7_SCALING',
            f'{model}: higher-power UV remainder does not scale as scalar mass2')
    exact = recomputed['matrices']['full_triangle']
    uv_error = float(np.linalg.norm(recomputed['matrices']['higher_dimension_remainder']) /
                     max(np.linalg.norm(exact), 1e-250))
    hard = recomputed['matrices']['hard_muF']
    return {
        'model': model, 'status': 'PORTAL_NUMERIC_CHECKS_PASS' if
        power_status == 'PASS_APPROX_DIM7_POWER' else 'PORTAL_NUMERIC_PASS_POWER_UNDETERMINED',
        'hard_input_json_sha256': sha(meta_path),
        'hard_input_npz_sha256': sha(npz_path),
        'saved_diagnostic_sha256': sha(calc_path),
        'source_fermion_boundary_sha256': provenance.get('fermion_source_sha256'),
        'source_running_sha256': provenance.get('running_sha256'),
        'F_matching_scale_gev': mf, 'S_matching_scale_gev': ms,
        'hard_frobenius_per_unknown_N_gamma_gev_inverse': float(np.linalg.norm(hard)),
        'dim5_frobenius_per_unknown_N_gamma_gev_inverse': float(np.linalg.norm(baseline)),
        'hard_to_dim5_norm_ratio': float(np.linalg.norm(hard) / max(np.linalg.norm(baseline), 1e-250)),
        'relative_exact_triangle_truncation_error': uv_error,
        'max_projection_residual': max(float(x) for x in recomputed['residuals'].values()),
        'matching_scale_muF_relative_drift_fixed_inputs': scale_f,
        'matching_scale_muS_relative_drift_fixed_inputs': scale_s,
        'scalar_mass2_power_exponent_1_to_1quarter': exponent_first,
        'scalar_mass2_power_exponent_1quarter_to_1sixteenth': exponent_second,
        'UV_power_scan_status': power_status,
        'n_heavy_flavours': len(z['fermion_masses_gev']),
        'n_scalar_real_components': len(z['scalar_mass_squared_gev2']),
        'n_preserved_HHSS_source_channels': len(z['channel_names']),
        'saved_to_recomputed_relative_errors': errors,
        'physical_C5_certified': False,
    }


def crossworkstream_evidence(portal_summary: Path | None,
                             gammacc_report: Path | None,
                             previous_scalar_summary: Path | None) -> dict:
    """Treat WS2/Workstream C/old scalar output as read-only, never certificates.

    In particular, the installed Matchete GammaCC *source locator* is not a
    signed normalization calibration, and one scalar portal span is not the
    complete HHSS+UV invariant inventory.
    """
    result = {'portal': {'status': 'NOT_PROVIDED', 'models': {}},
              'gammacc': {'status': 'NOT_PROVIDED', 'signed_factor_certified': False},
              'legacy_scalar': {'status': 'NOT_PROVIDED', 'models': {}}}
    if portal_summary is not None and portal_summary.is_file():
        p = json_dict(portal_summary)
        require(p.get('format') == 't3-ws2-exact-portal-span-audit-v2' and
                p.get('all16_audit_passed') is True,
                'Workstream 2 portal-span summary format or audit status changed')
        rows = p.get('models')
        require(isinstance(rows, list) and len(rows) == 16 and
                {r.get('model') for r in rows} == set(ALL_MODELS),
                'Workstream 2 portal-span model inventory not all 16')
        result['portal'] = {
            'status': 'SOURCE_SPAN_AUDITED_NOT_GENERAL_UV_CLOSURE',
            'sha256': sha(portal_summary),
            'models': {
                r['model']: {
                    'rgbeta_directions_missing_from_native_HHSS_seed':
                        r.get('RGBeta_portal_directions_missing_from_native_seed'),
                    'extra_symmetry_allowed_HHSS_candidates': r.get('extra_HHSS_candidate_count'),
                    'S1_span_status': r.get('S1_span_status'),
                    'S2_span_status': r.get('S2_span_status'),
                    'remaining_upstream_blockers': r.get('blockers', []),
                } for r in rows},
        }
    if gammacc_report is not None and gammacc_report.is_file():
        g = json_dict(gammacc_report)
        require(g.get('format') == 't3-bc-conditional-majorana-mass-normalization-v1' and
                g.get('native_GammaCC_to_Warsaw_prefactor_certified') is False,
                'Workstream C existing conditional-normalization report changed; manual audit required')
        result['gammacc'] = {
            'status': 'GRASSMANN_AND_VEV_VERIFIED_SIGNED_MAP_OPEN',
            'sha256': sha(gammacc_report),
            'source_probe_status': g.get('installed_source_probe', {}).get('status'),
            'signed_factor_certified': False,
        }
    if previous_scalar_summary is not None and previous_scalar_summary.is_file():
        l = json_dict(previous_scalar_summary)
        require(l.get('format') == 't3-exploratory-all16-conditional-c5-benchmark-v6',
                'Earlier v6 scalar audit summary format changed')
        rows = l.get('models')
        require(isinstance(rows, list) and len(rows) == 16 and
                {r.get('model') for r in rows} == set(ALL_MODELS),
                'Earlier scalar audit model inventory changed')
        result['legacy_scalar'] = {
            'status': 'INDEPENDENT_EARLIER_SCALAR_CALCULATION_WITH_UNRESOLVED_FAILURES',
            'sha256': sha(previous_scalar_summary),
            'models': {r['model']: {'previous_status': r.get('status'),
                                  'previous_projection_failed':
                                  'Weinberg Ward projection' in r.get('error', '')}
                       for r in rows},
        }
    return result


def assemble(rows: list[dict], evidence: dict) -> dict:
    require({r.get('model') for r in rows} == set(ALL_MODELS) and len(rows) == 16,
            'Require a complete all-16 numerical consistency audit')
    combined = []
    for r in rows:
        model = r['model']
        portal = evidence['portal']['models'].get(model)
        prior = evidence['legacy_scalar']['models'].get(model)
        blockers = [
            'Signed Matchete GammaCC-to-master Weyl map not certified',
            'Native Matchete/RGBeta scalar quartic aliases not physically normalized',
            'All independent dimension-five operators and finite hard counterterms not certified',
        ]
        if portal is None:
            blockers.append('Workstream 2 portal-span evidence unavailable')
        else:
            missing = portal['rgbeta_directions_missing_from_native_HHSS_seed']
            extra = portal['extra_symmetry_allowed_HHSS_candidates']
            if missing is None or (isinstance(missing, int) and missing > 0):
                blockers.append('Native HHSS source lacks some independently counted portal directions')
            if extra is None or (isinstance(extra, int) and extra > 0):
                blockers.append('Extra symmetry-allowed HHSS candidates not source-normalized')
        if prior is not None and prior['previous_projection_failed']:
            blockers.append('Earlier independent scalar-threshold Weinberg projection failed')
        if r['status'] != 'PORTAL_NUMERIC_CHECKS_PASS':
            blockers.append('Numerical UV power scan not positively verified')
        combined.append({**r, 'portal_span': portal,
                         'earlier_scalar_comparison': prior,
                         'physical_blockers': blockers, 'physical_C5_certified': False})
    numeric_ok = sum(r['status'] == 'PORTAL_NUMERIC_CHECKS_PASS' for r in rows)
    return {
        'format': 't3-lo-conditional-portal-crossworkstream-audit-v1',
        'model_count': 16, 'portal_numeric_consistency_passed': numeric_ok,
        'portal_numeric_incomplete_or_failed': 16 - numeric_ok,
        'physical_full_one_loop_certified': False,
        'physical_C5_certified': False,
        'source_evidence': {k: {kk: vv for kk, vv in v.items() if kk != 'models'}
                            for k, v in evidence.items()},
        'previous_scalar_projection_failures': [r['model'] for r in combined
            if (r['earlier_scalar_comparison'] or {}).get('previous_projection_failed')],
        'models_with_missing_native_portal_directions': [r['model'] for r in combined
            if r['portal_span'] is not None and
               isinstance(r['portal_span']['rgbeta_directions_missing_from_native_HHSS_seed'], int) and
               r['portal_span']['rgbeta_directions_missing_from_native_HHSS_seed'] > 0],
        'models': combined,
        'eft_sequence': 'UV -> all F at ONE muF -> EFT1(S1,S2,SM) -> both S at ONE muS -> SMEFT',
        'numerical_scope': 'source-matched scalar-portal triangle, fixed-input matching-scale test',
        'UV_power_test_scope': 'same analytic UV triangle and its own dim5 expansion, not independent UV derivation',
        'normalization': 'per unknown signed N_GammaCC; physical sign not assumed',
    }
