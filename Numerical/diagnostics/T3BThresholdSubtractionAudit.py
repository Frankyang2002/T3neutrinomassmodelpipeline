"""T3-B benchmark: verify full-flavour direct-C5 normalization and no-double-counting.

Re-runs the existing renormalisable benchmark and production hard/direct C5
breakdown, then compares the direct matrix with the independent bridge's
frozen-MF leading log. Only if they align does it report a diagnostic
higher-order off-diagonal-only correction; it NEVER replaces production C5.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from Numerical.diagnostics.T3BFlavorMixingBridge import matrix_from_json, matrix_to_json


def calibrated_subtraction(control: dict, bridge: dict, hard: np.ndarray,
                           direct: np.ndarray, *, relative_tolerance=1e-5) -> dict:
    if control.get('status') != 'FixedOrderControlComplete':
        raise ValueError('Expected completed fixed-order control')
    if bridge.get('status') != 'OffDiagonalFullFlavorDiagnosticComplete':
        raise ValueError('Expected completed full-flavour off-diagonal bridge')
    if control.get('model') != 'T3-B alpha=-1' or bridge.get('model') != 'T3-B alpha=-1':
        raise ValueError('Only T3-B alpha=-1 supported')
    for scale in ('mu_f_gev', 'mu_s_gev'):
        if not np.isclose(float(control[scale]), float(bridge[scale]), rtol=1e-10):
            raise ValueError(f'Mismatched benchmark scale: {scale}')
    fixed = matrix_from_json(control['fixed_order_c5_at_ms'])
    bridge_fixed = matrix_from_json(bridge['fixed_order_c5_at_ms'])
    hard = np.asarray(hard, complex)
    direct = np.asarray(direct, complex)
    if hard.shape != (3, 3) or direct.shape != (3, 3):
        raise ValueError('Require 3x3 physical Majorana hard and direct matrices')
    if not np.all(np.isfinite(hard)) or not np.all(np.isfinite(direct)):
        raise ValueError('Non-finite production matrices')
    norm_fixed = max(np.linalg.norm(fixed), 1e-30)
    decomposition_residual = float(np.linalg.norm(fixed - hard - direct) / norm_fixed)
    bridge_residual = float(np.linalg.norm(fixed - bridge_fixed) / norm_fixed)
    if decomposition_residual > relative_tolerance or bridge_residual > relative_tolerance:
        raise ValueError('Production fixed-order C5 disagrees with hard+direct or bridge control')
    ll = matrix_from_json(bridge['c5_at_ms_frozen_mf_leading_log'])
    evolved = matrix_from_json(bridge['c5_at_ms_coupled_offdiagonal_only'])
    forward = matrix_from_json(bridge['c5_at_ms_forward_only_running_lambda'])
    ll_norm = np.linalg.norm(ll)
    direct_norm = np.linalg.norm(direct)
    if ll_norm < 1e-30 or direct_norm < 1e-30:
        raise ValueError('Cannot calibrate vanishing leading-log/direct matrices')
    # Match the complete complex 3x3 matrix, not just its Frobenius norm.
    factor = np.vdot(ll, direct) / np.vdot(ll, ll)
    alignment_residual = float(np.linalg.norm(direct - factor * ll) / direct_norm)
    if alignment_residual > relative_tolerance:
        raise ValueError(f'Production direct C5 is not proportional to bridge LL: residual={alignment_residual:.6g}')
    # Subtract exactly what the production result already contains.
    extra_forward = factor * (forward - ll)
    extra_reverse = factor * (evolved - forward)
    extra_total = extra_forward + extra_reverse
    proxy = fixed + extra_total
    return {
        'status': 'OffDiagonalThresholdSubtractionDiagnosticComplete',
        'model': 'T3-B alpha=-1',
        'scope': 'higher-order OFF-DIAGONAL-only corrections, not complete physical C5',
        'direct_vs_leading_log_complex_factor': [float(factor.real), float(factor.imag)],
        'production_hard_plus_direct_residual': decomposition_residual,
        'bridge_fixed_order_residual': bridge_residual,
        'direct_vs_leading_log_alignment_residual': alignment_residual,
        'production_fixed_c5': matrix_to_json(fixed),
        'production_hard_c5': matrix_to_json(hard),
        'production_direct_c5': matrix_to_json(direct),
        'extra_due_to_running_lambdaT3_offdiagonal_only': matrix_to_json(extra_forward),
        'extra_due_to_reverse_mixing_offdiagonal_only': matrix_to_json(extra_reverse),
        'additional_offdiagonal_only_c5': matrix_to_json(extra_total),
        'comparison_proxy_c5_not_physical_prediction': matrix_to_json(proxy),
        'relative_additional_offdiagonal_to_fixed_c5': float(np.linalg.norm(extra_total)/norm_fixed),
        'neutrino_mass_prediction_available': False,
        'full_flavor_anomalous_dimensions_verified': False,
        'warning': 'Comparison proxy omits diagonal flavour evolution and operator closure; do not use as a corrected neutrino-mass prediction.'
    }


def alignment_failure_diagnostic(control: dict, bridge: dict, hard: np.ndarray,
                                 direct: np.ndarray) -> dict:
    """Report a failed alignment without claiming a physical correction."""
    fixed = matrix_from_json(control['fixed_order_c5_at_ms'])
    bridge_fixed = matrix_from_json(bridge['fixed_order_c5_at_ms'])
    ll = matrix_from_json(bridge['c5_at_ms_frozen_mf_leading_log'])
    hard, direct = np.asarray(hard, complex), np.asarray(direct, complex)
    direct_norm = float(np.linalg.norm(direct))
    ll_norm = float(np.linalg.norm(ll))
    if direct_norm <= 1e-30 or ll_norm <= 1e-30:
        raise ValueError('Cannot diagnose vanishing direct/leading-log matrix')
    factor = np.vdot(ll, direct) / np.vdot(ll, ll)
    discrepancy = direct - factor * ll
    norm_fixed = max(float(np.linalg.norm(fixed)), 1e-30)
    return {
        'status': 'AlignmentMismatch',
        'model': 'T3-B alpha=-1',
        'scope': 'diagnostic only; no corrected C5 or mass prediction',
        'production_hard_plus_direct_residual': float(np.linalg.norm(fixed - hard - direct) / norm_fixed),
        'bridge_fixed_order_residual': float(np.linalg.norm(fixed - bridge_fixed) / norm_fixed),
        'direct_vs_leading_log_complex_factor': [float(factor.real), float(factor.imag)],
        'direct_vs_leading_log_alignment_residual': float(np.linalg.norm(discrepancy) / direct_norm),
        'direct_c5_frobenius_norm': direct_norm,
        'bridge_ll_frobenius_norm': ll_norm,
        'alignment_discrepancy_matrix': matrix_to_json(discrepancy),
        'heavy_masses_gev': bridge.get('heavy_masses_gev'),
        'possible_explanation_not_verified': 'Production hierarchical heavy-mass logarithms differ from a common MF logarithm; check source-level matching conventions.',
        'subtraction_accepted': False,
        'physical_c5_prediction_available': False,
        'neutrino_mass_prediction_available': False,
    }


def rerun_production_breakdown(config: Path, uv_path: Path, eft1_path: Path, final_path: Path):
    from Numerical.core.RGBetaEvaluator import load_rgbeta_payload
    from Numerical.fitting.ScanCLI import build_uv_state_from_config
    from Numerical.running.T3Trajectory import run_renormalisable_t3_trajectory
    from Numerical.diagnostics.FinalC5ContributionDiagnostics import evaluate_final_c5_contributions
    cfg = json.loads(config.read_text(encoding='utf-8'))
    raw = dict(cfg)
    raw.setdefault('scan', {'bindings': {}})
    state = build_uv_state_from_config(SimpleNamespace(raw=raw), {}).validated()
    rep = state.representation
    if (rep.d_s1, rep.d_s2, rep.d_f, rep.alpha, rep.shared_scalar) != (2,2,1,-1,False):
        raise ValueError('Requires T3-B alpha=-1')
    scales = cfg['scales']
    numerical = cfg.get('numerical', {})
    traj = run_renormalisable_t3_trajectory(
        state, load_rgbeta_payload(uv_path), load_rgbeta_payload(eft1_path),
        mu_fermion_threshold_gev=float(scales['mu_fermion_threshold_gev']),
        mu_scalar_threshold_gev=float(scales['mu_scalar_threshold_gev']),
        rtol=float(numerical.get('rtol', 1e-8)),
        atol=float(numerical.get('atol', 1e-11)))
    breakdown = evaluate_final_c5_contributions(final_path, traj)
    return breakdown.hard, breakdown.direct_running


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('config', 'uv-rgbeta', 'eft1-rgbeta', 'final-c5', 'control', 'bridge', 'output'):
        parser.add_argument('--' + name, required=True, type=Path)
    args = parser.parse_args()
    control = json.loads(args.control.read_text(encoding='utf-8'))
    bridge = json.loads(args.bridge.read_text(encoding='utf-8'))
    hard, direct = rerun_production_breakdown(args.config, args.uv_rgbeta,
                                              args.eft1_rgbeta, args.final_c5)
    try:
        result = calibrated_subtraction(control, bridge, hard, direct)
    except ValueError as exc:
        if 'not proportional to bridge LL' not in str(exc):
            raise
        result = alignment_failure_diagnostic(control, bridge, hard, direct)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    summary_keys = ('status', 'production_hard_plus_direct_residual',
                    'direct_vs_leading_log_alignment_residual',
                    'direct_vs_leading_log_complex_factor',
                    'neutrino_mass_prediction_available')
    print(json.dumps({**{key: result[key] for key in summary_keys},
                      **({'relative_additional_offdiagonal_to_fixed_c5':
                          result['relative_additional_offdiagonal_to_fixed_c5']}
                         if 'relative_additional_offdiagonal_to_fixed_c5' in result
                         else {})}, indent=2))


if __name__ == '__main__':
    main()
