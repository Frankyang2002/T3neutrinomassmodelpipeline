"""Diagnose T3-B direct-Weinberg matrix alignment with heavy-resolved logarithms.

Runs the existing renormalisable trajectory and production hard/direct evaluator.
Compares the *shape* of two independent leading-log kernels, allowing a best-fit
complex overall normalisation for diagnosis only. Never reports corrected C5.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from types import SimpleNamespace
import numpy as np

def matrix_from_json(p):
    a = np.asarray(p['real'], float) + 1j*np.asarray(p['imag'], float)
    if a.shape != (3,3) or not np.all(np.isfinite(a)):
        raise ValueError('Expected finite 3x3 complex matrix')
    return a


def matrix_to_json(a):
    a = np.asarray(a, complex)
    return {'real': a.real.tolist(), 'imag': a.imag.tolist(),
            'frobenius_norm': float(np.linalg.norm(a))}


def flavor_kernel_by_heavy(y1, y2, masses):
    """Symmetrised C12_r, matching the tree full-flavour convention."""
    a, b = np.asarray(y1, complex), np.asarray(y2, complex)
    m = np.asarray(masses, float)
    if a.shape != b.shape or a.ndim != 2 or a.shape[0] != 3 or a.shape[1] != len(m):
        raise ValueError('Yukawa matrices must be matching 3 x n_heavy arrays')
    if np.any(m <= 0) or not np.all(np.isfinite(m)):
        raise ValueError('Heavy masses must be positive and finite')
    out = []
    for r, mr in enumerate(m):
        u, v = np.conj(a[:, r]), np.conj(b[:, r])
        out.append((np.outer(u, v) + np.outer(v, u)) / (2.0 * mr))
    return np.asarray(out)


def fit_shape(target, candidate):
    t, c = np.asarray(target, complex), np.asarray(candidate, complex)
    if t.shape != (3, 3) or c.shape != (3, 3):
        raise ValueError('Expected full 3x3 flavour matrices')
    if np.linalg.norm(t) <= 1e-30 or np.linalg.norm(c) <= 1e-30:
        raise ValueError('Cannot align vanishing matrices')
    alpha = np.vdot(c, t) / np.vdot(c, c)
    residual = np.linalg.norm(t - alpha*c) / np.linalg.norm(t)
    return complex(alpha), float(residual)


def analyze(direct, y1, y2, masses, ms, mf_common, *, bridge_ll=None):
    masses = np.asarray(masses, float)
    ms, mf_common = float(ms), float(mf_common)
    if not (0 < ms < min(masses) and mf_common > ms):
        raise ValueError('Threshold ordering invalid')
    tensors = flavor_kernel_by_heavy(y1, y2, masses)
    resolved_logs = np.log(ms / masses)
    common_log = np.log(ms / mf_common)
    common = np.sum(tensors, axis=0) * common_log
    resolved = np.einsum('r,rij->ij', resolved_logs, tensors)
    f_com, r_com = fit_shape(direct, common)
    f_res, r_res = fit_shape(direct, resolved)
    prod_norm = float(np.linalg.norm(direct))
    result = {
        'status': 'HeavyLogAlignmentDiagnosticComplete',
        'scope': 'matrix-shape comparison only; best-fit factors are NOT Wilson normalisation fixes',
        'heavy_masses_gev': masses.tolist(),
        'ms_gev': ms, 'mf_common_gev': mf_common,
        'individual_log_ms_over_mr': resolved_logs.tolist(),
        'common_log_ms_over_mf': float(common_log),
        'common_log': {'factor': [f_com.real, f_com.imag], 'alignment_residual': r_com},
        'heavy_resolved_log': {'factor': [f_res.real, f_res.imag], 'alignment_residual': r_res},
        'alignment_improvement': float(r_com-r_res),
        'heavy_logs_explain_mismatch': bool(r_res < 1e-5),
        'production_direct_c5_norm': prod_norm,
        'production_direct_c5': matrix_to_json(direct),
        'heavy_resolved_kernel': matrix_to_json(resolved),
        'heavy_common_kernel': matrix_to_json(common),
        'neutrino_mass_prediction_available': False,
    }
    if bridge_ll is not None:
        f_b, r_b = fit_shape(direct, bridge_ll)
        result['existing_bridge_common_ll'] = {
            'factor': [f_b.real, f_b.imag], 'alignment_residual': r_b,
        }
    return result


def run_production(config, uv_path, eft1_path, final_path):
    from Numerical.core.RGBetaEvaluator import load_rgbeta_payload
    from Numerical.fitting.ScanCLI import build_uv_state_from_config
    from Numerical.running.T3Trajectory import run_renormalisable_t3_trajectory
    from Numerical.diagnostics.FinalC5ContributionDiagnostics import evaluate_final_c5_contributions
    from Numerical.running.FinalC5TrajectoryAdapter import extract_final_c5_inputs
    cfg = json.loads(Path(config).read_text(encoding='utf-8'))
    raw = dict(cfg)
    raw.setdefault('scan', {'bindings': {}})
    state = build_uv_state_from_config(SimpleNamespace(raw=raw), {}).validated()
    rep = state.representation
    if (rep.d_s1, rep.d_s2, rep.d_f, rep.alpha, rep.shared_scalar) != (2,2,1,-1,False):
        raise ValueError('Requires ordinary T3-B alpha=-1')
    scales, numerical = cfg['scales'], cfg.get('numerical', {})
    trajectory = run_renormalisable_t3_trajectory(
        state, load_rgbeta_payload(uv_path), load_rgbeta_payload(eft1_path),
        mu_fermion_threshold_gev=float(scales['mu_fermion_threshold_gev']),
        mu_scalar_threshold_gev=float(scales['mu_scalar_threshold_gev']),
        rtol=float(numerical.get('rtol', 1e-8)),
        atol=float(numerical.get('atol', 1e-11)))
    direct = evaluate_final_c5_contributions(final_path, trajectory).direct_running
    inputs = extract_final_c5_inputs(trajectory)
    return direct, inputs, scales


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('config','uv-rgbeta','eft1-rgbeta','final-c5','bridge','output'):
        p.add_argument('--'+key, type=Path, required=True)
    args = p.parse_args()
    direct, inputs, scales = run_production(args.config, args.uv_rgbeta, args.eft1_rgbeta, args.final_c5)
    bridge = json.loads(args.bridge.read_text(encoding='utf-8'))
    if bridge.get('status') != 'OffDiagonalFullFlavorDiagnosticComplete':
        raise ValueError('Bridge diagnostic has unexpected status')
    result = analyze(direct, inputs.y1, inputs.y2, inputs.heavy_masses_gev,
                     inputs.grouped_scalar_matching_scale_gev,
                     scales['mu_fermion_threshold_gev'],
                     bridge_ll=matrix_from_json(bridge['c5_at_ms_frozen_mf_leading_log']))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n', encoding='utf-8')
    print(json.dumps({k:result[k] for k in ('status','common_log','heavy_resolved_log',
         'existing_bridge_common_ll','alignment_improvement','heavy_logs_explain_mismatch',
         'neutrino_mass_prediction_available')},indent=2))

if __name__=='__main__':
    main()
