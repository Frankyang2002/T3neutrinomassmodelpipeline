"""A(-2) reduced-scalar UV mass-boundary feasibility pilot, read-only inputs.

Replays actual zero-additional-quartic candidate one-loop UV/EFT1 beta equations.
At FIXED historical (or explicitly selected frozen) Yukawas, reconstructs the
2x2 affine map from UV mS1Sq/mS2Sq to EFT1 mS1Sq/mS2Sq at muS. Tests whether
positive UV inputs can give 0.5-5 TeV scalar masses without relaxing the
fermion-first hierarchy or other existing screens.

Does not refit Yukawas, change lambdaT3, match a new C5 or alter old results.
Negative UV roots are mathematical extrapolations, NEVER accepted inputs.

Example:
  python -u -m scripts.ProbeT3ReducedUVScalarBoundaries --run-id a_m2_historical_v1
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'output/studies/reduced_scalar_uv_boundary_pilot'
MODEL = 'A_-2'
CASE = 'small_uv_running_5e-9'
LOW_TEV = .5
HIGH_TEV = 5.
UV_GEV = 1e7
MU_F_GEV = 1e5
MU_S_GEV = 1e3
AFFINE_TOLERANCE_TEV2 = .05
FORMAT = 't3-reduced-scalar-a2-uv-boundary-feasibility-v1'


def require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError('A(-2) UV scalar boundary pilot: ' + message)


def array2(x, name: str) -> np.ndarray:
    arr = np.asarray(x, dtype=float)
    require(arr.shape == (2,) and np.all(np.isfinite(arr)),
            f'{name} must be a finite length-two vector')
    return arr


def build_affine_map(
    x0_tev2: np.ndarray,
    y0_tev2: np.ndarray,
    y_shift1_tev2: np.ndarray,
    y_shift2_tev2: np.ndarray,
    step_tev2: float,
) -> tuple[np.ndarray, np.ndarray]:
    """y = A x + b, where x and y are scalar mass squares in TeV^2."""
    x0 = array2(x0_tev2, 'x0')
    y0 = array2(y0_tev2, 'y0')
    y1 = array2(y_shift1_tev2, 'shift1')
    y2 = array2(y_shift2_tev2, 'shift2')
    require(math.isfinite(step_tev2) and step_tev2 > 0,
            'finite positive probe step required')
    a = np.column_stack(((y1 - y0) / step_tev2,
                         (y2 - y0) / step_tev2))
    b = y0 - a @ x0
    require(np.isfinite(a).all() and np.isfinite(b).all(),
            'nonfinite affine scalar mass map')
    return a, b


def check_affinity(a, b, samples, tolerance=AFFINE_TOLERANCE_TEV2):
    """Fail closed if actual RG flows disagree with the proposed affine map."""
    a = np.asarray(a, dtype=float)
    b = array2(b, 'offset')
    require(a.shape == (2, 2) and np.isfinite(a).all(), 'invalid map matrix')
    worst = 0.
    for x, y in samples:
        xx = array2(x, 'verification input')
        yy = array2(y, 'verification output')
        worst = max(worst, float(np.max(np.abs(a @ xx + b - yy))))
    require(math.isfinite(tolerance) and tolerance > 0 and worst <= tolerance,
            f'affine map not numerically validated; max residual={worst:.6g} TeV^2')
    return worst


def signed_target_solution(a, b, target_tev):
    require(math.isfinite(target_tev) and LOW_TEV <= target_tev <= HIGH_TEV,
            'target mass outside retained 0.5-5 TeV window')
    a = np.asarray(a, dtype=float)
    b = array2(b, 'offset')
    require(a.shape == (2, 2) and np.isfinite(a).all(), 'bad map matrix')
    condition = float(np.linalg.cond(a))
    if not math.isfinite(condition) or condition > 1e10:
        return {'status': 'IllConditionedMap', 'condition_number': condition}
    target_sq = np.array([target_tev**2, target_tev**2])
    x = np.linalg.solve(a, target_sq - b)
    return {'status': 'SolvedAlgebraically',
            'uv_mass_sq_tev2': x.tolist(),
            'positive_uv_mass_squared': bool(np.all(x >= 0)),
            'condition_number': condition,
            'note': 'Signed algebraic extrapolation only; negative inputs are not physical benchmark candidates.'}


def find_mass_window(a, b, min_uv_mass_tev=0., max_uv_mass_tev=None):
    """Linear-program feasibility of 0.5..5 TeV EFT masses with UV mass bounds.

    Affine map and linear constraints apply to mass SQUARED, in TeV^2.
    This does not by itself verify physical threshold/hierarchy/perturbativity.
    """
    a = np.asarray(a, dtype=float)
    b = array2(b, 'offset')
    require(a.shape == (2, 2) and np.isfinite(a).all(), 'invalid map matrix')
    require(math.isfinite(min_uv_mass_tev) and min_uv_mass_tev >= 0,
            'UV lower bound must be nonnegative')
    if max_uv_mass_tev is not None:
        require(math.isfinite(max_uv_mass_tev) and max_uv_mass_tev >= min_uv_mass_tev,
                'UV upper bound smaller than lower bound')
    lo = np.full(2, LOW_TEV**2)
    hi = np.full(2, HIGH_TEV**2)
    bounds = [(min_uv_mass_tev**2,
               None if max_uv_mass_tev is None else max_uv_mass_tev**2)] * 2
    # A x + b <= hi  AND -A x - b <= -lo.
    result = linprog(np.ones(2), A_ub=np.vstack((a, -a)),
                     b_ub=np.concatenate((hi-b, b-lo)),
                     bounds=bounds, method='highs')
    record = {'feasible': bool(result.success),
              'solver_status': int(result.status),
              'solver_message': str(result.message),
              'uv_min_mass_tev': float(min_uv_mass_tev),
              'uv_max_mass_tev': None if max_uv_mass_tev is None else float(max_uv_mass_tev)}
    if result.success:
        x = np.asarray(result.x, dtype=float)
        y = a @ x + b
        require(np.all(x >= -1e-10) and np.all(y >= lo-1e-4) and np.all(y <= hi+1e-4),
                'solver result violates the scalar window')
        record.update({'candidate_uv_mass_sq_tev2': x.tolist(),
                       'predicted_eft_mass_sq_tev2': y.tolist(),
                       'predicted_eft_mass_tev': np.sqrt(np.maximum(0., y)).tolist()})
    return record


def yukawa_input(hist, source: str, *, frozen_fit=None):
    """Select a source explicitly; frozen Yukawas are NOT identified as UV-fitted."""
    from scripts import RunT3GeneralScalar48Campaign as bulk
    if source == 'historical':
        return np.asarray(hist['seed'], dtype=float), {
            'source': 'historical_running_uv_Yukawas',
            'fitted_config_path': str(hist['cfg_path'].relative_to(ROOT)).replace('\\', '/'),
            'fitted_config_sha256': bulk.sha(hist['cfg_path']),
            'historical_fit_sha256': bulk.sha(hist['fit_path'])}, []
    require(source == 'frozen-newfit', 'unknown Yukawa source')
    path = frozen_fit or (ROOT / 'output/studies/reduced_scalar_sobol/'
                          'reduced_frozen_overnight_v2/small_uv_frozen_5e-9/A_m2/fit.json')
    path = Path(path).resolve()
    require(path.is_file(), 'missing frozen fit record: ' + str(path))
    data = bulk.read(path)
    require(data.get('model') == MODEL and data.get('case') == 'small_uv_frozen_5e-9' and
            data.get('status') == 'Accepted', 'not an accepted A(-2) frozen-low fit')
    best = (data.get('best_mass_safe') or {}).get('parameters') or {}
    require(all(k in best for k in bulk.NAMES), 'missing fitted Yukawa parameters')
    x = np.array([best[k] for k in bulk.NAMES], dtype=float)
    require(np.isfinite(x).all() and x.shape == (18,), 'invalid fitted Yukawas')
    return x, {'source': 'frozen_low_optimised_Yukawas_used_as_hypothetical_UV_input',
               'fit_path': str(path), 'fit_sha256': bulk.sha(path)}, [path]


def physical_rge_oracle(*, source, run_dir):
    """Prepare native candidate one-loop equations once and return RG(x)."""
    from scripts import RunT3GeneralScalar48Campaign as bulk
    from scripts import RunT3ReducedScalarSobolComparison as reduced
    from scripts import CoupledUVRunning as s31

    reg = bulk.read(bulk.REGISTRY)
    require(len(reg.get('models', {})) == 16, 'unexpected quartic registry')
    entry, new, reused = bulk.inventory(reg, MODEL)
    require(not reused and len(new) == 5 and entry['quadratics_not_implemented'] == 1,
            'unexpected A(-2) scalar inventory')
    hist = bulk.historical(CASE, MODEL, entry)
    require(hist['uv_on'] and hist['lambda_fixed'] == 5e-9,
            'historical UV-running lambdaT3 configuration differs')
    y, metadata, extra_paths = yukawa_input(hist, source)
    cfg = bulk.update_yukawas(hist, y)
    x_original = np.array([cfg['base_state']['ordinary'][k]
                           for k in ('mS1Sq', 'mS2Sq')], dtype=float) / 1e6
    # The RGBeta systems are source-verified *candidate* equations; no
    # Matchete coupling normalization or new scalar threshold matching claimed.
    _, layout = s31.prepare_state(cfg, entry, reused, {})
    uv_beta, uv_beta_path = bulk.uv_betas(MODEL, entry, new, layout)
    eft_payload, eft_beta_path = reduced.verify_existing_beta(
        MODEL, entry, new, reused, run_dir,
        wolframscript='wolframscript', timeout=3600)
    # Fail fast on missing beta keys with a real local state.
    from Numerical.core.StateVector import unpack_uv_state
    initial_vec, lay = s31.prepare_state(cfg, entry, reused, {})
    from Numerical.running.FermionThresholdBoundary import build_intermediate_scalar_boundary
    from scripts.RunT3GeneralScalar48Campaign import validate_eft_betas
    temporary_uv = unpack_uv_state(initial_vec, lay, mu_gev=MU_F_GEV)
    validate_eft_betas(eft_payload,
        build_intermediate_scalar_boundary(temporary_uv, matching_scale_gev=MU_F_GEV),
        MODEL, new, reused)

    def forward(x_tev2, label):
        x = array2(x_tev2, 'UV scalar boundary')
        require(np.all(x >= 0), 'physical RG probe requires nonnegative UV mS^2')
        candidate = bulk.update_yukawas(hist, y)
        ordinary = candidate['base_state']['ordinary']
        ordinary['mS1Sq'] = float(x[0] * 1e6)
        ordinary['mS2Sq'] = float(x[1] * 1e6)
        vec, lay = s31.prepare_state(candidate, entry, reused, {})
        uv, q_f, uv_nfev = bulk.run_uv(vec, lay, uv_beta, new, (0j,) * len(new))
        eft, physical_masses, q_s, eft_nfev = bulk.run_eft(
            uv, q_f, eft_payload, MODEL, new, reused)
        y_sq = np.array([eft.mS1Sq, eft.mS2Sq], dtype=float) / 1e6
        require(np.all(np.isfinite(y_sq)), 'nonfinite scalar endpoint mass squared')
        chk = bulk.testscreen(uv, eft, physical_masses, q_f, q_s)
        return {
            'label': label,
            'uv_scalar_msq_tev2': x.tolist(),
            'fermion_threshold_scalar_msq_tev2':
                [float(uv.mS1Sq)/1e6, float(uv.mS2Sq)/1e6],
            'scalar_threshold_msq_tev2': y_sq.tolist(),
            'scalar_threshold_mass_tev': (np.asarray(physical_masses, dtype=float)/1e3).tolist(),
            'screen': chk,
            'uv_nfev': int(uv_nfev), 'eft1_nfev': int(eft_nfev),
        }

    paths = [hist['cfg_path'], hist['fit_path'], bulk.REGISTRY,
             uv_beta_path, eft_beta_path, ROOT/'scripts/RunT3ReducedScalarSobolComparison.py',
             ROOT/'scripts/RunT3GeneralScalar48Campaign.py', *extra_paths]
    files = {str(p.relative_to(ROOT)).replace('\\', '/') if p.is_relative_to(ROOT) else str(p):
             bulk.sha(p) for p in paths}
    return forward, x_original, metadata, files


def run(*, run_id: str, source: str, target_tev: float, probe_shift_tev: float):
    require(bool(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}', run_id)),
            'invalid run ID')
    require(math.isfinite(target_tev) and LOW_TEV <= target_tev <= HIGH_TEV,
            'target mass must be 0.5-5 TeV')
    require(math.isfinite(probe_shift_tev) and .5 <= probe_shift_tev <= 10,
            'probe shift must be 0.5-10 TeV')
    base = OUT / run_id
    require(not base.exists(), 'output already exists; choose new run ID (no overwrite)')
    base.mkdir(parents=True, exist_ok=False)
    forward, x0, source_data, files = physical_rge_oracle(source=source, run_dir=base)
    # Distinct one-loop RG solutions with ALL additional scalar couplings zero
    # at the UV boundary. No mass-squared boundary is silently changed.
    shift = probe_shift_tev**2
    ex = np.array([shift, 0.])
    ey = np.array([0., shift])
    probe_def = [('original', x0), ('shift_S1', x0+ex),
                 ('shift_S2', x0+ey), ('zero_mass_squared', np.zeros(2)),
                 ('combined_shift', x0+ex+ey)]
    probes = []
    for label, x in probe_def:
        print(f'RUN {MODEL} {label}: UV mS^2={x.tolist()} TeV^2', flush=True)
        probes.append(forward(x, label))
    a, b = build_affine_map(x0, probes[0]['scalar_threshold_msq_tev2'],
                             probes[1]['scalar_threshold_msq_tev2'],
                             probes[2]['scalar_threshold_msq_tev2'], shift)
    check_resid = check_affinity(a, b, [
        (np.zeros(2), probes[3]['scalar_threshold_msq_tev2']),
        (x0+ex+ey, probes[4]['scalar_threshold_msq_tev2'])])
    exact = signed_target_solution(a, b, target_tev)
    unrestricted = find_mass_window(a, b, 0., None)
    light = find_mass_window(a, b, LOW_TEV, HIGH_TEV)
    verified = {}
    for label, rec in [('positive_uv', unrestricted), ('light_positive_uv', light)]:
        if rec['feasible']:
            print('VERIFY', label, flush=True)
            point = forward(np.maximum(0., rec['candidate_uv_mass_sq_tev2']), f'verify_{label}')
            sq = array2(point['scalar_threshold_msq_tev2'], label + ' observed output')
            predicted = array2(rec['predicted_eft_mass_sq_tev2'], label + ' predicted')
            resid = float(np.max(np.abs(sq - predicted)))
            pass_window = bool(np.all(sq >= LOW_TEV**2 - 1e-4) and
                               np.all(sq <= HIGH_TEV**2 + 1e-4))
            verified[label] = {'rg_result': point, 'max_affine_residual_tev2': resid,
                'scalar_window_verified': bool(pass_window and resid <= AFFINE_TOLERANCE_TEV2),
                'all_original_screens_pass': bool(point['screen']['screen'])}
    any_pass = any(x['scalar_window_verified'] and x['all_original_screens_pass']
                   for x in verified.values())
    report = {
        'format': FORMAT, 'model': MODEL, 'scenario': CASE,
        'status': 'PositiveUVBoundaryAndScreensFound' if any_pass else
                  'NoViablePositiveUVBoundaryAtTheseFixedYukawas',
        'fit_optimisation_performed': False,
        'yukawa_boundary': source_data,
        'lambdaT3_input': 5e-9,
        'scales_gev': {'uv': UV_GEV, 'fermion_threshold': MU_F_GEV,
                       'scalar_threshold': MU_S_GEV},
        'fixed_scalar_window_tev': [LOW_TEV, HIGH_TEV],
        'original_uv_scalar_msq_tev2': x0.tolist(),
        'original_baseline_fermion_rhoF': probes[0]['screen'].get('rhoF'),
        'probe_shift_mass_sq_tev2': shift,
        'affine_map': {'A': a.tolist(), 'offset_b_tev2': b.tolist(),
                       'y_tev2_equals': 'A @ x_tev2 + b',
                       'validation_max_abs_error_tev2': check_resid},
        'exact_equal_scalar_mass_target_tev': target_tev,
        'signed_exact_target': exact,
        'nonnegative_uv_mass_squared_feasibility': unrestricted,
        'light_positive_uv_mass_feasibility': light,
        'forward_verifications': verified,
        'forward_probes': probes,
        'source_sha256': files,
        'limits': [
          'Candidate RGBeta one-loop equations, not a physically certified Matchete/RGBeta basis.',
          'Fixed Yukawa, MF and quartic inputs; only two UV scalar squared-masses varied.',
          'Affineness applies to this mass-independent one-loop RG truncation, not arbitrary EFT matching.',
          'The endpoint mass test is distinct from complete trajectory vacuum-stability checks.',
          'Signed negative squared-mass extrapolations are not valid positive-UV benchmarks.',
          'No neutrino refitting, new C5 matching, or proof of impossibility after varying other UV inputs.',
        ],
    }
    out = base/'result.json'
    with out.open('x', encoding='utf-8') as f:
        json.dump(report, f, indent=2, allow_nan=False)
        f.write('\n')
    summary = [
      '# A(-2) running-UV scalar boundary pilot', '',
      f"Fixed Yukawa input: `{source_data['source']}`; topology coupling 5e-9 at 1e7 GeV.",
      f"Original UV mS^2 (TeV^2): {x0.tolist()}.",
      f"Original running scalar masses at 1 TeV: {probes[0]['scalar_threshold_mass_tev']} TeV.",
      f'Affine map closure error: {check_resid:.6g} TeV^2.',
      f'Equal-mass target: {target_tev:g} TeV for S1 and S2.',
      f"Algebraic UV target (TeV^2): {exact.get('uv_mass_sq_tev2', 'ill-conditioned')}.",
      f"Some nonnegative UV mass-squared choice reaches 0.5-5 TeV: **{unrestricted['feasible']}**.",
      f"Some UV masses also within 0.5-5 TeV: **{light['feasible']}**.",
      f'Physically screened feasible solution found: **{any_pass}**.',
      '', f"Status: `{report['status']}`.", '',
      'This conclusion is restricted to the chosen fixed Yukawas, MF, quartics and candidate one-loop RGEs.',
      'Negative UV mass-squared solutions are reported diagnostically, not adopted as benchmarks.',
      'No zero-extra-bilinear threshold corrections or general scalar mass mixing are implemented.',
    ]
    (base/'summary.md').write_text('\n'.join(summary)+'\n', encoding='utf-8')
    print('RESULT', report['status'], flush=True)
    print('REPORT', out, flush=True)
    return report


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run-id', required=True)
    ap.add_argument('--yukawa-source', default='historical',
                    choices=('historical','frozen-newfit'))
    ap.add_argument('--target-tev', type=float, default=2.)
    ap.add_argument('--probe-shift-tev', type=float, default=2.)
    args = ap.parse_args()
    run(run_id=args.run_id, source=args.yukawa_source,
        target_tev=args.target_tev, probe_shift_tev=args.probe_shift_tev)


if __name__ == '__main__':
    main()
