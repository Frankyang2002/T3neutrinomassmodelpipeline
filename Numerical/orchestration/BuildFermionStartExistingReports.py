"""Fermion-start data adapters for the repository's EXISTING HTML report builders.

RGE: Numerical.plotting.InteractiveModelComparison.write_dashboard
Sensitivity ovals: Numerical.plotting.BuildFittedModelEnvelopes.build

Only the 011 (UV frozen, EFT1+Weinberg running) numerical prescription is used.
Existing fit/config/physics code is not modified. New saved data are namespaced.
"""
from __future__ import annotations

import argparse
import copy
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from Numerical.orchestration import RunFermionStartLambdaStudy as study
from Numerical.orchestration import RunSmallLambdaRefit as base
from Numerical.fitting.BenchmarkSensitivity import gamma_grid, scaled_config_payload
from Numerical.fitting.ParameterScan import neutrino_mass_from_c5, charged_lepton_mass_basis_matrix
from Numerical.fitting.OscillationFit import evaluate_mass_matrix_fit
from physics.NeutrinoObservables import calculate_neutrino_observables
from Numerical.plotting.InteractiveModelComparison import write_dashboard
from Numerical.plotting.BuildFittedModelEnvelopes import build as build_ovals
from Numerical.running.ScalarThresholdBoundary import build_sm_weinberg_initial_conditions
from Numerical.running.SMWeinbergEvolution import _pack, _unpack, _beta
from Numerical.diagnostics.IntermediateWeinbergDiagnostics import sample_intermediate_direct_weinberg
from Numerical.diagnostics.FinalC5ContributionDiagnostics import evaluate_final_c5_contributions

MODES = ('optimised', 'shared_a_m2')
CASES = ('5e-9', '5e-10')
DIRECTIONS = ('T1', 'T2', 'T3')
OBS = tuple(study.FIELDS)
DESCRIPTIONS = {
    'T1': 'lambdaT3 -> gamma lambdaT3; y1,y2 -> gamma^(-1/2) y1,y2',
    'T2': 'lambdaT3 -> gamma lambdaT3; y2 -> gamma^(-1) y2',
    'T3': 'y1 -> gamma y1; y2 -> gamma^(-1) y2',
}
ROOT = base.ROOT


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    tmp.replace(path)


def source_path(case, mode, model):
    return study.ROOT / case / 'comparison' / mode / 'models' / model / 'fermion_start_config.json'


def evaluator(case, mode, model, mapping):
    path = source_path(case, mode, model)
    if not path.is_file():
        raise FileNotFoundError(f'Run compare stage first: {path}')
    cfg = base.load(path)
    return cfg, study.FermionStartEvaluator(cfg, study.artifacts_for(model, mapping))


def selection(case, mode, model):
    return {'case': case, 'mode': mode, 'model': model,
            'boundary': 'mu_F', 'uv_running': False,
            'fermion_matching': True, 'eft1_running': True,
            'scalar_matching': True, 'weinberg_running': True,
            'source_config': str(source_path(case, mode, model).relative_to(ROOT))}


def prediction(fit):
    return dict(zip(fit.observable_names, map(float, fit.prediction)))


def assert_saved_endpoint(case, mode, model, fitted):
    original = study.ROOT / case / 'comparison' / mode / 'models' / model / 'evaluation.json'
    ref = base.load(original)
    if ref.get('status') != 'Success':
        raise ValueError(f'{model}: saved comparison not successful')
    maximum = max(abs(float(fitted[k]) - float(ref['prediction'][k])) for k in OBS)
    if maximum > 1e-7 * max(1.0, *(abs(float(ref['prediction'][k])) for k in OBS)):
        raise RuntimeError(f'{model}: endpoint differs from saved comparison by {maximum}')
    return ref


def rge_one(case, mode, model, mapping, n_scales):
    cfg, ev = evaluator(case, mode, model, mapping)
    fit, renorm, full = ev.evaluate(study.vector_values(cfg), return_trajectory=True)
    ref = assert_saved_endpoint(case, mode, model, prediction(fit))
    mu_s, mu_low = ev.mu_s, ev.mu_low
    mu = np.geomspace(mu_s, mu_low, n_scales)
    initial = build_sm_weinberg_initial_conditions(renorm.final_sm_boundary,
                                                    full.c5_at_scalar_threshold)
    solution = solve_ivp(_beta, (math.log(mu_s), math.log(mu_low)), _pack(initial),
                         t_eval=np.log(mu), method='DOP853', rtol=ev.rtol, atol=ev.atol)
    if not solution.success:
        raise RuntimeError(f'Weinberg trajectory failed: {solution.message}')
    curves = {k: [] for k in OBS}
    c5_abs = []
    masses_ev = []
    for col in solution.y.T:
        _, _, _, _, ye, _, _, K = _unpack(col)
        K = (K + K.T) / 2
        mass = neutrino_mass_from_c5(K, vev_gev=ev.vev)
        flavour = charged_lepton_mass_basis_matrix(mass, ye)
        fitted = evaluate_mass_matrix_fit(flavour, ev.target)
        # Use the same Takagi / mass-ordering routine as the canonical pipeline.
        physical = calculate_neutrino_observables(flavour, ordering=ev.target.ordering)
        mass_eigenvalues = np.asarray(physical.masses_ev, dtype=float)
        if mass_eigenvalues.shape != (3,) or not np.all(np.isfinite(mass_eigenvalues)):
            raise ValueError(f"{model}: invalid neutrino mass eigenvalues at sampled scale")
        masses_ev.append(mass_eigenvalues.tolist())
        p = prediction(fitted)
        for k in OBS:
            curves[k].append(p[k])
        c5_abs.append(np.abs(K).tolist())
    for k in OBS:
        actual = float(ref['prediction'][k]); predicted = float(curves[k][-1])
        if abs(actual - predicted) > 3e-5 * max(abs(actual), 1e-10):
            raise RuntimeError(f'{model}: sampled Weinberg endpoint mismatch in {k}: {predicted} vs {actual}')
    # EFT1 running was saved by the actual fermion-start integrator.
    eft = renorm.eft1
    intermediate = {'mu_gev': [float(x) for x in eft.mu_gev],
                    'state_vectors': np.asarray(eft.y.T, dtype=float).tolist(),
                    'layout_blocks': [{'name': b.name, 'kind': b.kind,
                                       'start': int(b.start), 'stop': int(b.stop)}
                                      for b in eft.layout.blocks]}
    # Original fixed-order LLSS -> O5 contribution, not a coupled EFT RGE.
    # Sample the same symbolic K_running expression used by the original pipeline.
    artifacts = study.artifacts_for(model, mapping)
    intermediate_mu = np.geomspace(ev.mu_f, ev.mu_s, n_scales)
    direct = sample_intermediate_direct_weinberg(artifacts[2], renorm, intermediate_mu)
    contributions = evaluate_final_c5_contributions(artifacts[2], renorm)
    combined = np.asarray(contributions.combined, dtype=complex)
    direct_at_s = np.asarray(direct.delta_c5[-1], dtype=complex)
    hard_at_s = np.asarray(contributions.direct_running, dtype=complex)
    # Use amplitude-scaled residuals: absolute 1e-14 tolerances would mask
    # relative differences in neutrino dimension-five coefficients.
    scale = max(float(np.max(np.abs(combined))), 1e-40)
    if np.max(np.abs(combined - full.c5_at_scalar_threshold)) / scale > 1e-7:
        raise RuntimeError(f'{model}: symbolic C5 decomposition does not match fit threshold')
    if np.max(np.abs(direct_at_s - hard_at_s)) / scale > 1e-7:
        raise RuntimeError(f'{model}: intermediate direct C5 endpoint does not match symbolic threshold')
    intermediate_direct = {
        'mu_gev': direct.mu_gev.tolist(),
        'delta_c5_abs': direct.delta_c5_abs.tolist(),
        'delta_c5_real': direct.delta_c5.real.tolist(),
        'delta_c5_imag': direct.delta_c5.imag.tolist(),
        'scheme': 'fixed_order_one_loop_direct_LLSS_to_Weinberg',
        'note': 'Direct fixed-order LLSS -> O5 contribution; not a full intermediate Wilson RGE integration.',
    }
    # Matching at mu_F and mu_S remains intact. No fabricated UV curve.
    diagnostic = {'status': 'Success', 'source': 'FermionStartExistingReports',
                  'fermion_start': selection(case, mode, model),
                  'scales_gev': {'mu_fermion_threshold': ev.mu_f,
                                 'mu_scalar_threshold': mu_s, 'mu_low': mu_low},
                  'intermediate_running': intermediate,
                  'intermediate_direct_weinberg': intermediate_direct,
                  'c5_threshold_contributions': contributions.as_json_dict(),
                  'scalar_threshold': {'c5': {'abs': np.abs(full.c5_at_scalar_threshold).tolist()}},
                  'final_running': {'mu_gev': list(map(float, mu)),
                                    'masses_ev': masses_ev,
                                    'c5_abs': c5_abs, **curves},
                  'low_energy': {'chi2': float(fit.chi2), 'prediction': prediction(fit)},
                  'notes': {'uv_running': 'Disabled: no UV curve exists',
                            'intermediate_c5': 'Fixed-order LLSS -> O5 mixing from symbolic K_running; not coupled Wilson RGE evolution'}}
    out = study.ROOT / case / 'reports_existing' / mode / 'models' / model / 'running_diagnostics.json'
    write_json(out, diagnostic)
    print(case, mode, model, 'RGE saved', flush=True)


def sensitivity_one(case, mode, model, mapping, gamma_values):
    cfg, ev = evaluator(case, mode, model, mapping)
    baseline = ev.evaluate(study.vector_values(cfg))
    assert_saved_endpoint(case, mode, model, prediction(baseline))
    directions = {}
    for direction in DIRECTIONS:
        points = []
        for gamma in gamma_values:
            point = {'gamma': float(gamma)}
            try:
                varied = scaled_config_payload(cfg, direction, float(gamma))
                current = study.FermionStartEvaluator(varied, study.artifacts_for(model, mapping))
                fit = current.evaluate(study.vector_values(varied))
                point.update(status='Success', chi2=float(fit.chi2), prediction=prediction(fit))
            except Exception as exc:
                point.update(status='Failed', message=f'{type(exc).__name__}: {exc}')
            points.append(point)
        directions[direction] = {'definition': DESCRIPTIONS[direction], 'points': points}
    centers = [directions[d]['points'][len(gamma_values)//2] for d in DIRECTIONS]
    if any(p['status'] != 'Success' or abs(p['gamma'] - 1) > 1e-12 for p in centers):
        raise RuntimeError(f'{model}: missing gamma=1 centre')
    for d, point in zip(DIRECTIONS, centers):
        for key in OBS:
            if abs(point['prediction'][key] - prediction(baseline)[key]) > 1e-7 * max(1.0, abs(prediction(baseline)[key])):
                raise RuntimeError(f'{model}: gamma=1 differs for {d}/{key}')
    target = ev.target
    scan = {'status': 'Success', 'gamma_min': float(gamma_values[0]),
            'gamma_max': float(gamma_values[-1]), 'n_points': len(gamma_values),
            'spacing': 'log', 'n_failed_points': sum(p['status'] != 'Success'
                         for d in DIRECTIONS for p in directions[d]['points']),
            'directions': directions, 'selection': selection(case, mode, model),
            'target': {'ordering': target.ordering,
                       'observable_names': list(target.observable_names),
                       'central_values': target.central_values.tolist(),
                       'covariance': target.covariance.tolist(),
                       'source': target.source,
                       'note': 'Gaussian approximation, not official joint NuFIT contours.'}}
    out = study.ROOT / case / 'reports_existing' / mode / 'sensitivity' / model / 'sensitivity_diagnostics.json'
    write_json(out, scan)
    print(case, mode, model, 'sensitivity saved, failed points=', scan['n_failed_points'], flush=True)
    return scan


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case', choices=(*CASES, 'all'), default='all')
    p.add_argument('--mode', choices=(*MODES, 'all'), default='all')
    p.add_argument('--kind', choices=('rge', 'sensitivity', 'all'), default='all')
    p.add_argument('--model', help='Short or full model key (for pilot tests)')
    p.add_argument('--gamma-min', type=float, default=.1)
    p.add_argument('--gamma-max', type=float, default=10)
    p.add_argument('--points', type=int, default=11)
    p.add_argument('--rge-points', type=int, default=41)
    p.add_argument('--artifact-map', type=Path)
    p.add_argument('--force', action='store_true', help='Recompute saved diagnostics, including intermediate LLSS to O5 terms')
    args = p.parse_args()
    if args.rge_points < 3:
        p.error('--rge-points must be >= 3')
    gammas = gamma_grid(args.gamma_min, args.gamma_max, args.points)
    mapping = base.load(args.artifact_map) if args.artifact_map else None
    cases = CASES if args.case == 'all' else (args.case,)
    modes = MODES if args.mode == 'all' else (args.mode,)
    kinds = ('rge', 'sensitivity') if args.kind == 'all' else (args.kind,)
    rows = base.original_rows()
    if args.model:
        rows = [r for r in rows if args.model in (r['model'], study.shorten(r['model']))]
        if not rows:
            p.error(f'Unknown model {args.model}')
    errors = []
    for case in cases:
        for mode in modes:
            source = study.ROOT / case / 'reports_existing' / mode
            for kind in kinds:
                for row in rows:
                    model = row['model']
                    diagnostic = (source / 'models' / model / 'running_diagnostics.json' if kind == 'rge'
                                  else source / 'sensitivity' / model / 'sensitivity_diagnostics.json')
                    if diagnostic.exists() and not args.force:
                        print('Reuse', diagnostic, flush=True)
                        continue
                    try:
                        if kind == 'rge':
                            rge_one(case, mode, model, mapping, args.rge_points)
                        else:
                            sensitivity_one(case, mode, model, mapping, gammas)
                    except Exception as exc:
                        errors.append((case, mode, model, kind, f'{type(exc).__name__}: {exc}'))
                        print('FAILED', *errors[-1], flush=True)
            report_dest = ROOT / 'Reports' / 'output' / 'display_fermion_start' / case / 'reused' / mode
            report_dest.mkdir(parents=True, exist_ok=True)
            if 'rge' in kinds:
                try:
                    write_dashboard(source / 'models', report_dest / 'interactive_rge.html',
                                    f'T3 Fermion Start — UV Off — {case} — {mode}')
                except Exception as exc:
                    errors.append((case,mode,'ALL','RGE dashboard',str(exc)))
                    print('DASHBOARD FAILED', errors[-1],flush=True)
            if 'sensitivity' in kinds:
                sens = source / 'sensitivity'
                status = []
                for row in base.original_rows():
                    model = row['model']
                    item = sens / model / 'sensitivity_diagnostics.json'
                    if not item.exists():
                        status.append({'model': model, 'status': 'MissingSavedScan'});continue
                    scan = base.load(item)
                    center = scan['directions']['T1']['points'][scan['n_points']//2]
                    status.append({'model':model, 'status':'Success', 'verified_chi2':center['chi2'],
                                   'fit_stage':f'fermion_start_{mode}'})
                write_json(sens / 'comparison_summary.json', status)
                try:
                    build_ovals(sens, report_dest / 'sensitivity_ovals.html')
                except Exception as exc:
                    errors.append((case,mode,'ALL','oval dashboard',str(exc)))
                    print('OVAL DASHBOARD FAILED', errors[-1],flush=True)
    if errors:
        print('ERRORS:', len(errors), flush=True)
        for e in errors:print(' ',e,flush=True)
        raise SystemExit(1)
    print('Generated reused HTML reports under Reports/output/display_fermion_start/*/reused', flush=True)

if __name__ == '__main__':
    main()
