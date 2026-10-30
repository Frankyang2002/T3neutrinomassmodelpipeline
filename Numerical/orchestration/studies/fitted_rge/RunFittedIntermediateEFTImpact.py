"""Fixed-benchmark 2x2 ablation of intermediate EFT and final SMEFT running.

All configurations retain UV running and both threshold-matching steps.
The intermediate EFT state and/or final Weinberg coefficient are frozen in
selected branches. A common fully evolved charged-lepton Yukawa basis is used
for ALL four predictions, following RunFittedStagedRGComparison conventions.

Cases: 00 (both frozen), 10 (intermediate on), 01 (SM Weinberg on),
11 (both on). These are numerical RGE switch-off diagnostics, not separate
perturbative matching prescriptions. Outputs are written separately.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import csv
import json
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from Numerical.orchestration.support.ResumePipelineNumericalResults import PROJECT_ROOT, parse_model_key
from Numerical.orchestration.studies.fitted_rge.RunFittedSensitivityComparison import best_candidate
from Numerical.orchestration.studies.fermion_scales.RunFittedFermionScaleScan import (
    OBS, SCENARIO, discover_fitted_models, locate_artifacts, read_json,
)
from Numerical.orchestration.studies.fitted_rge.RunFittedStagedRGComparison import (
    frozen_result, with_uv_eft_states,
)
from Numerical.orchestration.studies.fitted_rge.RunFittedWeinbergRunningComparison import predictions
from Numerical.core.RGBetaEvaluator import load_rgbeta_payload
from Numerical.fitting.ScanCLI import build_uv_state_from_config
from Numerical.fitting.OscillationFit import OscillationFitTarget
from Numerical.running.FermionThresholdBoundary import build_intermediate_scalar_boundary
from Numerical.running.FinalC5TrajectoryAdapter import FinalC5TrajectoryEvaluator
from Numerical.running.T3Trajectory import (
    run_renormalisable_t3_trajectory, continue_with_weinberg_running,
)

CASES = ('00', '10', '01', '11')
LABELS = {
    '00': 'Intermediate frozen; SMEFT Weinberg frozen',
    '10': 'Intermediate running; SMEFT Weinberg frozen',
    '01': 'Intermediate frozen; SMEFT Weinberg running',
    '11': 'Intermediate running; SMEFT Weinberg running',
}
ANGLES = ('sin2_theta12', 'sin2_theta13', 'sin2_theta23')
MASS_SPLITTINGS = ('delta_m21_sq_ev2', 'delta_m3l_sq_ev2')


def evaluate_four(payload, uv_beta, eft_beta, c5_evaluator, target):
    scales = payload['scales']
    mu_f = float(scales['mu_fermion_threshold_gev'])
    mu_s = float(scales['mu_scalar_threshold_gev'])
    mu_low = float(scales['mu_low_gev'])
    numerics = payload.get('numerical', {})
    vev = float(numerics.get('vev_gev', 246.22))
    rtol = float(numerics.get('rtol', 1e-8))
    atol = float(numerics.get('atol', 1e-11))
    cfg = deepcopy(payload)
    cfg.setdefault('scan', {'bindings': {}})
    uv_initial = build_uv_state_from_config(SimpleNamespace(raw=cfg), {})
    mu_uv = float(uv_initial.mu_gev)
    if not mu_uv > mu_f > mu_s > mu_low > 0:
        raise ValueError('Require mu_UV > mu_F > mu_S > mu_low > 0')

    full = run_renormalisable_t3_trajectory(
        uv_initial, uv_beta, eft_beta,
        mu_fermion_threshold_gev=mu_f, mu_scalar_threshold_gev=mu_s,
        rtol=rtol, atol=atol,
    )
    c5_full = np.asarray(c5_evaluator({}, full), dtype=complex)
    low_full = continue_with_weinberg_running(
        full, c5_full, mu_low, rtol=rtol, atol=atol,
    )
    common_ye = low_full.weinberg.ye

    # Keep the *same* UV trajectory; freeze only the intermediate EFT.
    eft_start = build_intermediate_scalar_boundary(
        full.uv.final_state, matching_scale_gev=mu_f,
    )
    eft_frozen = frozen_result(eft_start, mu_s, uv=False)
    frozen_intermediate = with_uv_eft_states(
        full.uv, eft_frozen, mu_uv, mu_f, mu_s,
    )
    c5_frozen = np.asarray(c5_evaluator({}, frozen_intermediate), dtype=complex)
    low_frozen_intermediate = continue_with_weinberg_running(
        frozen_intermediate, c5_frozen, mu_low, rtol=rtol, atol=atol,
    )

    return {
        '00': predictions(c5_frozen, common_ye, target, vev),
        '10': predictions(c5_full, common_ye, target, vev),
        '01': predictions(low_frozen_intermediate.weinberg.K, common_ye, target, vev),
        '11': predictions(low_full.weinberg.K, common_ye, target, vev),
    }


def angle_degrees(sin2):
    if not -1e-12 <= sin2 <= 1 + 1e-12:
        raise ValueError('sin^2(theta) outside physical interval')
    return math.degrees(math.asin(math.sqrt(min(1.0, max(0.0, sin2)))))


def compare_cases(cases):
    result = {}
    for case in CASES:
        values = cases[case]['prediction']
        derived = dict(values)
        for observable in ANGLES:
            derived['theta_' + observable.removeprefix('sin2_theta') + '_deg'] = angle_degrees(values[observable])
        result[case] = derived

    deltas = {}
    for case in ('10', '01', '11'):
        differences = {}
        for observable in MASS_SPLITTINGS:
            ref, new = result['00'][observable], result[case][observable]
            differences[observable + '_percent'] = (100*(new/ref-1)) if ref != 0 else None
        for observable in ANGLES:
            key = 'theta_' + observable.removeprefix('sin2_theta') + '_deg'
            differences[key + '_shift_deg'] = result[case][key] - result['00'][key]
        deltas[case] = differences
    # Exact interaction term: 11 - 10 - 01 + 00, in absolute observable units.
    interaction = {key: result['11'][key] - result['10'][key] - result['01'][key] + result['00'][key]
                   for key in result['00']}
    return result, deltas, interaction


def write_outputs(rows, folder):
    folder.mkdir(parents=True, exist_ok=True)
    (folder/'results.json').write_text(json.dumps(rows, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    fields = ['model', 'representation', 'fit_stage', 'artifact_set', 'config_path', 'status', 'message']
    obs = list(MASS_SPLITTINGS) + list(ANGLES)
    angles = ['theta_12_deg', 'theta_13_deg', 'theta_23_deg']
    fields += [f'{case}_{o}' for case in CASES for o in obs + angles]
    fields += [f'{case}_{o}_percent' for case in ('10','01','11') for o in MASS_SPLITTINGS]
    fields += [f'{case}_{a}_shift_deg' for case in ('10','01','11') for a in angles]
    fields += [f'interaction_{o}' for o in obs + angles]
    with (folder/'results.csv').open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            record = {key: row.get(key, '') for key in fields}
            if 'representation' in row:
                record['representation'] = json.dumps(row['representation'], separators=(',',':'))
            if row['status'] == 'Success':
                for case in CASES:
                    for o in obs + angles:
                        record[f'{case}_{o}'] = row['derived'][case][o]
                for case in ('10','01','11'):
                    for key, value in row['deltas'][case].items():
                        record[f'{case}_{key}'] = '' if value is None else value
                for o in obs + angles:
                    record[f'interaction_{o}'] = row['interaction'][o]
            writer.writerow(record)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--scenario', default=SCENARIO)
    ap.add_argument('--limit', type=int, default=0,
                    help='Maximum discovered fitted models (0: all)')
    ap.add_argument('--run', action='store_true', help='Compute; otherwise check input availability')
    ap.add_argument('--output', type=Path, default=None)
    args = ap.parse_args(argv)
    if args.limit < 0:
        ap.error('--limit must be >= 0')
    output = args.output or (PROJECT_ROOT/'Reports'/'output'/'full'/'intermediate_eft_impact'/args.scenario)
    models = discover_fitted_models(args.scenario)
    if not models:
        print('No fitted models found for scenario:', args.scenario)
        return 2
    target_path = PROJECT_ROOT/'configs'/'nufit_6_1_ic24_no.json'
    rows = []
    for model in models[:args.limit or None]:
        row = {'model': model, 'status': '', 'message': ''}
        rows.append(row)
        selected = best_candidate(args.scenario, model)
        if selected is None:
            row.update(status='MissingFit', message='No usable fitted config')
            continue
        _, stage, config_path, _ = selected
        row.update(fit_stage=stage, config_path=str(config_path))
        try:
            rep = parse_model_key(model)
            row['representation'] = list(rep)
            artifact, problem = locate_artifacts(model, rep)
            if artifact is None:
                row.update(status='MissingArtifacts', message=str(problem))
                continue
            group, paths = artifact
            row['artifact_set'] = group
            payload = read_json(config_path)
            actual = tuple(int(payload['representation'][k]) for k in ('d_s1','d_s2','d_f','alpha'))
            if actual != rep:
                raise ValueError(f'Representation mismatch: {actual} != {rep}')
            s = payload['scales']
            scales = [float(payload['base_state']['mu_gev'])] + [float(s[k]) for k in (
                'mu_fermion_threshold_gev','mu_scalar_threshold_gev','mu_low_gev')]
            if not all(a > b > 0 for a,b in zip(scales, scales[1:])):
                raise ValueError('Invalid EFT scale ordering')
            row['status'] = 'Ready'
            if args.run:
                if not target_path.is_file():
                    raise FileNotFoundError(target_path)
                uv_beta = load_rgbeta_payload(paths[0])
                eft_beta = load_rgbeta_payload(paths[1])
                c5_evaluator = FinalC5TrajectoryEvaluator(paths[2])
                target = OscillationFitTarget.from_json(target_path)
                cases = evaluate_four(payload, uv_beta, eft_beta, c5_evaluator, target)
                derived, deltas, interaction = compare_cases(cases)
                row.update(status='Success', cases=cases, derived=derived,
                           deltas=deltas, interaction=interaction)
        except Exception as exc:
            row.update(status='Failed', message=f'{type(exc).__name__}: {exc}')
        print(model, row['status'], row['message'][:160], flush=True)
    if not args.run:
        print('Dry-run only. Pass --run to evaluate. Ready:', sum(x['status']=='Ready' for x in rows))
        return 0
    write_outputs(rows, output)
    counts = {s:sum(r['status']==s for r in rows) for s in sorted({r['status'] for r in rows})}
    print('Status counts:', counts)
    print('Wrote:', output)
    return 0 if counts.get('Success',0) else 1


if __name__ == '__main__':
    raise SystemExit(main())
