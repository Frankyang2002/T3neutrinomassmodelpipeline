"""Generate genuine SM+Weinberg/observable RG trajectories for the frozen-64 campaign.

Numerically integrates the repository's existing SM+Weinberg RGE from 1 TeV to
100 GeV using each *saved* C5(1 TeV), charged-lepton, and EFT1 endpoint state.
Does not invent C12/LLSS running or redo general-scalar matching.
All source files and derived outputs remain in the local working tree.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from scripts import T3Frozen64Sources as policy

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'output/studies/general_scalar_extended/full_scalar_frozen_64'
FORMAT = 't3-frozen-64-saved-sm-weinberg-trajectory-v1'
FIELDS = tuple(policy.FIELDS)


def require(condition, message):
    policy.require(bool(condition), message)


def cell_inputs(cell: dict, run_id: str):
    """Fail-closed on exact input provenance. No filesystem discovery/globbing."""
    require(cell['status'] == 'NEUTRINO_COMPLETE', 'no saved central neutrino result')
    bench, model = cell['benchmark'], cell['model']
    base = f'output/studies/general_scalar_extended/full_scalar_frozen_64/{run_id}/{bench}'
    require(cell['trajectory'] == f'{base}/{policy.tag(model)}/continuous_trajectory.json',
            'wrong scalar trajectory path')
    require(cell['neutrino'] == f'{base}/neutrino_v1/{policy.tag(model)}_neutrino.json',
            'wrong neutrino path')
    tp = policy.validate_source(ROOT, cell['trajectory'], cell['trajectory_sha256'])
    npth = policy.validate_source(ROOT, cell['neutrino'], cell['neutrino_sha256'])
    t, n = policy.read_json(tp), policy.read_json(npth)
    require(t.get('model') == model and t.get('boundary_mode') == 'FROZEN_UV'
            and t.get('settings', {}).get('benchmark') == bench
            and t.get('uv_ode_evaluations') == 0,
            'trajectory does not have frozen-UV provenance')
    require(n.get('model') == model and n.get('source_continuous_sha256') == cell['trajectory_sha256']
            and n.get('full_general_scalar_physical_c5_certified') is False
            and n.get('sm_weinberg_rge_performed') is True,
            'neutrino matching provenance/status inconsistent')
    require(t.get('scalar_endpoint_physically_admissible') is True,
            'tachyonic scalar endpoint')
    return t, n


def calculate(t: dict, n: dict, *, points: int = 41) -> dict:
    """Reuse exact low-energy solver physics and verify against central saved endpoint."""
    import numpy as np
    from scripts import RunT3FullScalarNeutrino as nu
    from Numerical.running.ScalarThresholdBoundary import (
        build_final_sm_boundary, build_sm_weinberg_initial_conditions,
    )
    from Numerical.running.WeinbergTrajectory import run_weinberg_trajectory
    from Numerical.fitting.ParameterScan import (
        neutrino_mass_from_c5, charged_lepton_mass_basis_matrix,
    )
    from Numerical.fitting.OscillationFit import prediction_from_mass_matrix

    require(3 <= points <= 201, 'points must be 3..201')
    require(t['model'] == n['model'], 'model differs across central sources')
    eft = nu.decode_eft(t['eft1_scalar_endpoint_state'])
    c5 = nu.validated_c5(nu.complex_matrix(n['c5_at_1tev_gev_inverse'], 'saved C5(1TeV)'))
    boundary = build_final_sm_boundary(eft, matching_scale_gev=1000.0)
    initial = build_sm_weinberg_initial_conditions(boundary, c5)
    scales = np.geomspace(1000., 100., points)
    flow = run_weinberg_trajectory(initial, 1000., 100., save_scales_gev=scales,
                                   rtol=1e-8, atol=1e-11)
    rows = []
    for i in range(flow.n_points):
        p = flow.point_at_index(i)
        k = nu.validated_c5(p.c5)
        flavour = charged_lepton_mass_basis_matrix(
            neutrino_mass_from_c5(k, vev_gev=246.22), p.ye)
        obs = prediction_from_mass_matrix(flavour, ordering='NO').as_dict()
        masses = np.linalg.svd(flavour, compute_uv=False) * 1e9
        require(np.isfinite(masses).all() and (masses >= 0).all(), 'invalid neutrino masses')
        require(all(math.isfinite(float(obs[key])) for key in FIELDS),
                'nonfinite neutrino observables')
        item = {'mu_gev': float(p.mu_gev),
                'c5_norm': float(np.linalg.norm(k)),
                'c5_components': {f'{a+1}{b+1}': float(abs(k[a,b]))
                                  for a in range(3) for b in range(a,3)},
                'masses_ev_descending': [float(x) for x in masses],
                'observables': {key: float(obs[key]) for key in FIELDS}}
        rows.append(item)
    require(len(rows) == points and math.isclose(rows[0]['mu_gev'], 1000.)
            and math.isclose(rows[-1]['mu_gev'], 100.), 'wrong grid endpoints')

    # Check actual complex Wilson coefficient, not just absolute values.
    k_final = nu.validated_c5(flow.final_point.c5)
    saved_final = nu.complex_matrix(n['c5_at_100gev_gev_inverse'], 'saved C5(100GeV)')
    rel = float(np.linalg.norm(k_final - saved_final)
                / max(np.linalg.norm(saved_final), 1e-300))
    require(rel < 1e-4, f'Weinberg endpoint C5 mismatch: {rel:.3g}')
    for key in FIELDS:
        a, b = rows[-1]['observables'][key], float(n['oscillation_observables_normal_ordering'][key])
        require(math.isclose(a, b, rel_tol=3e-4, abs_tol=1e-9),
                f'100 GeV saved observable parity failure: {key}: {a} != {b}')
    for a, b in zip(rows[-1]['masses_ev_descending'], n['takagi_masses_ev_descending']):
        require(math.isclose(a, float(b), rel_tol=3e-4, abs_tol=1e-8),
                'saved Takagi mass parity failure')
    return {'status': 'VerifiedSMWeinbergEndPointParity',
            'model': n['model'], 'matching_status': n['status'],
            'eft1_operator_c12_status': 'NOT_SAVED_OR_INTEGRATED_FOR_THIS_FROZEN_RUN',
            'sm_weinberg_c5_status': 'NumericallyIntegratedAgainstSavedEndpoint',
            'physical_full_general_scalar_matching_certified': False,
            'points': rows, 'c5_endpoint_relative_residual': rel,
            'solver_nfev': int(flow.nfev)}


def generate(run_id: str, output_id: str, *, points: int = 41,
             resume: bool = False, bench_filter: str = 'all',
             model_filter: str = 'all', plan: bool = False) -> dict:
    require(bool(policy.RUN_ID.fullmatch(run_id)) and bool(policy.RUN_ID.fullmatch(output_id)),
            'unsafe ID')
    require(3 <= points <= 201, 'points must be 3..201')
    src = DATA / run_id
    mp = src / 'manifest.json'
    manifest = policy.validate_manifest(policy.read_json(mp))
    generator_hash = policy.digest(Path(__file__))
    selected = [cell for cell in manifest['cells']
                if (bench_filter == 'all' or cell['benchmark'] == bench_filter)
                and (model_filter == 'all' or cell['model'] in model_filter.split(','))]
    target = src / output_id
    tally = {'requested': len(selected), 'complete': 0, 'blocked': 0,
             'not_calculated': 0, 'verified_resumed': 0}
    for cell in selected:
        b, m = cell['benchmark'], cell['model']
        if cell['status'] != 'NEUTRINO_COMPLETE':
            tally['not_calculated'] += 1
            print(f'{b}/{m}: no central neutrino result')
            continue
        try:
            t, n = cell_inputs(cell, run_id)
            if plan:
                tally['complete'] += 1
                continue
            output = target / b / (policy.tag(m) + '.json')
            inputs = {'trajectory_sha256': cell['trajectory_sha256'],
                      'neutrino_sha256': cell['neutrino_sha256'],
                      'generator_sha256': generator_hash,
                      'points_requested': points}
            if output.exists():
                require(resume, 'output exists: ' + str(output))
                old = policy.read_json(output)
                require(old.get('source') == inputs and
                        old.get('status') == 'VerifiedSMWeinbergEndPointParity',
                        'cached trajectory incompatible; use new output ID')
                tally['verified_resumed'] += 1
                continue
            item = calculate(t, n, points=points)
            item['benchmark'] = b
            item['source'] = inputs
            output.parent.mkdir(parents=True, exist_ok=True)
            from scripts import RunT3FullScalarNeutrino as nu
            nu.atomic_json(output, item)
            tally['complete'] += 1
            print(f'{b}/{m}: verified {points}-point Weinberg + neutrino trajectory')
        except (OSError, ValueError, KeyError, TypeError, RuntimeError, ArithmeticError) as error:
            tally['blocked'] += 1
            print(f'BLOCKED {b}/{m}: {type(error).__name__}: {error}', flush=True)
    if not plan:
        from scripts import RunT3FullScalarNeutrino as nu
        nu.atomic_json(target / 'summary.json', {
            'format': FORMAT, 'run_id': run_id, 'output_id': output_id,
            'manifest_sha256': policy.digest(mp), 'source_generator_sha256': generator_hash,
            'scope': 'SM+Weinberg from 1 TeV to 100 GeV only; C12 not solved',
            'tally_current_command': tally,
        })
    print(json.dumps({'plan': plan, 'counts': tally}, indent=2))
    return tally


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-id', default='frozen_general64_v1')
    p.add_argument('--output-id', default='weinberg_running_v1')
    p.add_argument('--points', type=int, default=41)
    p.add_argument('--benchmark', choices=('all', *policy.BENCHMARKS), default='all')
    p.add_argument('--models', default='all', help='Comma-separated model IDs, or all')
    p.add_argument('--resume', action='store_true')
    p.add_argument('--plan', action='store_true')
    a = p.parse_args(argv)
    if a.models != 'all':
        ms = a.models.split(',')
        require(len(ms) == len(set(ms)) and all(x in policy.MODELS for x in ms),
                'invalid/duplicate models')
    tally = generate(a.run_id, a.output_id, points=a.points,
                     resume=a.resume, bench_filter=a.benchmark,
                     model_filter=a.models, plan=a.plan)
    return 2 if tally['blocked'] else 0


if __name__ == '__main__':
    raise SystemExit(main())
