"""Compute 16 x 4 candidate *general-scalar* frozen-UV T3 backgrounds.

This is a new frozen boundary on the existing FULL scalar EFT1 engine. It is
NOT the old bulk48 general-quartic evaluator. It uses historical fitted Yukawa
matrices as named seeds (NOT new full-general optimisations) and a clearly
specified illustrative Stage29 extra-scalar boundary at 100 TeV.

For each frozen configuration: full ordinary + all permitted extra quartics
and allowed nonzero holomorphic scalar bilinear are evolved through EFT1.
Optional C5/SM-Weinberg continuation reuses RunT3FullScalarNeutrino and retains
all its explicitly conditional physics flags. No UV ODE is run.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from dataclasses import replace
import json
import math
from pathlib import Path
import numpy as np

from scripts import T3Frozen64Sources as policy
from scripts import RunT3FullScalarContinuous as full
from scripts import RunT3FullScalarNeutrino as neutrino
from scripts import RunT3GeneralScalar48Campaign as historical

ROOT = full.ROOT
OUT = ROOT / 'output/studies/general_scalar_extended/full_scalar_frozen_64'


def inputs(benchmark: str, model: str, *, stage22: Path, betas: Path) -> dict:
    strategy, case = policy.BENCHMARKS[benchmark]
    registry = full.load_registry(full.REGISTRY)
    entry, new, reused = historical.inventory(registry, model)
    original = historical.historical(case, model, entry)
    donor = historical.historical(case, 'A_-2', registry['models']['A_-2']) if strategy == 'shared_a_m2' else None
    cfg = deepcopy(original['cfg'])
    if donor is not None:
        # A(-2) donor from the SAME high/low frozen scalar scenario.
        for key in ('y1', 'y2'):
            cfg['base_state']['ordinary'][key] = deepcopy(donor['cfg']['base_state']['ordinary'][key])
    stage29_file = stage22 / 'stage29_convention_uv_handoff.json'
    stage29 = full.read(stage29_file)
    full.s30.verify_stage29_registry_provenance(stage29)
    values = full.s31.stage29_values(stage29, entry, new, reused)
    vec, layout = full.s31.prepare_state(cfg, entry, reused, values)
    base = full.unpack_uv_state(vec, layout, mu_gev=full.s30.START).validated()
    # Frozen UV: impose the SAME saved historical state at muF; no UV RG.
    base = replace(base, mu_gev=full.s30.STOP).validated()
    scalar = full.GeneralT3ScalarRGState(base, {n: values[n] for n in new}).validated(registry)
    kind = full.quadratic_kind(full.load_manifest(full.MANIFEST), model)
    beta_file = betas / f'{full.model_tag(model)}_eft1_quadratic_live_rgbeta.json'
    sources = {
        'registry': full.REGISTRY,
        'manifest': full.MANIFEST,
        'stage29_boundary': stage29_file,
        'historical_fitted_config': original['cfg_path'],
        'historical_fit': original['fit_path'],
        'eft1_mass_beta': beta_file,
        'quadratic_rgbeta_source': ROOT / 'RGE/running/rgbeta/RunGeneralT3QuadraticRGBeta.wl',
        'continuous_runner': Path(full.__file__),
        'neutrino_runner': Path(neutrino.__file__),
        'frozen_runner': Path(__file__).resolve(),
    }
    if donor is not None:
        sources['same_scenario_a_m2_donor_config'] = donor['cfg_path']
        sources['same_scenario_a_m2_donor_fit'] = donor['fit_path']
    hashes = {policy.relative(ROOT, p): policy.digest(p) for p in sources.values()}
    return {'registry': registry, 'scalar': scalar, 'kind': kind,
            'beta_file': beta_file, 'hashes': hashes, 'case': case,
            'strategy': strategy, 'original_config': original['cfg_path'],
            'donor_config': donor['cfg_path'] if donor else None}


def settings_for(model: str, benchmark: str, mu2_fraction: float, points: int) -> dict:
    return {'model': model, 'benchmark': benchmark, 'frozen_uv': True,
            'scalar_boundary_kind': 'historical_lambda_plus_stage29_illustrative_extras',
            'historical_yukawas_reoptimised_for_full_scalar': False,
            'mu2_fraction': float(mu2_fraction), 'mu2_phase_rad': 0.0,
            'saved_points_per_segment': int(points),
            'uv_start_gev': float(full.s30.START),
            'fermion_threshold_gev': float(full.s30.STOP),
            'scalar_endpoint_gev': float(full.s34.STOP)}


def ordinary_quartics(state):
    doc = full.s33.json_state(state)
    return {name: value for name, value in doc.items()
            if name.startswith("lambda") and value is not None}


def frozen_rows(base, extras, mu2: complex, points: int):
    scales = np.geomspace(float(full.s30.START), float(full.s30.STOP), points)
    row = {'mS1Sq_gev2': float(base.mS1Sq), 'mS2Sq_gev2': float(base.mS2Sq),
           'extra_quartics': {k: full.cjson(v) for k, v in extras.items()},
           'ordinary_quartics': ordinary_quartics(base),
           'mu2_gev2': full.cjson(mu2)}
    return [{'mu_gev': float(scale), **row} for scale in scales]


def calculate(benchmark: str, model: str, dest: Path, *, mu2_fraction: float,
              points: int, stage22: Path, betas: Path, resume: bool) -> dict:
    spec = inputs(benchmark, model, stage22=stage22, betas=betas)
    settings = settings_for(model, benchmark, mu2_fraction, points)
    target = dest / 'continuous_trajectory.json'
    if target.exists():
        policy.require(resume, 'existing output; use --resume: ' + str(target))
        policy.require(full.result_is_resumable(target, expected_inputs=spec['hashes'],
                      expected_settings=settings), 'incompatible existing result')
        return full.read(target)
    kind = spec['kind']
    registry = spec['registry']
    frozen = spec['scalar']
    uv = frozen.baseline
    mu2 = (mu2_fraction * min(uv.mS1Sq, uv.mS2Sq) + 0j) if kind != 'none' else 0j
    policy.require(0 <= mu2_fraction <= 0.45 and math.isfinite(mu2_fraction),
                   'invalid quadratic fraction')
    policy.require(min(uv.mS1Sq, uv.mS2Sq) > 0, 'nonpositive frozen scalar mass')
    intermediate = full.build_intermediate_scalar_boundary(uv, matching_scale_gev=full.s30.STOP)
    eft_initial = full.GeneralT3ScalarRGState(intermediate, frozen.extra).validated(registry)
    eft_y, eft_layout = full.pack_general_scalar_state(eft_initial, registry)
    if kind != 'none':
        eft_y = np.concatenate((eft_y, [mu2.real, mu2.imag]))
    beta_payload = full.read(spec['beta_file'])
    rhs = full.beta_inventory(beta_payload, model=model, stage='EFT1',
                              representation=intermediate.representation,
                              layout=eft_layout, registry=registry, kind=kind)
    eft = full.integrate_segment(eft_y, mu_start=full.s34.START,
                                 mu_end=full.s34.STOP, layout=eft_layout,
                                 registry=registry, stage='EFT1', kind=kind,
                                 betas=rhs, points=points)
    endpoint = full.unpack_general_scalar_state(eft.y[:eft_layout.size, -1], eft_layout,
                                                 registry, uv=False, mu_gev=full.s34.STOP)
    mu2_final = (complex(float(eft.y[eft_layout.size, -1]),
                         float(eft.y[eft_layout.size+1, -1])) if kind != 'none' else 0j)
    spectrum = full.quadratic_mass_hessian(full.load_manifest(full.MANIFEST), model,
                          endpoint.baseline.mS1Sq, endpoint.baseline.mS2Sq, mu2_final)
    policy.require(spectrum['analytic_spectrum_verified'] is True,
                   'quadratic spectrum internal check failed')
    dest.mkdir(parents=True, exist_ok=True)
    sidecar = dest / f'{full.model_tag(model)}_eft1_quadratic_running.json'
    side = {
        'format': full.qeft.FORMAT,
        'status': 'CandidateQuadraticAndFullQuarticEFT1Integrated',
        'model': model, 'quadratic_kind': kind,
        'mu_initial_gev': float(full.s34.START), 'mu_final_gev': float(full.s34.STOP),
        'mu2_beta_included': kind != 'none',
        'mu2_initial_gev2': full.cjson(mu2), 'mu2_final_gev2': full.cjson(mu2_final),
        'initial_baseline': full.s33.json_state(eft_initial.baseline),
        'final_baseline': full.s33.json_state(endpoint.baseline),
        'initial_new_quartics': {k: full.cjson(v) for k, v in eft_initial.extra.items()},
        'final_new_quartics': {k: full.cjson(v) for k, v in endpoint.extra.items()},
        'physical_matchete_rgbeta_mapping_certified': False,
        'scalar_mass_threshold_matching_performed': False,
        'wilson_coefficient_matching_performed': False,
        'physical_neutrino_mass_prediction': False,
        'continuous_uv_quadratic_running_performed': False,
        'continuous_uv_fermion_threshold_eft1': True,
        'input_sha256': spec['hashes'], 'settings': settings,
    }
    policy.require(not sidecar.exists(), 'orphaned pre-existing sidecar: ' + str(sidecar))
    full.write_json_atomic(sidecar, side)
    eft_rows = full.trajectory_rows(eft, layout=eft_layout, registry=registry,
                                    stage='EFT1', kind=kind)
    for i, row in enumerate(eft_rows):
        state = full.unpack_general_scalar_state(eft.y[:eft_layout.size, i],
                    eft_layout, registry, uv=False, mu_gev=row['mu_gev'])
        row['ordinary_quartics'] = ordinary_quartics(state.baseline)
    payload = {
        'format': full.FORMAT, 'status': full.STATUS,
        'model': model, 'settings': settings, 'input_sha256': spec['hashes'],
        'boundary_mode': 'FROZEN_UV', 'yukawa_mode': spec['strategy'],
        'yukawa_seed_status': 'historical_frozen_fit_not_full_general_refit',
        'additional_quartic_boundary_status': 'Stage29 illustrative values imposed at muF',
        'quadratic_kind': kind,
        'quadratic_mu2_uv_gev2': full.cjson(mu2),
        'quadratic_mu2_at_fermion_threshold_gev2': full.cjson(mu2),
        'quadratic_mu2_at_scalar_endpoint_gev2': full.cjson(mu2_final),
        'new_extra_quartics_evolved': len(frozen.extra),
        'reused_extra_quartics_evolved': len(registry['models'][model]['preexisting_rgbeta_symbols']),
        'extra_quadratic_terms_evolved': int(kind != 'none'),
        'uv_ode_evaluations': 0, 'eft1_ode_evaluations': int(eft.nfev),
        'heavy_fermion_threshold_masses_gev': full.fermion_threshold_masses(uv).masses_gev.tolist(),
        'uv_trajectory': frozen_rows(uv, frozen.extra, mu2, points),
        'eft1_trajectory': eft_rows,
        'uv_fermion_threshold_state': full.s33.json_state(uv),
        'eft1_scalar_endpoint_state': full.s33.json_state(endpoint.baseline),
        'scalar_endpoint_physically_admissible': not spectrum['tachyonic'],
        'scalar_endpoint_spectrum': {k: v for k, v in spectrum.items()
                                    if k != 'real_scalar_hessian_gev2'},
        'saved_candidate_eft1_endpoint': sidecar.name,
        'mass_bilinear_tree_threshold_prescription': 'frozen at muF, then candidate EFT1 running',
        'physical_matchete_rgbeta_mapping_certified': False,
        'fermion_threshold_wilson_matching_performed': False,
        'scalar_threshold_wilson_matching_performed': False,
        'full_sequential_neutrino_mass_prediction': False,
        'limitations': [
            'Stage29 extra-scalar couplings are illustrative and reimposed at muF, not fitted',
            'Historical optimized/shared A(-2) Yukawas are seeds, not full-general refits',
            'Frozen UV segment is a held boundary, not a solved UV RGE',
            'Full Matchete/RGBeta scalar normalization and sequential Wilson matching uncertified',
            'Six A/D/E self-bilinear C5 corrections absent from existing neutrino kernel',
        ],
    }
    full.write_json_atomic(target, payload)
    return payload


def run(args) -> int:
    policy.require(bool(policy.RUN_ID.fullmatch(args.run_id)), 'unsafe run ID')
    policy.require(4 <= args.points <= 200, 'points outside range 4..200')
    selected_benchmarks = list(policy.BENCHMARKS) if args.benchmark == 'all' else [args.benchmark]
    selected_models = list(policy.MODELS) if args.models == 'all' else args.models.split(',')
    policy.require(len(set(selected_models)) == len(selected_models) and
                   all(m in policy.MODELS for m in selected_models), 'unknown/duplicate model selection')
    root = OUT / args.run_id
    if root.exists() and not args.resume and not args.plan:
        raise ValueError('Output run directory already exists: ' + str(root))
    old = None
    if args.resume and (root / 'manifest.json').is_file():
        old = policy.validate_manifest(policy.read_json(root / 'manifest.json'))
    cells = []
    errors = 0
    for benchmark, (strategy, case) in policy.BENCHMARKS.items():
        for model in policy.MODELS:
            output = root / benchmark / full.model_tag(model)
            cell = {'benchmark': benchmark, 'model': model, 'case': case,
                    'yukawa_mode': strategy, 'status': 'NOT_RUN',
                    'trajectory': None, 'trajectory_sha256': None,
                    'neutrino': None, 'neutrino_sha256': None,
                    'matching_status': None, 'error': None}
            if old:
                prev = next(x for x in old['cells'] if (x['benchmark'], x['model']) == (benchmark,model))
                cell.update(prev)
            if benchmark in selected_benchmarks and model in selected_models:
                try:
                    spec = inputs(benchmark, model, stage22=args.stage22, betas=args.quadratic_betas)
                    cell['source_sha256'] = spec['hashes']
                    cell['source_config'] = policy.relative(ROOT, spec['original_config'])
                    cell['donor_config'] = policy.relative(ROOT, spec['donor_config']) if spec['donor_config'] else None
                    if not args.plan:
                        result = calculate(benchmark, model, output, mu2_fraction=args.mu2_fraction,
                            points=args.points, stage22=args.stage22, betas=args.quadratic_betas,
                            resume=args.resume)
                        path = output / 'continuous_trajectory.json'
                        cell['status'] = 'SCALAR_COMPLETE'
                        cell['trajectory'] = policy.relative(ROOT, path)
                        cell['trajectory_sha256'] = policy.digest(path)
                        cell['scalar_endpoint_admissible'] = result['scalar_endpoint_physically_admissible']
                        if not args.skip_neutrino and result['scalar_endpoint_physically_admissible']:
                            try:
                                nroot = root / benchmark / 'neutrino_v1'
                                n = neutrino.compute_one(model, source_root=root/benchmark,
                                                           destination_root=nroot, resume=args.resume)
                                npth = nroot / (full.model_tag(model) + '_neutrino.json')
                                cell['status'] = 'NEUTRINO_COMPLETE'
                                cell['neutrino'] = policy.relative(ROOT, npth)
                                cell['neutrino_sha256'] = policy.digest(npth)
                                cell['matching_status'] = n['status']
                            except (KeyError, ValueError, TypeError, RuntimeError, OSError,
                                    OverflowError, ArithmeticError) as exc:
                                # Preserve valid scalar RGE data if only C5 is blocked.
                                errors += 1
                                cell['neutrino_error'] = f'{type(exc).__name__}: {exc}'[:500]
                                print(f'NEUTRINO BLOCKED {benchmark}/{model}: '
                                      f'{cell["neutrino_error"]}', flush=True)
                except (KeyError, ValueError, TypeError, RuntimeError, OSError,
                        OverflowError, ArithmeticError) as exc:
                    errors += 1
                    cell['status'] = 'BLOCKED'
                    cell['error'] = f'{type(exc).__name__}: {exc}'[:500]
                    print(f'BLOCKED {benchmark}/{model}: {cell["error"]}', flush=True)
                print(f'{benchmark}/{model}: {cell["status"]}', flush=True)
            cells.append(cell)
    doc = {'format': policy.FORMAT, 'run_id': args.run_id,
           'case_definitions': {k: list(v) for k,v in policy.BENCHMARKS.items()},
           'physics_scope': 'candidate full-scalar EFT1; historical Yukawa seeds; conditional C5',
           'stage29_extra_scalar_at_muF': True,
           'models': list(policy.MODELS), 'cells': cells}
    policy.validate_manifest(doc)
    if not args.plan:
        root.mkdir(parents=True, exist_ok=True)
        full.write_json_atomic(root / 'manifest.json', doc)
        print('Manifest:', root / 'manifest.json')
    else:
        print(json.dumps({'cells':len(cells),'selected':len(selected_benchmarks)*len(selected_models),
                          'blocked':errors,'status':'READ_ONLY_SOURCE_PREFLIGHT'},indent=2))
    return 2 if errors else 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-id', default='frozen_general64_v1')
    p.add_argument('--benchmark', choices=('all',*policy.BENCHMARKS), default='all')
    p.add_argument('--models', default='all', help='Comma-separated exact model IDs, or all')
    p.add_argument('--stage22', type=Path, default=full.STAGE22)
    p.add_argument('--quadratic-betas', type=Path, default=full.BETA_DIR)
    p.add_argument('--mu2-fraction', type=float, default=0.01)
    p.add_argument('--points', type=int, default=16)
    p.add_argument('--plan', action='store_true', help='Validate exactly named sources; no writes')
    p.add_argument('--resume', action='store_true')
    p.add_argument('--skip-neutrino', action='store_true')
    return run(p.parse_args(argv))


if __name__ == '__main__':
    raise SystemExit(main())
