"""All-16 T3 zero-additional-scalar control, on the historical frozen boundary.

Zero ALL model-registry additional quartic coordinates (both new RGBeta fields
and reused fields in the baseline state) and all allowed holomorphic quadratic
coefficients, WITHOUT zeroing original lambdaT3, diagonal masses, historical
Yukawa matrices or original ordinary scalar quartics.

Two deliberately different calculations:
  held_zero: remove the extra coordinates' beta sources as well; a restricted
             potential diagnostic, NOT the RG flow of the full general theory.
  boundary_zero: impose zero extra couplings at muF=100 TeV but evolve ALL
             full candidate beta functions; allowed terms can be regenerated.

Frozen UV (no 10 PeV -> 100 TeV integration); EFT1 100 TeV -> 1 TeV.
Reuses the existing conditional C5 and SM-Weinberg numerical evaluators.
This is NOT an independently certified sequential EFT prediction.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import replace
from hashlib import sha256
import json
from pathlib import Path
from typing import Any

import numpy as np

from scripts import RunT3FrozenFullScalar64 as frozen
from scripts import RunT3FullScalarContinuous as full
from scripts import RunT3FullScalarNeutrino as neutrino
from scripts import RunT3GeneralScalar48Campaign as historical
from scripts import T3Frozen64Sources as policy

ROOT = full.ROOT
OUTPUT = ROOT / 'output/studies/general_scalar_extended/zero_extra_all16'
FORMAT = 't3-all16-zero-extra-scalar-control-v1'
MODES = ('held_zero', 'boundary_zero')
FIELDS = tuple(policy.FIELDS)


def require(test: bool, message: str) -> None:
    if not test:
        raise ValueError('All16 zero-extra scalar control: ' + message)


def select_csv(text: str, choices: tuple[str, ...]) -> tuple[str, ...]:
    selection = tuple(x.strip() for x in text.split(','))
    require(bool(selection) and all(x in choices for x in selection) and
            len(set(selection)) == len(selection), 'invalid/duplicate selection: ' + text)
    return selection


def restrict_scalar(spec: dict, model: str) -> tuple[Any, tuple[str, ...], tuple[str, ...]]:
    """Remove only registry-classed additional interactions from a *copy*.

    The original topology coupling lambdaT3 is explicitly protected.
    """
    registry = spec['registry']
    entry, new, reused = historical.inventory(registry, model)
    initial = spec['scalar'].validated(registry)
    require(set(new) == set(initial.extra), 'new quartic inventory mismatch')
    require('lambdaT3' not in new and 'lambdaT3' not in reused,
            'lambdaT3 must not be zeroed')
    base = initial.baseline
    replacement = {}
    for name in reused:
        value = getattr(base, name, None)
        require(value is not None, 'missing reused scalar quartic ' + name)
        replacement[name] = 0j  # registry reused coordinates are complex
    copied_base = replace(base, **replacement).validated() if replacement else base
    reduced = full.GeneralT3ScalarRGState(
        copied_base, {name: 0j for name in new}).validated(registry)
    require(complex(reduced.baseline.lambdaT3) == complex(base.lambdaT3),
            'original lambdaT3 changed')
    require(reduced.baseline.mS1Sq == base.mS1Sq and
            reduced.baseline.mS2Sq == base.mS2Sq,
            'original scalar diagonal masses changed')
    return reduced, tuple(new), tuple(reused)


def zero_extra_betas(betas: dict[str, str], new: tuple[str, ...],
                     reused: tuple[str, ...], quadratic: str | None) -> dict[str, str]:
    inactive = set(new) | set(reused) | ({quadratic} if quadratic else set())
    require(inactive.issubset(betas), 'missing beta coordinates ' +
            repr(sorted(inactive - set(betas))))
    require('lambdaT3' not in inactive, 'must retain original topology beta')
    return {k: ('0' if k in inactive else v) for k, v in betas.items()}


def prepare(benchmark: str, model: str, stage22: Path, beta_dir: Path) -> dict:
    require(benchmark in ('optimised_low', 'optimised_high'),
            'use original optimised high/low fitted Yukawas')
    spec = frozen.inputs(benchmark, model, stage22=stage22, betas=beta_dir)
    reduced, new, reused = restrict_scalar(spec, model)
    hashes = dict(spec['hashes'])
    self_path = Path(__file__).resolve()
    hashes[policy.relative(ROOT, self_path)] = sha256(self_path.read_bytes()).hexdigest()
    return {**spec, 'scalar': reduced, 'new': new, 'reused': reused, 'hashes': hashes}


def parameters(benchmark: str, model: str, mode: str, points: int) -> dict:
    return {
        'model': model, 'benchmark': benchmark, 'frozen_uv': True,
        'scalar_boundary_kind': 'historical_fit_all_registry_extras_zero_at_muF',
        'zero_extra_control_mode': mode,
        'historical_yukawas_reoptimised_for_full_scalar': False,
        'mu2_fraction': 0.0, 'mu2_phase_rad': 0.0,
        'saved_points_per_segment': points,
        'uv_start_gev': float(full.s30.START),
        'fermion_threshold_gev': float(full.s30.STOP),
        'scalar_endpoint_gev': float(full.s34.STOP),
    }


def extra_endpoint_values(state: Any, new: tuple[str, ...],
                          reused: tuple[str, ...]) -> dict[str, float]:
    return {name: abs(complex(state.extra[name])) for name in new} | {
        name: abs(complex(getattr(state.baseline, name))) for name in reused}


def calculate_scalar(spec: dict, *, benchmark: str, model: str,
                     mode: str, points: int, folder: Path, resume: bool) -> dict:
    require(mode in MODES, 'invalid mode')
    cfg = parameters(benchmark, model, mode, points)
    target = folder / 'continuous_trajectory.json'
    if target.exists():
        require(resume and full.result_is_resumable(
            target, expected_inputs=spec['hashes'], expected_settings=cfg),
            'existing output/incompatible inputs; use --resume or new --run-id')
        return full.read(target)
    require(not folder.exists() or not any(folder.iterdir()),
            'partial/orphaned output, use a fresh --run-id')
    original = spec['scalar']
    uv = original.baseline
    registry = spec['registry']
    kind = spec['kind']
    mu2 = 0j
    require(min(uv.mS1Sq, uv.mS2Sq) > 0, 'non-positive original scalar masses')
    eft_base = full.build_intermediate_scalar_boundary(
        uv, matching_scale_gev=full.s30.STOP)
    initial = full.GeneralT3ScalarRGState(eft_base, original.extra).validated(registry)
    y0, layout = full.pack_general_scalar_state(initial, registry)
    if kind != 'none':
        y0 = np.concatenate((y0, [0., 0.]))
    betas = full.beta_inventory(full.read(spec['beta_file']), model=model,
            stage='EFT1', representation=eft_base.representation,
            layout=layout, registry=registry, kind=kind)
    if mode == 'held_zero':
        betas = zero_extra_betas(betas, spec['new'], spec['reused'],
                                full.MODEL_MASS[kind])
    sol = full.integrate_segment(y0, mu_start=full.s34.START,
            mu_end=full.s34.STOP, layout=layout, registry=registry,
            stage='EFT1', kind=kind, betas=betas, points=points)
    final = full.unpack_general_scalar_state(
        sol.y[:layout.size, -1], layout, registry, uv=False,
        mu_gev=full.s34.STOP)
    mu2_final = (complex(float(sol.y[layout.size, -1]),
                         float(sol.y[layout.size+1, -1]))
                 if kind != 'none' else 0j)
    extra_end = extra_endpoint_values(final, spec['new'], spec['reused'])
    max_extra = max(extra_end.values(), default=0.0)
    if mode == 'held_zero':
        require(max_extra < 1e-13 and abs(mu2_final) < 1e-13,
                'extra scalar interactions did not remain zero')
    spectrum = full.quadratic_mass_hessian(
        full.load_manifest(full.MANIFEST), model,
        final.baseline.mS1Sq, final.baseline.mS2Sq, mu2_final)
    require(spectrum['analytic_spectrum_verified'] is True,
            'scalar eigenmass check failed')

    folder.mkdir(parents=True, exist_ok=False)
    side_path = folder / (full.model_tag(model) + '_eft1_quadratic_running.json')
    sidecar = {
        'format': full.qeft.FORMAT,
        'status': 'CandidateQuadraticAndFullQuarticEFT1Integrated',
        'model': model, 'quadratic_kind': kind,
        'mu_initial_gev': float(full.s34.START),
        'mu_final_gev': float(full.s34.STOP),
        'mu2_beta_included': kind != 'none',
        'mu2_beta_restricted_to_zero': mode == 'held_zero' and kind != 'none',
        'mu2_initial_gev2': full.cjson(0j),
        'mu2_final_gev2': full.cjson(mu2_final),
        'initial_baseline': full.s33.json_state(initial.baseline),
        'final_baseline': full.s33.json_state(final.baseline),
        'initial_new_quartics': {n: full.cjson(initial.extra[n]) for n in spec['new']},
        'final_new_quartics': {n: full.cjson(final.extra[n]) for n in spec['new']},
        'physical_matchete_rgbeta_mapping_certified': False,
        'scalar_mass_threshold_matching_performed': False,
        'wilson_coefficient_matching_performed': False,
        'physical_neutrino_mass_prediction': False,
        'continuous_uv_quadratic_running_performed': False,
        'continuous_uv_fermion_threshold_eft1': True,
        'input_sha256': spec['hashes'], 'settings': cfg,
    }
    full.write_json_atomic(side_path, sidecar)
    rows = full.trajectory_rows(sol, layout=layout, registry=registry,
                                stage='EFT1', kind=kind)
    payload = {
        'format': full.FORMAT, 'status': full.STATUS,
        'model': model, 'settings': cfg, 'input_sha256': spec['hashes'],
        'boundary_mode': 'FROZEN_UV', 'yukawa_mode': 'optimised',
        'yukawa_seed_status': 'historical_frozen_fit_not_full_general_refit',
        'additional_quartic_boundary_status': 'ALL REGISTRY EXTRAS ZERO AT muF',
        'quadratic_kind': kind,
        'quadratic_mu2_uv_gev2': full.cjson(0j),
        'quadratic_mu2_at_fermion_threshold_gev2': full.cjson(0j),
        'quadratic_mu2_at_scalar_endpoint_gev2': full.cjson(mu2_final),
        'new_extra_quartics_evolved': len(spec['new']),
        'reused_extra_quartics_evolved': len(spec['reused']),
        'extra_quadratic_terms_evolved': int(kind != 'none'),
        'uv_ode_evaluations': 0, 'eft1_ode_evaluations': int(sol.nfev),
        'heavy_fermion_threshold_masses_gev':
            full.fermion_threshold_masses(uv).masses_gev.tolist(),
        'uv_trajectory': frozen.frozen_rows(uv, original.extra, mu2, points),
        'eft1_trajectory': rows,
        'uv_fermion_threshold_state': full.s33.json_state(uv),
        'eft1_scalar_endpoint_state': full.s33.json_state(final.baseline),
        'scalar_endpoint_physically_admissible': not spectrum['tachyonic'],
        'scalar_endpoint_spectrum': {k: v for k, v in spectrum.items()
                                     if k != 'real_scalar_hessian_gev2'},
        'saved_candidate_eft1_endpoint': side_path.name,
        'mass_bilinear_tree_threshold_prescription':
            'zero at 100 TeV; candidate EFT1 running or held-zero RG projection',
        'physical_matchete_rgbeta_mapping_certified': False,
        'fermion_threshold_wilson_matching_performed': False,
        'scalar_threshold_wilson_matching_performed': False,
        'full_sequential_neutrino_mass_prediction': False,
        'zero_extra_control_mode': mode,
        'extra_quartic_initial_magnitudes': {name: 0.0 for name in spec['new']+spec['reused']},
        'extra_quartic_final_magnitudes': extra_end,
        'max_extra_quartic_at_scalar_endpoint': max_extra,
        'limitations': [
            'Original historical fitted Yukawas are unchanged and not refitted',
            'Original lambdaT3 and diagonal scalar masses are retained',
            'Frozen UV segment has no UV ODE integration',
            'Full sequential Wilson threshold matching not implemented',
            'Matchete/RGBeta general scalar action normalization not certified',
            ('Extra scalar beta coordinates are forcibly held at zero; '
             'not a full general scalar RG solution') if mode == 'held_zero' else
            'Extra couplings initially zero may be regenerated by candidate beta functions',
            'For self-bilinear sectors, any regenerated nonzero bilinear is OMITTED from C5',
        ],
    }
    full.write_json_atomic(target, payload)
    return payload


def historic_result(benchmark: str, model: str, registry: dict) -> dict:
    case = policy.BENCHMARKS[benchmark][1]
    src = historical.historical(case, model, registry['models'][model])
    raw = src['fit']
    prediction = raw.get('prediction')
    require(isinstance(prediction, dict) and all(f in prediction for f in FIELDS),
            'historical fit prediction absent for ' + model)
    return {
        'fit_status': raw.get('status'),
        'prediction': {f: float(prediction[f]) for f in FIELDS},
    }


def legacy_kernel_control(model: str, payload: dict) -> dict:
    """Separate original-kernel control, SAME held-zero endpoint (no self mixing)."""
    uv = neutrino.decode_uv(payload['uv_fermion_threshold_state'])
    eft = neutrino.decode_eft(payload['eft1_scalar_endpoint_state'])
    c5, source, kind = neutrino.matched_nonmixed_c5(model, uv, eft)
    result = neutrino.low_energy(model, c5, eft)
    return {
        'status': 'OriginalKernelAtHeldZeroEndpoint',
        'matching_kernel_kind': kind,
        'matching_source': policy.relative(ROOT, source),
        'neutrino': result,
    }


def format_row(benchmark: str, model: str, mode: str,
               payload: dict, result: dict, historical: dict,
               legacy: dict | None) -> dict:
    current = result['oscillation_observables_normal_ordering']
    orig = historical['prediction']
    uv_state = payload['uv_fermion_threshold_state']
    eft_state = payload['eft1_scalar_endpoint_state']
    row = {
        'benchmark': benchmark, 'model': model, 'mode': mode,
        'status': result['status'], 'original_fit_status': historical['fit_status'],
        'quadratic_kind': payload['quadratic_kind'],
        'max_extra_quartic_endpoint': payload['max_extra_quartic_at_scalar_endpoint'],
        'abs_mu2_endpoint_gev2': abs(complex(
            payload['quadratic_mu2_at_scalar_endpoint_gev2']['real'],
            payload['quadratic_mu2_at_scalar_endpoint_gev2']['imag'])),
        'lambdaT3_input_abs': abs(complex(uv_state['lambdaT3']['real'],
                                          uv_state['lambdaT3']['imag'])),
        'lambdaT3_endpoint_abs': abs(complex(eft_state['lambdaT3']['real'],
                                             eft_state['lambdaT3']['imag'])),
        'largest_neutrino_mass_ev': result['takagi_masses_ev_descending'][0],
        **{f: float(current[f]) for f in FIELDS},
        **{'historical_' + f: orig[f] for f in FIELDS},
        **{'difference_' + f: float(current[f]) - orig[f] for f in FIELDS},
    }
    if legacy:
        for f in FIELDS:
            row['original_kernel_' + f] = float(
                legacy['neutrino']['oscillation_observables_normal_ordering'][f])
    return row


def write_report(root: Path, settings: dict, rows: list[dict],
                 failures: list[dict]) -> None:
    full.write_json_atomic(root/'summary.json', {
        'format': FORMAT, 'settings': settings,
        'rows': rows, 'failures': failures,
        'full_sequential_matching_certified': False,
        'held_zero_is_full_general_rg_solution': False,
    })
    columns = (['benchmark','model','mode','status','original_fit_status',
                'quadratic_kind','max_extra_quartic_endpoint','abs_mu2_endpoint_gev2',
                'lambdaT3_input_abs','lambdaT3_endpoint_abs',
                'largest_neutrino_mass_ev'] + list(FIELDS) +
               ['historical_' + f for f in FIELDS] +
               ['difference_' + f for f in FIELDS] +
               ['original_kernel_' + f for f in FIELDS])
    tmp = root/'summary.csv.tmp'
    with tmp.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    tmp.replace(root/'summary.csv')


def run(args: argparse.Namespace) -> int:
    require(bool(policy.RUN_ID.fullmatch(args.run_id)), 'invalid --run-id')
    require(4 <= args.points <= 200, 'points outside 4..200')
    models = policy.MODELS if args.models == 'all' else select_csv(args.models, policy.MODELS)
    modes = select_csv(args.modes, MODES)
    root = OUTPUT/args.run_id
    stamp = sha256(Path(__file__).read_bytes()).hexdigest()
    settings = {
        'benchmark': args.benchmark, 'models': list(models), 'modes': list(modes),
        'points': args.points, 'script_sha256': stamp,
        'scope': 'frozen-UV candidate EFT1; zero all registry extras and quadratic terms',
    }
    control_file = root/'control_settings.json'
    if not args.plan:
        if control_file.exists():
            require(args.resume and full.read(control_file) == settings,
                    'settings/source changed; choose new --run-id')
        else:
            require(not root.exists(), 'run directory already exists without settings')
            root.mkdir(parents=True)
            full.write_json_atomic(control_file, settings)

    rows: list[dict] = []
    failures: list[dict] = []
    for model in models:
        try:
            spec = prepare(args.benchmark, model, args.stage22, args.quadratic_betas)
            reference = historic_result(args.benchmark, model, spec['registry'])
        except (OSError, ValueError, KeyError, TypeError, RuntimeError) as exc:
            if args.plan:
                print('BLOCKED PREFLIGHT', model, type(exc).__name__, str(exc), flush=True)
                failures.append({'model': model, 'error':str(exc)})
                continue
            failures.append({'model':model,'mode':'preflight',
                             'error':f'{type(exc).__name__}: {exc}'})
            print('BLOCKED PREFLIGHT',model,exc,flush=True)
            continue
        if args.plan:
            print(f'{model}: zero {len(spec["new"])} new + {len(spec["reused"])} reused '
                  f'extra quartics; quadratic={spec["kind"]}; fit={reference["fit_status"]}',flush=True)
            continue
        for mode in modes:
            source_root = root/args.benchmark/mode
            folder = source_root/full.model_tag(model)
            try:
                data = calculate_scalar(spec, benchmark=args.benchmark, model=model,
                        mode=mode, points=args.points, folder=folder, resume=args.resume)
                require(data['scalar_endpoint_physically_admissible'],
                        'tachyonic scalar endpoint; no neutrino prediction')
                candidate = neutrino.compute_one(model, source_root=source_root,
                        destination_root=source_root/'neutrino_v1', resume=args.resume)
                legacy = None
                if mode == 'held_zero':
                    legacy_file = source_root/'original_kernel_v1'/(
                        full.model_tag(model)+'_original_kernel.json')
                    if legacy_file.exists():
                        require(args.resume, 'original kernel output exists without --resume')
                        legacy = full.read(legacy_file)
                    else:
                        try:
                            legacy = legacy_kernel_control(model, data)
                            full.write_json_atomic(legacy_file, legacy)
                        except (OSError, ValueError, KeyError, TypeError, RuntimeError,
                                ArithmeticError, np.linalg.LinAlgError) as exc:
                            print('ORIGINAL KERNEL UNAVAILABLE',model,exc,flush=True)
                row = format_row(args.benchmark, model, mode, data,
                                 candidate, reference, legacy)
                rows.append(row)
                print(f'{model} {mode}: lambdaT3={row["lambdaT3_endpoint_abs"]:.6g}, '
                      f'sin2theta13={row["sin2_theta13"]:.7g}, '
                      f'max extra={row["max_extra_quartic_endpoint"]:.5g}',flush=True)
            except (OSError, ValueError, KeyError, TypeError, RuntimeError,
                    ArithmeticError, np.linalg.LinAlgError) as exc:
                failures.append({'model': model,'mode':mode,
                                 'error':f'{type(exc).__name__}: {exc}'[:600]})
                print('BLOCKED',model,mode,exc,flush=True)
        write_report(root, settings, rows, failures)
    if args.plan:
        return 2 if failures else 0
    if not rows:
        write_report(root, settings, rows, failures)
    print('Summary:',root/'summary.csv')
    print('Successful:',len(rows),'blocked:',len(failures))
    return 2 if failures else 0


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-id', default='zero_extra_all16_v1')
    p.add_argument('--benchmark', choices=('optimised_low','optimised_high'),
                   default='optimised_low')
    p.add_argument('--models', default='all',
                   help='Comma-separated exact model IDs, or all 16')
    p.add_argument('--modes', default='held_zero,boundary_zero',
                   help='held_zero,boundary_zero or just one')
    p.add_argument('--points', type=int, default=16)
    p.add_argument('--stage22', type=Path, default=full.STAGE22)
    p.add_argument('--quadratic-betas', type=Path, default=full.BETA_DIR)
    p.add_argument('--plan', action='store_true',help='Read-only input validation')
    p.add_argument('--resume', action='store_true',help='Reuse identical saved runs')
    return run(p.parse_args(argv))


if __name__ == '__main__':
    raise SystemExit(main())
