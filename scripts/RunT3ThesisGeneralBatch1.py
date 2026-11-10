#!/usr/bin/env python3
"""A(-4)/D(-2) thesis-benchmark *full inventory* renormalisable RG stage.

Runs the repository's existing one-loop RGBeta UV and two-scalar EFT1 solvers;
does not run Wilson RGEs or either one-loop matching stage. The exact RG
endpoints, heavy-fermion parameters and all live scalar couplings are exported
for the subsequent RunSequentialT3Matching-based matching adapter.

PowerShell, from project root:
  python -X utf8 -m scripts.RunT3ThesisGeneralBatch1 --all \
    --symbolic-run-id thesis_full_general_sym_v1 --run-id thesis_batch1_v1 \
    --generate-symbolic

The symbolic run is cached by the repository's existing canonical tool, not
regenerated for every numerical Yukawa point.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
from datetime import datetime, timezone
from hashlib import sha256
import json
import math
from pathlib import Path
import re
import subprocess
import sys
from typing import Any, Mapping

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Numerical.parallel_t3.background.general_batch1_api import (
    FORMAT, NativeBatch1Trajectory, decode_native, encode_native, read_json,
)

MODELS = ('A_-4', 'D_-2')
BENCHMARK = ROOT / 'configs/thesis_numerics/t3_thesis_A4_Dm2_pilot.json'
REGISTRY = ROOT / 'RGE/general/generated/general_rgbeta_quartic_registry.json'
OUTPUT_ROOT = ROOT / 'output/parallel_t3/full_general_batch1'
HEAVY = frozenset(('y1', 'y2', 'y1Conj', 'y2Conj', 'MF'))
SM = frozenset(('gY', 'g2', 'g3', 'yu', 'yd', 'ye', 'lambdaH'))
CONVENTION = '16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g'


def require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError('Batch1 full-general RG: ' + message)


def hashed(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def dump_atomic(path: Path, data: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.partial')
    tmp.write_text(json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + '\n',
                   encoding='utf-8')
    tmp.replace(path)


def safe_id(value: str) -> str:
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}', value):
        raise argparse.ArgumentTypeError('ID must match [A-Za-z0-9][A-Za-z0-9_-]{0,63}')
    return value


def read_benchmark(path: Path) -> dict:
    cfg = read_json(path)
    require(cfg.get('schema') == 't3-thesis-benchmark-two-model-pilot-v1',
            'expected source thesis benchmark JSON schema')
    require(cfg.get('scales_GeV') == {'UV': 1e7, 'F': 1e5, 'S': 1e3, 'low': 100.},
            'expected thesis 1e7 -> 1e5 -> 1e3 -> 100 GeV scales')
    require(cfg.get('lambdaT3_UV') == 0.1 and cfg.get('lambda_other_UV') == 0.01,
            'thesis scalar coupling inputs unexpectedly changed')
    require(cfg.get('class_E_extra_half_factor') is False,
            'historical class E half factor cannot enter this benchmark')
    require(set(MODELS) <= set(cfg.get('models', {})), 'missing model thesis Yukawa inputs')
    require(len(cfg.get('heavy_fermion_masses_UV_GeV', [])) == 3,
            'need three heavy Dirac masses')
    return cfg


def real_yukawas(cfg: dict, model: str) -> tuple[np.ndarray, np.ndarray]:
    scale = float(cfg['fitting']['input_Yukawa_scale'])
    require(scale == 1e-4, 'unexpected original Yukawa input scaling')
    mats = []
    for key in ('y1_UV_times_1e4', 'y2_UV_times_1e4'):
        raw = np.asarray(cfg['models'][model][key], dtype=float)
        require(raw.shape == (3, 3) and np.isfinite(raw).all(),
                model + ': invalid real Yukawa input ' + key)
        mats.append(raw * scale)
    return mats[0], mats[1]


def encode_matrix(array: np.ndarray) -> list:
    a = np.asarray(array, dtype=np.complex128)
    require(a.shape == (3, 3) and np.isfinite(a).all(), 'expected finite 3x3 matrix')
    return encode_native(a)


def source_inventory(uv: dict, eft: dict, registry: dict, model: str) -> dict:
    """Enforce complete general native RGBeta inventories; reject projections."""
    entry = registry.get('models', {}).get(model)
    require(isinstance(entry, dict), 'model absent from general quartic registry')
    for stage, payload in (('UV', uv), ('EFT1', eft)):
        require(payload.get('format') == 'general-t3-quadratic-rgbeta-experimental-v1'
                and payload.get('status') == 'SuccessCandidateMassBeta',
                stage + ': must supply complete mass-and-quartic RGBeta JSON')
        meta = payload.get('metadata')
        require(isinstance(meta, dict), stage + ': missing metadata')
        for key, value in {'GeneralScalarModel': model,
                           'GeneralScalarStage': stage,
                           'ReportBetaConvention': CONVENTION,
                           'GeneralScalarQuarticBetaBasis': 'RGBetaFullQuarticProjection',
                           'GeneralScalarAdditionalQuartics': entry['total_additional_quartics'],
                           'GeneralScalarQuadraticMassMixing': 'IncludedCandidate',
                           'PhysicalMatcheteRGBetaMassBasisCertified': False}.items():
            require(meta.get(key) == value, stage + ': inconsistent '+key)
        betas = payload.get('report_betas')
        require(isinstance(betas, dict) and betas, stage + ': empty beta inventory')
        require(set(payload.get('beta_symbols', [])) == set(betas),
                stage + ': beta symbol list and expressions differ')
        require(all(isinstance(e, str) and e.strip() for e in betas.values()),
                stage + ': empty beta expression')
        expected_extra = {row['name'] for row in entry['new_rgbeta_symbols']}
        require(expected_extra <= set(betas), stage + ': missing additional quartic beta(s)')
        require(set(meta.get('GeneralScalarNewBetaKeys', [])) == expected_extra,
                stage + ': metadata/registry quartic mismatch')
    uv_names, eft_names = set(uv['report_betas']), set(eft['report_betas'])
    require(eft_names == uv_names - HEAVY,
            'EFT1 must retain every non-fermionic native coupling; missing='
            + repr(sorted((uv_names - HEAVY) - eft_names))
            + '; extra=' + repr(sorted(eft_names - (uv_names - HEAVY))))
    expected_y = {'A_-4': set(), 'D_-2': {'y2Conj'}}[model]
    require((uv_names & {'y1Conj', 'y2Conj'}) == expected_y,
            model + ': incomplete/incorrect conjugate Yukawa inventory')
    expected_new = {'A_-4': 0, 'D_-2': 5}[model]
    expected_mu = {'A_-4': None, 'D_-2': 'muT3GeneralS2Sq'}[model]
    require(entry['total_additional_quartics'] == expected_new,
            model + ': source registry quartic inventory changed')
    for stage, payload in (('UV', uv), ('EFT1', eft)):
        names = set(payload['report_betas'])
        require((expected_mu in names) if expected_mu else
                not any(x.startswith('muT3General') for x in names),
                stage + ': holomorphic scalar bilinear coverage mismatch')
        require(payload['metadata'].get('QuadraticMassSymbol') == expected_mu,
                stage + ': general quadratic symbol mismatch')
    require({'y1', 'y2', 'MF', 'mS1Sq', 'mS2Sq', 'lambdaT3'} <= uv_names,
            'UV missing required mass or topology interactions')
    return {'uv': sorted(uv_names), 'eft1': sorted(eft_names),
            'additional_quartics': sorted(x['name'] for x in entry['new_rgbeta_symbols']),
            'holomorphic_mass': expected_mu,
            'extra_yukawa': sorted(expected_y)}


def thesis_overrides(cfg: dict, model: str, names: list[str], options: Any) -> dict:
    """Reuse Workstream-B SM flavour conventions; replace all thesis BSM inputs."""
    from Numerical.parallel_t3.background.runner import _bench_overrides
    base = _bench_overrides(model, names, options)
    require(set(base) <= set(names), 'development initializer added unknown beta variables')
    values = dict(base)
    y1, y2 = real_yukawas(cfg, model)
    values.update({
        'MF': encode_matrix(np.diag(np.asarray(cfg['heavy_fermion_masses_UV_GeV'], dtype=float))),
        'y1': encode_matrix(y1), 'y2': encode_matrix(y2),
        'lambdaH': float(cfg['lambdaH_UV']),
        'lambdaT3': float(cfg['lambdaT3_UV']),
    })
    values.update({n: float(v) for n, v in cfg['gauge_UV'].items()})
    values.update({n: float(v) for n, v in cfg['scalar_masses_squared_UV_GeV2'].items()})
    for name in names:
        if name.startswith('lambda') and name not in ('lambdaH', 'lambdaT3'):
            values[name] = float(cfg['lambda_other_UV'])
        elif name.startswith('muT3General'):
            values[name] = {'real': 0.0, 'imag': 0.0}
        elif name in ('y1Conj', 'y2Conj'):
            values[name] = encode_matrix(np.zeros((3, 3), dtype=complex))
    require(set(values) <= set(names), 'UV overrides contain unknown couplings')
    return values


def _shifted_boundary(boundary: Any, file: Path | None, model: str, mu_f: float):
    """Optional *explicit* EFT1 renormalisable corrections from a future F match.

    JSON form: {"format":"t3-batch1-eft1-shifts-v1", "model":"D_-2",
                "mu_F_GeV":100000, "deltas":{"lambdaH":0.001}}
    Does not imply that any one-loop matching was performed here.
    """
    if file is None:
        return boundary, {}, None
    from Numerical.core.GeneralCouplingLayout import spec_for
    from Numerical.running.IntermediateScalarRGBetaEvaluator import intermediate_scalar_state_environment
    content = read_json(file)
    require(content.get('format') == 't3-batch1-eft1-shifts-v1'
            and content.get('model') == model
            and content.get('mu_F_GeV') == mu_f,
            'fermion boundary shift manifest differs from the model/scale')
    deltas = content.get('deltas')
    require(isinstance(deltas, dict), 'boundary deltas must be an object')
    native = intermediate_scalar_state_environment(boundary)
    require(set(deltas) <= set(native) and not (set(deltas) & HEAVY),
            'cannot shift unknown or integrated-out parameters')
    updated = dict(native)
    parsed = {}
    for key, json_value in deltas.items():
        delta = decode_native(json_value)
        spec = spec_for(key)
        delta = spec.validate(delta, majorana=False)
        target = spec.validate(np.asarray(native[key]) + delta, majorana=False)
        updated[key] = target
        parsed[key] = delta
    # Preserve the tested repository data classes rather than bypassing them.
    new_sm = replace(boundary.sm, **{k: updated[k] for k in deltas if k in SM})
    general = dict(boundary.general_couplings)
    for key in deltas:
        if key in general:
            general[key] = updated[key]
    ordinary = {k: updated[k] for k in deltas
                if k not in SM and k not in general}
    shifted = replace(boundary, sm=new_sm, general_couplings=general,
                      **ordinary).validated()
    check = intermediate_scalar_state_environment(shifted)
    for key in updated:
        require(np.allclose(np.asarray(check[key]), np.asarray(updated[key]),
                            rtol=1e-12, atol=1e-12),
                'boundary update lost scalar parameter ' + key)
    return shifted, parsed, hashed(file)


def snapshots(run: Any, stage: str, names: list[str]) -> dict[str, np.ndarray]:
    from Numerical.core.RGBetaEvaluator import state_environment
    from Numerical.running.IntermediateScalarRGBetaEvaluator import intermediate_scalar_state_environment
    env_fn = state_environment if stage == 'uv' else intermediate_scalar_state_environment
    values: dict[str, list[np.ndarray]] = {n: [] for n in names}
    for i in range(len(run.mu_gev)):
        environment = env_fn(run.state_at_index(i))
        require(set(environment) == set(names),
                f'{stage}: source state keys differ from full beta inventory; '
                f'missing={sorted(set(names)-set(environment))}, '
                f'extra={sorted(set(environment)-set(names))}')
        for name in names:
            arr = np.asarray(environment[name], dtype=np.complex128)
            require(np.isfinite(arr).all(), f'{stage}: nonfinite {name}')
            values[name].append(arr)
    return {name: np.stack(series).astype(np.complex128)
            for name, series in values.items()}


def _git_head() -> str | None:
    try:
        result = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=ROOT,
                                check=True, capture_output=True, timeout=5, text=True)
        return result.stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return None


def run_one(model: str, cfg_path: Path, symbolic_run_id: str, run_id: str,
            *, generate_symbolic: bool, samples_uv: int, samples_eft1: int,
            rtol: float, atol: float, shifts: Path | None) -> dict:
    from Numerical.orchestration.GeneralRunningBenchmark import (
        initial_state, load_canonical_general_exports,
    )
    from Numerical.parallel_t3.background.runner import Options
    from Numerical.running.FermionThresholdBoundary import (
        fermion_threshold_masses, build_intermediate_scalar_boundary,
    )
    from Numerical.running.UVRunner import run_uv_segment
    from Numerical.running.IntermediateScalarRunner import run_intermediate_scalar_segment
    from Numerical.running.IntermediateScalarRGBetaEvaluator import intermediate_scalar_state_environment
    from Numerical.core.RGBetaEvaluator import state_environment

    cfg = read_benchmark(cfg_path)
    scales = cfg['scales_GeV']
    uv_payload, eft_payload, provenance = load_canonical_general_exports(
        model, symbolic_run_id, generate=generate_symbolic,
    )
    registry = read_json(REGISTRY)
    inventory = source_inventory(uv_payload, eft_payload, registry, model)
    opts = Options(root=ROOT, mu_uv=scales['UV'], mu_f=scales['F'], mu_s=scales['S'])
    uv_overrides = thesis_overrides(cfg, model, inventory['uv'], opts)
    source_model, state = initial_state(
        uv_payload, mu_uv=scales['UV'], mu_f=scales['F'], mu_s=scales['S'],
        overrides=uv_overrides,
    )
    require(source_model == model, 'initial state model mismatch')
    require(set(state_environment(state)) == set(inventory['uv']),
            'initial UV state is not closed under complete general RGEs')
    uv = run_uv_segment(
        state, uv_payload, scales['F'],
        save_scales_gev=np.geomspace(scales['UV'], scales['F'], samples_uv),
        rtol=rtol, atol=atol,
    )
    at_f_uv = uv.final_state
    fermions = fermion_threshold_masses(at_f_uv)
    eft_raw = build_intermediate_scalar_boundary(
        at_f_uv, matching_scale_gev=scales['F'],
    )
    require(set(intermediate_scalar_state_environment(eft_raw)) == set(inventory['eft1']),
            'complete EFT1 native coupling inventory not preserved at F threshold')
    eft_boundary, deltas, shifts_sha = _shifted_boundary(eft_raw, shifts, model, scales['F'])
    eft = run_intermediate_scalar_segment(
        eft_boundary, eft_payload, scales['S'],
        save_scales_gev=np.geomspace(scales['F'], scales['S'], samples_eft1),
        rtol=rtol, atol=atol,
    )
    uv_snap = snapshots(uv, 'uv', inventory['uv'])
    eft_snap = snapshots(eft, 'eft1', inventory['eft1'])
    for name in inventory['eft1']:
        expected = uv_snap[name][-1] + np.asarray(deltas.get(name, 0), dtype=complex)
        require(np.allclose(expected, eft_snap[name][0], rtol=1e-10, atol=1e-9),
                'fermion-boundary state discontinuity (not accounted for by explicit shifts): ' + name)
    require(set(eft_snap) == set(uv_snap)-HEAVY, 'heavy fermion coordinates survived in EFT1')
    dest = OUTPUT_ROOT / run_id / model
    require(not dest.exists(), 'existing output (choose a fresh --run-id): ' + str(dest))
    dest.mkdir(parents=True)
    npz = dest / 'trajectories.npz'
    arrays: dict[str, np.ndarray] = {
        'uv__mu_gev': np.asarray(uv.mu_gev, dtype=np.float64),
        'eft1__mu_gev': np.asarray(eft.mu_gev, dtype=np.float64),
    }
    arrays.update({'uv__'+key: value for key, value in uv_snap.items()})
    arrays.update({'eft1__'+key: value for key, value in eft_snap.items()})
    np.savez_compressed(npz, **arrays)
    uv_first = state_environment(uv.initial_state)
    f_uv = state_environment(at_f_uv)
    f_eft_raw = intermediate_scalar_state_environment(eft_raw)
    f_eft = intermediate_scalar_state_environment(eft_boundary)
    s_eft = intermediate_scalar_state_environment(eft.final_state)
    f_input = {
        'format': FORMAT + '-fermion-matching-input',
        'model': model, 'scale_GeV': scales['F'],
        'all_heavy_fermions_active_above': True,
        'all_heavy_fermions_to_integrate_out_together': True,
        'both_heavy_scalar_multiplets_remain_light': True,
        'native_uv_rgbeta_couplings_at_muF': encode_native(f_uv),
        'native_eft1_tree_projection_at_muF': encode_native(f_eft_raw),
        'native_eft1_boundary_used_at_muF': encode_native(f_eft),
        'explicit_eft1_boundary_shifts': encode_native(deltas),
        'heavy_fermion_singular_masses_GeV': encode_native(fermions.masses_gev),
        'singular_mass_note': 'singular values only; no physical Matchete fermion mass-basis rotation applied',
        'native_coordinate_basis': 'RGBeta interaction and flavour coordinates, not Matchete tensor coefficients',
        'symbolic_matching_engine_to_reuse': 'Lagrangian/RunMatching.wl::RunSequentialT3Matching',
        'fermion_one_loop_wilson_matching_performed': False,
        'renormalisable_one_loop_boundary_shifts_derived_by_batch1': False,
    }
    from Numerical.experimental.GeneralT3ScalarQuadraticMasses import (
        load_manifest, quadratic_mass_hessian,
    )
    quadratic_source = load_manifest(ROOT / 'docs/general_scalar_trial_operator_manifest.json')
    hessian = quadratic_mass_hessian(
        quadratic_source, model, float(s_eft['mS1Sq']), float(s_eft['mS2Sq']),
        complex(s_eft[inventory['holomorphic_mass']])
        if inventory['holomorphic_mass'] else 0j,
    )
    s_input = {
        'format': FORMAT + '-scalar-matching-input',
        'model': model, 'scale_GeV': scales['S'],
        'both_heavy_scalar_multiplets_to_integrate_out_together': True,
        'no_heavy_fermions_active': True,
        'native_eft1_rgbeta_couplings_at_muS': encode_native(s_eft),
        'quadratic_holomorphic_mass_symbol': inventory['holomorphic_mass'],
        'scalar_quadratic_terms_GeV2': encode_native({
            name: s_eft[name] for name in ('mS1Sq', 'mS2Sq')
            + ((inventory['holomorphic_mass'],) if inventory['holomorphic_mass'] else ())
        }),
        'general_quartics': encode_native({name: s_eft[name]
            for name in inventory['additional_quartics']}),
        'scalar_mass_diagonalisation_required_before_matchete_projection':
            bool(inventory['holomorphic_mass']),
        'candidate_spherical_scalar_quadratic_spectrum': hessian,
        'mass_basis_mapping_performed': False,
        'fermion_inherited_wilson_coefficients_included': False,
        'scalar_one_loop_wilson_matching_performed': False,
    }
    f_path = dest / 'fermion_matching_input.json'
    s_path = dest / 'scalar_matching_input.json'
    dump_atomic(f_path, f_input)
    dump_atomic(s_path, s_input)
    provenance_paths = {
        'thesis_config': cfg_path,
        'general_registry': REGISTRY,
        'runner': Path(__file__),
        'api': ROOT/'Numerical/parallel_t3/background/general_batch1_api.py',
    }
    source_hashes = {k: hashed(v) for k, v in provenance_paths.items()}
    for stage, path in (('uv', Path(provenance['uv_json'])),
                        ('eft1', Path(provenance['eft1_json']))):
        # Symbolic provenance can return relative paths on some installations.
        actual = path if path.is_absolute() else ROOT / path
        source_hashes[stage+'_rgbeta'] = hashed(actual)
    meta = {
        'format': FORMAT, 'status': 'CandidateGeneralRenormalisableRGIntegrated',
        'model': model, 'run_id': run_id,
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'local_git_head_read_only': _git_head(),
        'scales_gev': scales,
        'coupling_order': {'uv': inventory['uv'], 'eft1': inventory['eft1']},
        'coupling_shapes': {
            'uv': {n: list(uv_snap[n].shape[1:]) for n in uv_snap},
            'eft1': {n: list(eft_snap[n].shape[1:]) for n in eft_snap},
        },
        'native_parameter_convention': CONVENTION,
        'symbolic_run_id': symbolic_run_id,
        'symbolic_provenance': provenance,
        'source_sha256': source_hashes,
        'explicit_fermion_boundary_shift_source_sha256': shifts_sha,
        'eft1_boundary_deltas': encode_native(deltas),
        'inventory': inventory,
        'thesis_parameters_at_uv': encode_native(uv_first),
        'numerical_endpoints': {
            'UV_1e7': encode_native(uv_first),
            'UV_at_F_1e5': encode_native(f_uv),
            'EFT1_at_F_1e5': encode_native(f_eft),
            'EFT1_at_S_1e3': encode_native(s_eft),
        },
        'solvers': {'uv': {'success': bool(uv.solver_success), 'nfev': int(uv.nfev)},
                    'eft1': {'success': bool(eft.solver_success), 'nfev': int(eft.nfev)}},
        'trajectory_sha256': hashed(npz),
        'fermion_matching_input_sha256': hashed(f_path),
        'scalar_matching_input_sha256': hashed(s_path),
        'fermion_threshold_wilson_matching_performed': False,
        'scalar_threshold_wilson_matching_performed': False,
        'EFT1_wilson_running_performed': False,
        'mass_basis_mapping_certified': False,
        'neutrino_observables_computed': False,
        'heavy_fermion_threshold_mode': 'all_F_simultaneously_at_1e5_GeV',
        'scalar_threshold_mode': 'both_S_simultaneously_at_1e3_GeV',
        'renormalisable_projection_at_F': ('tree_level_plus_explicit_shifts' if deltas
                                          else 'tree_level_only'),
    }
    dump_atomic(dest / 'manifest.json', meta)
    # Decode and verify the actual saved/hashed interface before reporting success.
    NativeBatch1Trajectory.load(dest)
    return {'model': model, 'status': meta['status'],
            'uv_nfev': int(uv.nfev), 'eft1_nfev': int(eft.nfev),
            'output_directory': str(dest.relative_to(ROOT)),
            'state_variables_uv': len(uv_snap),
            'state_variables_eft1': len(eft_snap),
            'D2_y2Conj_RGE_evolved': 'y2Conj' in uv_snap,
            'D2_quadratic_mu2_RGE_evolved': inventory['holomorphic_mass'] in eft_snap,
            'D2_general_quartics_RGE_evolved': len(inventory['additional_quartics'])}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    choice = parser.add_mutually_exclusive_group(required=True)
    choice.add_argument('--model', choices=MODELS)
    choice.add_argument('--all', action='store_true')
    parser.add_argument('--run-id', type=safe_id, required=True)
    parser.add_argument('--symbolic-run-id', type=safe_id, required=True)
    parser.add_argument('--generate-symbolic', action='store_true',
                        help='Run/reuse canonical Wolfram general exports; not for every Yukawa evaluation')
    parser.add_argument('--benchmark', type=Path, default=BENCHMARK)
    parser.add_argument('--saved-uv-points', type=int, default=17)
    parser.add_argument('--saved-eft1-points', type=int, default=25)
    parser.add_argument('--rtol', type=float, default=1e-8)
    parser.add_argument('--atol', type=float, default=1e-11)
    parser.add_argument('--eft1-boundary-shifts', type=Path,
                        help='Optional Batch 2 native renormalisable matching deltas; do not infer them')
    args = parser.parse_args(argv)
    require(min(args.saved_uv_points, args.saved_eft1_points) >= 3,
            'need three or more saved points per RGE segment')
    require(0 < args.rtol < 1e-3 and 0 < args.atol < 1e-5,
            'positive controlled tolerances required')
    models = MODELS if args.all else (args.model,)
    results = {}
    for model in models:
        try:
            print('Running full-inventory renormalisable RG for', model, flush=True)
            row = run_one(model, args.benchmark.resolve(), args.symbolic_run_id, args.run_id,
                          generate_symbolic=args.generate_symbolic,
                          samples_uv=args.saved_uv_points,
                          samples_eft1=args.saved_eft1_points,
                          rtol=args.rtol, atol=args.atol,
                          shifts=args.eft1_boundary_shifts.resolve()
                          if args.eft1_boundary_shifts else None)
            results[model] = row
            print(model, row['status'], row['output_directory'], flush=True)
        except Exception as exc:
            results[model] = {'model': model, 'status': 'Blocked_'+type(exc).__name__,
                              'error': str(exc)}
            print(model, results[model]['status'], str(exc), file=sys.stderr, flush=True)
    out = OUTPUT_ROOT / args.run_id
    out.mkdir(parents=True, exist_ok=True)
    dump_atomic(out/'summary.json', {
        'format': FORMAT+'-summary', 'models': results,
        'successful': sum(row['status']=='CandidateGeneralRenormalisableRGIntegrated'
                          for row in results.values()),
        'total': len(models),
    })
    return 0 if all(row['status']=='CandidateGeneralRenormalisableRGIntegrated'
                    for row in results.values()) else 1


if __name__ == '__main__':
    raise SystemExit(main())
