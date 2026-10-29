"""Stage 36: use model-specific existing C5 matching artifacts and run SM+Weinberg.

Consumes the checked Stage 33 heavy-fermion state and Stage 35 scalar/SM
boundary, with no rerun of the completed UV or intermediate integrations.

Default: only validated hierarchical fixed-order coefficient artifacts.
Extreme-hypercharge models have only legacy one-loop C5 expressions in the
current repository. Their use requires --allow-legacy and is clearly reported.
Neither option certifies the new Stage29 quartic mapping or matches those
additional quartics. No tuning/refitting is performed.

  python -m scripts.FinalWeinbergRunning --model A_+0
  python -m scripts.FinalWeinbergRunning --all --allow-legacy
"""
from __future__ import annotations

import argparse
from dataclasses import fields
import json
import math
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np

from scripts import ScalarThresholdHandoff as s35

ROOT = s35.ROOT
DATA = s35.DATA
MODELS = s35.MODELS
MU_F = s35.START
MU_S = s35.STOP
MU_LOW = 100.0
VEV = 246.22
FORMAT = 'stage36-t3-existing-c5-sm-weinberg-v1'
SUMMARY_FORMAT = 'stage36-16-model-existing-c5-summary-v1'
STATUS = 'SuccessExistingC5WeinbergRunning'
SUMMARY_STATUS = 'Success16ExistingC5WeinbergRunsMixedMatchingSources'
LEGACY_MODELS = frozenset(('A_-4', 'B_-3', 'C_-3', 'E_-4'))


def fail(message: str) -> None:
    raise ValueError('Stage36 fail-closed: ' + message)


def source_path(model: str, *, root: Path = ROOT,
                allow_legacy: bool = False) -> tuple[Path, str]:
    """Never borrow a matching coefficient from another model/alpha."""
    if model not in MODELS:
        fail('Unsupported model: ' + str(model))
    group, alpha = model.split('_')
    sign = 'p' if alpha.startswith('+') else 'm'
    token = sign + alpha[1:]
    dirname = f'T3_{group}_alpha_{token}'
    if model in LEGACY_MODELS:
        if not allow_legacy:
            fail(model + ': only legacy C5 is available; rerun with --allow-legacy '
                 'to enable an explicitly nonuniform rough matching approximation')
        path = Path(root) / 'output' / 'hypercharge' / dirname / 'data' / 'c5_coefficient.txt'
        kind = 'legacy-one-loop-c5-no-hierarchical-hard-plus-direct'
    else:
        path = (Path(root) / 'output' / 'full' / 'hypercharge' / dirname
                / 'data' / 'final_weinberg_coefficient.json')
        kind = 'hierarchical-fixed-order-hard-plus-direct'
    if not path.is_file():
        fail(f'{model}: required model-specific matching artifact missing: {path}')
    return path, kind


def verify_hierarchical_artifact(path: Path) -> None:
    record = s35.read_json(path)
    if (record.get('status') != 'Success'
            or record.get('scheme') != 'fixed_order_one_loop_MSbar_hierarchical'
            or record.get('renormalized_hard_threshold', {}).get('status') != 'Success'
            or record.get('combined', {}).get('ready_for_physical_majorana_numerics') is not True
            or record.get('direct_running_contribution', {}).get('manifestly_symmetric_under_pq') is not True):
        fail('Hierarchical C5 matching artifact not validated/ready: ' + str(path))


def decode_uv_state(raw: Any):
    """Reconstruct the validated ordinary T3 UV state, including full flavor."""
    from Numerical.core.State import T3UVState, T3Representation, SMNumericalState
    from scripts import IntermediateEFTRunning as s34
    if not isinstance(raw, dict) or set(raw) != set(T3UVState.__dataclass_fields__):
        fail('Stage33 saved UV state inventory is incomplete')
    rep = raw.get('representation')
    sm = raw.get('sm')
    if (not isinstance(rep, dict)
            or set(rep) != set(T3Representation.__dataclass_fields__)
            or rep.get('shared_scalar') is not False):
        fail('Stage33 UV representation inventory invalid')
    if not isinstance(sm, dict) or set(sm) != set(SMNumericalState.__dataclass_fields__):
        fail('Stage33 SM inventory invalid')
    rep_obj = T3Representation(**rep).validated()
    sm_obj = SMNumericalState(
        gY=s34._scalar(sm['gY'], 'gY'),
        g2=s34._scalar(sm['g2'], 'g2'),
        g3=s34._scalar(sm['g3'], 'g3'),
        lambdaH=s34._scalar(sm['lambdaH'], 'lambdaH'),
        yu=s34._matrix(sm['yu'], 'yu'),
        yd=s34._matrix(sm['yd'], 'yd'),
        ye=s34._matrix(sm['ye'], 'ye')).validated()
    complex_names = {'lambdaT3'} | set(s34.REUSED_NAMES)
    others: dict[str, Any] = {}
    for f in fields(T3UVState):
        name = f.name
        if name in ('mu_gev', 'representation', 'sm'):
            continue
        value = raw[name]
        others[name] = (None if value is None else
                        s34._matrix(value, name) if name in ('y1', 'y2', 'MF') else
                        s34._complex(value, name) if name in complex_names else
                        s34._scalar(value, name))
    uv = T3UVState(mu_gev=s34._scalar(raw['mu_gev'], 'mu_gev'),
                   representation=rep_obj, sm=sm_obj, **others).validated()
    if not math.isclose(uv.mu_gev, MU_F, rel_tol=1e-12):
        fail('Stage33 UV state is not at the fermion threshold')
    return uv


def checked_inputs(model: str, folder: Path):
    """Check the Stage 33/34/35 hash chain before using any boundary values."""
    folder = Path(folder)
    tag = s35.s34.safe_model(model)
    summary_path = folder / 'stage35_16_model_scalar_threshold_summary.json'
    individual_path = folder / f'stage35_{tag}_scalar_threshold_handoff.json'
    stage33_path = folder / f'stage33_{tag}_fermion_threshold_handoff.json'
    summary = s35.read_json(summary_path)
    if (summary.get('format') != s35.SUMMARY_FORMAT
            or summary.get('status') != s35.SUMMARY_STATUS
            or summary.get('models_checked') != 16
            or summary.get('renormalisable_scalar_threshold_projections') != 16
            or summary.get('scalar_wilson_matching_performed') is not False
            or summary.get('physical_coupling_mapping_certified') is not False
            or set(summary.get('models', {})) != set(MODELS)):
        fail('Stage35 combined summary invalid')
    item = s35.read_json(individual_path)
    if (summary['models'][model].get('report_sha256') != s35.sha(individual_path)
            or item.get('format') != s35.FORMAT
            or item.get('status') != s35.STATUS
            or item.get('model') != model
            or item.get('mu_scalar_threshold_gev') != MU_S
            or item.get('scalar_wilson_matching_performed') is not False
            or item.get('fermion_wilson_matching_performed') is not False
            or item.get('physical_coupling_mapping_certified') is not False
            or item.get('new_quartics_matched_or_evolved_below_mu_f') is not False
            or item.get('source_stage34_summary_sha256') != summary.get('source_stage34_summary_sha256')):
        fail('Stage35 report/status/hash mismatch: ' + model)
    # Reuse Stage35's validated Stage34 proof and source-hash walk.
    upstream, _, upstream_path, upstream_summary = s35._checked_report(model, folder)
    if (item.get('source_stage34_report_sha256') != s35.sha(upstream_path)
            or item.get('source_stage34_summary_sha256') != s35.sha(upstream_summary)):
        fail('Stage35 upstream RGE provenance changed')
    if not stage33_path.is_file():
        fail('Stage33 full UV state unavailable: ' + str(stage33_path))
    stage33 = s35.read_json(stage33_path)
    if (stage33.get('format') != s35.s34.s33.FORMAT
            or stage33.get('status') != 'SuccessRenormalisableFermionThresholdProjection'
            or stage33.get('model') != model
            or stage33.get('mu_fermion_threshold_gev') != MU_F
            or stage33.get('wilson_coefficient_matching_performed') is not False
            or stage33.get('physical_coupling_mapping_certified') is not False):
        fail('Stage33 UV-state provenance invalid')
    # Hash was already checked by Stage35 against Stage34 source inventory.
    hashes = upstream.get('input_sha256', {})
    if s35.sha(stage33_path) not in set(hashes.values()):
        fail('Stage34 inputs do not reference this Stage33 handoff')
    uv = decode_uv_state(stage33.get('uv_threshold_state'))
    intermediate = s35.decode_final_state(item.get('intermediate_final_state_at_mu_s'))
    from Numerical.running.FermionThresholdBoundary import fermion_threshold_masses
    calculated_mf = fermion_threshold_masses(uv).masses_gev
    recorded_mf = np.asarray(stage33.get('heavy_fermion_singular_masses_gev'), dtype=float)
    if recorded_mf.shape != (3,) or not np.allclose(calculated_mf, recorded_mf, rtol=3e-10, atol=1e-7):
        fail('Stage33 heavy mass spectrum changed')
    from Numerical.running.ScalarThresholdBoundary import (
        scalar_threshold_masses, build_final_sm_boundary)
    recorded_s = np.asarray(item.get('scalar_threshold_masses_gev'), dtype=float)
    calculated_s = np.asarray(scalar_threshold_masses(intermediate).masses_gev)
    if recorded_s.shape != (2,) or not np.allclose(recorded_s, calculated_s, rtol=3e-10, atol=1e-7):
        fail('Stage35 scalar mass spectrum changed')
    final_boundary = build_final_sm_boundary(intermediate, matching_scale_gev=MU_S)
    s35._compare_json(s35.s34.s33.json_state(final_boundary),
                      item.get('final_sm_boundary_at_mu_s'), 'Stage35 SM boundary')
    return uv, intermediate, final_boundary, stage33_path, individual_path, summary_path


def matrix_json(array: np.ndarray) -> dict[str, Any]:
    arr = np.asarray(array, dtype=complex)
    if arr.shape != (3, 3) or not np.all(np.isfinite(arr.real)) or not np.all(np.isfinite(arr.imag)):
        fail('Expected finite 3x3 complex matrix')
    return {'real': arr.real.tolist(), 'imag': arr.imag.tolist()}


def _valid_c5(value: Any, label: str) -> np.ndarray:
    matrix = np.asarray(value, dtype=complex)
    if matrix.shape != (3, 3) or not (np.isfinite(matrix.real).all() and np.isfinite(matrix.imag).all()):
        fail(label + ': expected finite 3x3 C5')
    if not np.allclose(matrix, matrix.T, rtol=1e-10, atol=1e-17):
        fail(label + ': C5 must be symmetric')
    if np.linalg.norm(matrix) <= 0:
        fail(label + ': C5 matrix vanishes')
    return matrix


def run_model(model: str, *, folder: Path = DATA, root: Path = ROOT,
              allow_legacy: bool = False, output: Path | None = None) -> dict:
    if model not in MODELS:
        fail('Unknown model ' + str(model))
    uv, eft, sm_boundary, stage33_path, stage35_path, stage35_summary = checked_inputs(model, folder)
    artifact, kind = source_path(model, root=root, allow_legacy=allow_legacy)
    if kind.startswith('hierarchical'):
        verify_hierarchical_artifact(artifact)
    from Numerical.running.FinalC5TrajectoryAdapter import extract_final_c5_inputs
    from Numerical.running.SMWeinbergStage import (
        evaluate_final_weinberg_json, evaluate_symbolic_c5)
    from Numerical.running.ScalarThresholdBoundary import build_sm_weinberg_initial_conditions
    from Numerical.running.SMWeinbergEvolution import evolve_sm_weinberg
    from physics.NeutrinoMass import numerical_neutrino_mass_matrix
    from physics.NeutrinoTrajectory import charged_lepton_mass_basis_matrix
    from physics.NeutrinoObservables import calculate_neutrino_observables

    trajectory_view = SimpleNamespace(uv_threshold_state=uv,
                                      eft1_threshold_state=eft,
                                      mu_scalar_threshold_gev=MU_S)
    matching_input = extract_final_c5_inputs(trajectory_view)
    masses = np.asarray(matching_input.heavy_masses_gev, dtype=float)
    if not np.allclose(np.sort(masses), np.sort(np.linalg.svd(uv.MF, compute_uv=False)),
                       rtol=1e-9, atol=1e-7):
        fail('Mass-basis adapter heavy masses mismatch UV state')
    config = matching_input.as_config()
    # Never substitute fitted/production C5 or use one model's artifact for another.
    c5_at_s = (_valid_c5(evaluate_final_weinberg_json(artifact, config), 'matched C5')
               if kind.startswith('hierarchical') else
               _valid_c5(evaluate_symbolic_c5(artifact, config), 'legacy C5'))
    init = build_sm_weinberg_initial_conditions(sm_boundary, c5_at_s)
    result = evolve_sm_weinberg(init, MU_S, MU_LOW, rtol=1e-8, atol=1e-11)
    if not result.solver_success or result.mu_initial != MU_S or result.mu_final != MU_LOW:
        fail('Final SM+Weinberg solver returned invalid endpoint')
    c5_low = _valid_c5(result.K, 'low-scale C5')
    mass = numerical_neutrino_mass_matrix(c5_low, vev_gev=VEV)
    mass_flavor = charged_lepton_mass_basis_matrix(mass, result.ye)
    obs = calculate_neutrino_observables(mass_flavor, ordering='NO')
    if (not math.isfinite(obs.takagi_residual) or obs.takagi_residual > 1e-7
            or not np.isfinite(obs.masses_ev).all()):
        fail('Low-energy neutrino Takagi/observables invalid')
    record = {
        'format': FORMAT, 'status': STATUS, 'model': model,
        'mu_scalar_threshold_gev': MU_S, 'mu_low_gev': MU_LOW,
        'sm_weinberg_running_performed': True,
        'source_existing_c5_used': True,
        'existing_matching_artifact_kind': kind,
        'matching_source_is_hierarchical': kind.startswith('hierarchical'),
        'legacy_matching_approximation_opted_in': kind.startswith('legacy'),
        'new_general_quartics_included_in_matching': False,
        'c12_wilson_operators_evolved_from_stage33': False,
        'scalar_wilson_matching_recomputed': False,
        'physical_coupling_mapping_certified': False,
        'neutrino_results_are_conditional_benchmark_estimates': True,
        'matching_input_heavy_masses_gev': masses.tolist(),
        'matching_input_scalar_masses_gev': [matching_input.scalar_mass1_gev,
                                             matching_input.scalar_mass2_gev],
        'matching_input_lambdaT3': {'real': matching_input.lambda_t3.real,
                                    'imag': matching_input.lambda_t3.imag},
        'c5_at_scalar_threshold': matrix_json(c5_at_s),
        'c5_at_low_scale': matrix_json(c5_low),
        'neutrino_mass_matrix_low_scale_gev': matrix_json(mass),
        'final_sm_couplings': {'gY': float(result.gY), 'g2': float(result.g2),
                               'g3': float(result.g3), 'lambdaH': float(result.lambdaH)},
        'neutrino_observables': {
            'ordering': 'NO', 'masses_ev': obs.masses_ev.tolist(),
            'delta_m21_sq_ev2': float(obs.delta_m21_sq_ev2),
            'delta_m31_sq_ev2': float(obs.delta_m31_sq_ev2),
            'delta_m32_sq_ev2': float(obs.delta_m32_sq_ev2),
            'pmns_abs': np.asarray(obs.pmns_abs).tolist(),
            'takagi_residual': float(obs.takagi_residual),
        },
        'solver_nfev': int(result.nfev),
        'source_matching_artifact': str(artifact.relative_to(Path(root))) if artifact.is_relative_to(Path(root)) else str(artifact),
        'source_matching_artifact_sha256': s35.sha(artifact),
        'source_stage33_report_sha256': s35.sha(stage33_path),
        'source_stage35_report_sha256': s35.sha(stage35_path),
        'source_stage35_summary_sha256': s35.sha(stage35_summary),
        'limitations': [
            'Existing matching artifact only; not a newly derived C5 from Stage29 quartics',
            'New Stage29 quartics excluded from intermediate EFT and Wilson matching',
            'No new C12 Wilson evolution; existing fixed-order coefficient is reused',
            'Legacy C5 (if opted in) lacks the hierarchical hard-plus-direct matching scheme',
            'Physical Matchete/RGBeta quartic normalization map uncertified',
            'Common threshold scales are benchmark approximations',
        ],
    }
    target = Path(output) if output is not None else Path(folder) / f'stage36_{s35.s34.safe_model(model)}_low_energy_neutrinos.json'
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(record, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    return record


def write_summary(*, folder: Path = DATA, require_all: bool = True) -> dict:
    folder = Path(folder)
    upstream_path = folder / 'stage35_16_model_scalar_threshold_summary.json'
    upstream = s35.read_json(upstream_path)
    if (upstream.get('status') != s35.SUMMARY_STATUS or upstream.get('models_checked') != 16
            or set(upstream.get('models', {})) != set(MODELS)):
        fail('Missing/invalid Stage35 16-model summary')
    rows = {}
    for model in MODELS:
        p = folder / f'stage36_{s35.s34.safe_model(model)}_low_energy_neutrinos.json'
        if not p.is_file():
            if require_all:
                fail('Missing Stage36 result for ' + model)
            continue
        row = s35.read_json(p)
        source_match = Path(ROOT) / row.get('source_matching_artifact', '')
        stage35_single = folder / f'stage35_{s35.s34.safe_model(model)}_scalar_threshold_handoff.json'
        if (row.get('format') != FORMAT or row.get('status') != STATUS
                or row.get('model') != model or row.get('sm_weinberg_running_performed') is not True
                or row.get('physical_coupling_mapping_certified') is not False
                or row.get('new_general_quartics_included_in_matching') is not False
                or row.get('source_stage35_report_sha256') != s35.sha(stage35_single)
                or row.get('source_stage35_summary_sha256') != s35.sha(upstream_path)
                or upstream['models'][model].get('report_sha256') != s35.sha(stage35_single)
                or row.get('source_matching_artifact_sha256') != s35.sha(source_match)):
            fail('Stage36 result/source provenance invalid: ' + model)
        nu = row.get('neutrino_observables', {})
        if (not isinstance(nu.get('masses_ev'), list) or len(nu['masses_ev']) != 3
                or not all(math.isfinite(float(v)) and float(v) >= 0 for v in nu['masses_ev'])):
            fail('Invalid Stage36 neutrino eigenvalues: ' + model)
        rows[model] = {
            'report_sha256': s35.sha(p),
            'source_kind': row['existing_matching_artifact_kind'],
            'masses_ev': nu['masses_ev'],
            'delta_m21_sq_ev2': nu['delta_m21_sq_ev2'],
            'delta_m31_sq_ev2': nu['delta_m31_sq_ev2'],
            'lambdaH_final': row['final_sm_couplings']['lambdaH'],
            'solver_nfev': row['solver_nfev'],
        }
    count_h = sum(v['source_kind'].startswith('hierarchical') for v in rows.values())
    count_l = len(rows) - count_h
    out = {
        'format': SUMMARY_FORMAT,
        'status': SUMMARY_STATUS if len(rows) == 16 else 'PartialExistingC5WeinbergRuns',
        'models_checked': len(rows), 'hierarchical_sources': count_h,
        'legacy_sources': count_l, 'mu_initial_gev': MU_S, 'mu_low_gev': MU_LOW,
        'new_general_quartics_matched': False,
        'physical_coupling_mapping_certified': False,
        'source_stage35_summary_sha256': s35.sha(upstream_path),
        'models': rows,
    }
    p = folder / 'stage36_16_model_low_energy_neutrino_summary.json'
    p.write_text(json.dumps(out, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    options = ap.add_mutually_exclusive_group(required=True)
    options.add_argument('--model', choices=MODELS)
    options.add_argument('--all', action='store_true')
    options.add_argument('--summary', action='store_true')
    ap.add_argument('--allow-legacy', action='store_true',
                    help='Explicitly allow four old one-loop C5-only sources; not identical matching schemes')
    ap.add_argument('--input-dir', type=Path, default=DATA)
    ap.add_argument('--root', type=Path, default=ROOT)
    args = ap.parse_args()
    try:
        if args.summary:
            summary = write_summary(folder=args.input_dir)
            print('Stage36 summary: ' + str(args.input_dir / 'stage36_16_model_low_energy_neutrino_summary.json'))
            return
        choices = MODELS if args.all else (args.model,)
        for m in choices:
            res = run_model(m, folder=args.input_dir, root=args.root,
                            allow_legacy=args.allow_legacy)
            print(f"Stage36 {m}: {res['existing_matching_artifact_kind']}; "
                  f"m_nu[eV]={res['neutrino_observables']['masses_ev']}; "
                  f"nfev={res['solver_nfev']}", flush=True)
        if args.all:
            summary = write_summary(folder=args.input_dir)
            print(f"Stage36 {summary['models_checked']}/16: "
                  + str(args.input_dir / 'stage36_16_model_low_energy_neutrino_summary.json'))
    except (ValueError, KeyError, TypeError, RuntimeError, OSError, OverflowError,
            ArithmeticError, json.JSONDecodeError) as exc:
        ap.exit(2, str(exc) + '\n')


if __name__ == '__main__':
    main()
