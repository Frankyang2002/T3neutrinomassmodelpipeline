"""Stage 34: run supported post-fermion scalar-only EFT1 RGEs for all 16 T3 models.

Inputs: Stage 33 *individual* complete renormalisable boundaries and
Stage 34 live scalar-only RGBeta exports. Existing numerical EFT1 integrator.

No Wilson matching and no RGE or matching for the Stage 29 *new* quartic
sidecars: those interactions are outside the current EFT1 state/beta inventory.

  python -m scripts.Stage34IntermediateEFTRunning --model A_+0
  python -m scripts.Stage34IntermediateEFTRunning --all
  python -m scripts.Stage34IntermediateEFTRunning --summary
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from scripts import Stage33FermionThresholdHandoff as s33

ROOT = s33.ROOT
DATA = s33.DATA
MODELS = s33.MODELS
FORMAT = 'stage34-ordinary-t3-supported-intermediate-running-v1'
BETA_FORMAT = 'stage34-ordinary-t3-eft1-projected-rgbeta-v1'
START = 1.0e5
STOP = 1.0e3
POINTS = 24
REUSED_NAMES = s33.s31.COMPLEX_REUSED


def fail(message: str) -> None:
    raise ValueError('Stage34 fail-closed: ' + message)


def sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    obj = json.loads(path.read_text(encoding='utf-8'))
    if not isinstance(obj, dict):
        fail(f'Not a JSON object: {path}')
    return obj


def safe_model(model: str) -> str:
    if model not in MODELS:
        fail(f'Unknown ordinary T3 model: {model!r}')
    return model.replace('+', 'p').replace('-', 'm')


def _check_sha(actual: str, expected: Any, label: str) -> None:
    if not isinstance(expected, str) or len(expected) != 64 or actual != expected:
        fail('SHA-256 mismatch: ' + label)


def _complex(obj: Any, name: str) -> complex:
    if not isinstance(obj, dict) or set(obj) != {'real', 'imag'}:
        fail(name + ': complex value requires real/imag')
    try:
        value = complex(float(obj['real']), float(obj['imag']))
    except (ValueError, TypeError) as exc:
        fail(name + ': invalid complex value: ' + str(exc))
    if not (math.isfinite(value.real) and math.isfinite(value.imag)):
        fail(name + ': nonfinite complex value')
    return value


def _matrix(obj: Any, name: str) -> np.ndarray:
    if not isinstance(obj, dict) or set(obj) != {'real', 'imag'}:
        fail(name + ': complex matrix requires real/imag')
    try:
        a = np.asarray(obj['real'], dtype=float)
        b = np.asarray(obj['imag'], dtype=float)
    except (TypeError, ValueError) as exc:
        fail(name + ': invalid matrix: ' + str(exc))
    if a.shape != (3, 3) or b.shape != (3, 3) or not (
            np.isfinite(a).all() and np.isfinite(b).all()):
        fail(name + ': expected finite 3x3 real/imag matrices')
    return a + 1j * b


def _scalar(value: Any, name: str) -> float:
    if type(value) not in (float, int):
        fail(name + ': expected real number')
    v = float(value)
    if not math.isfinite(v):
        fail(name + ': nonfinite real number')
    return v


def decode_intermediate_state(raw: Any):
    """Reconstruct/validate canonical EFT1 state from Stage33's JSON schema."""
    from Numerical.core.State import T3Representation, SMNumericalState
    from Numerical.running.IntermediateScalarState import T3IntermediateScalarState

    if not isinstance(raw, dict) or set(raw) != set(T3IntermediateScalarState.__dataclass_fields__):
        fail('Stage33 intermediate scalar state has unexpected/missing fields')
    rep = raw.get('representation')
    sm = raw.get('sm')
    if not isinstance(rep, dict) or not isinstance(sm, dict):
        fail('Stage33 boundary representation/SM state absent')
    if set(rep) != set(T3Representation.__dataclass_fields__) or (
            rep.get('shared_scalar') is not False):
        fail('Stage33 boundary representation is not ordinary T3')
    # The numeric constructor will further validate representation consistency.
    repr_state = T3Representation(**rep).validated()
    if set(sm) != set(SMNumericalState.__dataclass_fields__):
        fail('Stage33 SM field inventory changed')
    sm_state = SMNumericalState(
        gY=_scalar(sm['gY'], 'sm.gY'),
        g2=_scalar(sm['g2'], 'sm.g2'),
        g3=_scalar(sm['g3'], 'sm.g3'),
        lambdaH=_scalar(sm['lambdaH'], 'sm.lambdaH'),
        yu=_matrix(sm['yu'], 'sm.yu'),
        yd=_matrix(sm['yd'], 'sm.yd'),
        ye=_matrix(sm['ye'], 'sm.ye')).validated()
    scalars = {}
    complex_names = {'lambdaT3'} | set(REUSED_NAMES)
    for name in T3IntermediateScalarState.__dataclass_fields__:
        if name in ('mu_gev', 'representation', 'sm'):
            continue
        val = raw[name]
        scalars[name] = (None if val is None else
                         _complex(val, name) if name in complex_names else
                         _scalar(val, name))
    state = T3IntermediateScalarState(mu_gev=_scalar(raw['mu_gev'], 'mu_gev'),
                                      representation=repr_state, sm=sm_state,
                                      **scalars).validated()
    if not math.isclose(state.mu_gev, START, rel_tol=1e-12):
        fail('Stage33 EFT1 boundary is not at the 1e5 GeV matching scale')
    return state


def verify_boundary(model: str, folder: Path) -> tuple[dict, Any, dict, tuple[str, ...], dict]:
    """Validate Stage33 summary, UV hashes, and scalar-only state for one model."""
    tag = safe_model(model)
    single_path = folder / f'stage33_{tag}_fermion_threshold_handoff.json'
    summary_path = folder / 'stage33_16_model_fermion_threshold_summary.json'
    uv_summary_path = folder / 'stage31_16_model_uv_summary.json'
    single, summary = read_json(single_path), read_json(summary_path)
    if (summary.get('format') != 'stage33-16-model-fermion-threshold-handoff-v1' or
            summary.get('status') != 'Success16RenormalisableThresholdProjections' or
            summary.get('models_checked') != 16 or
            summary.get('physical_coupling_mapping_certified') is not False or
            summary.get('wilson_coefficient_matching_performed') is not False or
            set(summary.get('models', {})) != set(MODELS)):
        fail('Stage33 combined summary is missing/inconsistent')
    if (single.get('format') != s33.FORMAT or
            single.get('status') != 'SuccessRenormalisableFermionThresholdProjection' or
            single.get('model') != model or
            single.get('uv_running_replayed_and_trajectory_verified') is not True or
            single.get('new_quartics_integrated_into_intermediate_eft') is not False or
            single.get('wilson_coefficient_matching_performed') is not False or
            single.get('intermediate_eft_running_performed') is not False or
            single.get('physical_coupling_mapping_certified') is not False or
            single.get('mu_fermion_threshold_gev') != START):
        fail('Stage33 individual handoff is missing/inconsistent')
    row = summary['models'][model]
    if row.get('status') != single['status'] or (
            row.get('report_sha256') != single.get('source_report_sha256')):
        fail('Stage33 individual result and summary disagree')
    mass_a = np.asarray(row.get('heavy_fermion_singular_masses_gev'), dtype=float)
    mass_b = np.asarray(single.get('heavy_fermion_singular_masses_gev'), dtype=float)
    if mass_a.shape != (3,) or mass_b.shape != (3,) or not np.allclose(
            mass_a, mass_b, rtol=1e-12, atol=1e-7):
        fail('Stage33 fermion mass diagnostics disagree')
    source_name = ('stage30_A_p0_coupled_uv_running.json' if model == 'A_+0'
                   else f'stage31_{tag}_coupled_uv_running.json')
    if single.get('source_report') != source_name:
        fail('Stage33 source UV report path is unexpected')
    uv_report_path = folder / source_name
    _check_sha(sha(uv_report_path), single.get('source_report_sha256'),
               'original UV running: ' + model)
    _check_sha(sha(uv_summary_path), single.get('source_summary_sha256'),
               'Stage31 UV summary: ' + model)
    uv_summary = read_json(uv_summary_path)
    s33._verify_summary(uv_summary, model, uv_report_path)

    extra = single.get('new_quartic_sidecar_at_threshold')
    reused = single.get('reused_quartics_in_intermediate_state')
    if not isinstance(extra, dict) or not isinstance(reused, list) or len(set(reused)) != len(reused):
        fail('Stage33 sidecar or reused quartic inventory invalid')
    if (len(extra) != row.get('new_sidecar_count') or
            len(reused) != row.get('reused_count') or
            not set(reused).issubset(REUSED_NAMES) or
            any(not name.startswith('lambdaT3GeneralQ') for name in extra)):
        fail('Stage33 quartic inventory does not match combined summary')
    for name, value in extra.items():
        _complex(value, 'new sidecar '+name)
    state = decode_intermediate_state(single.get('intermediate_scalar_boundary'))
    for key in reused:
        if getattr(state, key) is None:
            fail('Reused quartic missing from scalar-only EFT1 state: ' + key)
    return single, state, extra, tuple(reused), {
        str(single_path): sha(single_path), str(summary_path): sha(summary_path),
        str(uv_report_path): sha(uv_report_path),
        str(uv_summary_path): sha(uv_summary_path),
    }


def validate_beta_payload(payload: dict, state: Any) -> tuple[str, ...]:
    """Precheck exact EFT1 inventory before calling production evaluator."""
    from Numerical.running.IntermediateScalarStateVector import pack_intermediate_scalar_state
    from Numerical.running.IntermediateScalarRGBetaEvaluator import EXPECTED_CONVENTION
    vector, layout = pack_intermediate_scalar_state(state)
    rep = state.representation
    expected = {'SharedScalar': False, 'dS1': rep.d_s1,
                'dS2': rep.d_s2, 'dF': rep.d_f, 'alpha': rep.alpha,
                'IntegratedField': 'F', 'ActiveBSMFields': ['S1', 'S2'],
                'ReportBetaConvention': EXPECTED_CONVENTION}
    meta = payload.get('metadata')
    if (payload.get('format') != BETA_FORMAT or payload.get('status') != 'Success' or
            payload.get('physical_coupling_mapping_certified') is not False or
            payload.get('new_stage29_quartics_integrated') is not False or
            not isinstance(meta, dict) or
            any(meta.get(k) != v for k, v in expected.items())):
        fail('Stage34 scalar-only RGBeta beta export/metadata mismatch')
    betas = payload.get('report_betas')
    if (not isinstance(betas, dict) or set(betas) != set(layout.names) or
            set(payload.get('beta_symbols', ())) != set(layout.names) or
            any(not isinstance(expr, str) or not expr.strip() or
                any(x in expr for x in ('BetaTerm[', 'Missing[', '$Aborted', '$Failed'))
                for expr in betas.values())):
        fail('Stage34 scalar-only RGBeta beta symbol inventory mismatch')
    if not np.all(np.isfinite(vector)):
        fail('Nonfinite EFT1 initial vector')
    return tuple(layout.names)


def _snapshot(state: Any, reused: tuple[str, ...]) -> dict:
    result = {'mu_gev': state.mu_gev, 'gY': state.sm.gY, 'g2': state.sm.g2,
              'g3': state.sm.g3, 'lambdaH': state.sm.lambdaH,
              'mS1Sq': state.mS1Sq, 'mS2Sq': state.mS2Sq,
              'lambdaS1': state.lambdaS1, 'lambdaS2': state.lambdaS2,
              'lambdaH1': state.lambdaH1, 'lambdaH2': state.lambdaH2,
              'lambda12': state.lambda12,
              'lambdaT3': s33.json_state(state.lambdaT3),
              'reused_quartics': {x: s33.json_state(getattr(state, x)) for x in reused}}
    return result


def run_model(model: str, *, folder: Path = DATA, output: Path | None = None,
              points: int = POINTS) -> dict:
    from Numerical.running.IntermediateScalarRunner import run_intermediate_scalar_segment
    from Numerical.running.IntermediateScalarStateVector import pack_intermediate_scalar_state
    single, state, extra, reused, hashes = verify_boundary(model, folder)
    tag = safe_model(model)
    export_path = folder / f'stage34_{tag}_eft1_full_betas.json'
    payload = read_json(export_path)
    names = validate_beta_payload(payload, state)
    if type(points) is not int or points < 3 or points > 1000:
        fail('Requested trajectory point count must be integer in [3,1000]')
    grid = np.geomspace(START, STOP, points)
    trajectory = run_intermediate_scalar_segment(
        state, payload, STOP, save_scales_gev=grid,
        rtol=1e-8, atol=1e-11, method='DOP853')
    if not trajectory.solver_success or trajectory.mu_gev.size != points:
        fail('Scalar-only intermediate solver did not complete expected trajectory')
    if not np.all(np.isfinite(trajectory.y)) or not np.all(
            np.diff(trajectory.mu_gev) < 0):
        fail('Intermediate EFT trajectory contains nonfinite/incorrect scales')
    # Equality of input boundary and first solver point is required, not inferred.
    initial_vec, _layout = pack_intermediate_scalar_state(state)
    if not np.allclose(trajectory.initial_vector, initial_vec, rtol=1e-11, atol=1e-12):
        fail('EFT1 solver initial state differs from Stage33 threshold boundary')
    last = trajectory.final_state.validated()
    if not math.isclose(last.mu_gev, STOP, rel_tol=1e-12):
        fail('EFT1 solver endpoint differs from scalar threshold')
    from Numerical.running.ScalarThresholdBoundary import scalar_threshold_masses
    scalar_masses = scalar_threshold_masses(last)
    rows = [_snapshot(trajectory.state_at_index(i), reused) for i in range(points)]
    report = {
        'format': FORMAT,
        'status': 'SuccessSupportedRenormalisableIntermediateRunning',
        'model': model,
        'mu_initial_gev': START, 'mu_final_gev': STOP,
        'one_loop_scalar_only_rgbeta_evolution_performed': True,
        'native_quartic_class_projected': True,
        'eft1_field_content': ['SM', 'S1', 'S2'],
        'active_wilson_operator_running_performed': False,
        'fermion_wilson_matching_performed': False,
        'physical_coupling_mapping_certified': False,
        'new_stage29_quartics_evolved_in_eft1': False,
        'truncated_new_quartic_beta_system': bool(extra),
        'excluded_new_quartic_symbols': list(extra),
        'new_quartic_sidecar_at_mu_f_only': single['new_quartic_sidecar_at_threshold'],
        'reused_quartics_evolved_in_eft1': list(reused),
        'beta_symbol_count': len(names), 'beta_symbols': list(names),
        'solver': 'scipy.solve_ivp/DOP853',
        'solver_nfev': int(trajectory.nfev),
        'n_saved_points': points,
        'rtol': 1e-8, 'atol': 1e-11,
        'intermediate_initial_state': s33.json_state(state),
        'intermediate_final_state': s33.json_state(last),
        'scalar_threshold_singular_masses_gev': [float(m) for m in scalar_masses.masses_gev],
        'trajectory': rows,
        'input_sha256': {**hashes, str(export_path): sha(export_path)},
        'limitations': [
            'Supported renormalisable EFT1 subsystem only; no C5/C12 Wilson matching/running',
            'Stage29 new quartic sidecars not evolved below mu_F: no certified EFT1 map or beta system',
            'Stage33 no-finite-renormalisable-correction boundary convention retained',
            'Physical Matchete/RGBeta quartic map uncertified',
            'No final scalar matching or neutrino observable calculation',
        ],
    }
    if output is None:
        output = folder / f'stage34_{tag}_intermediate_eft_running.json'
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    return report


def write_summary(*, folder: Path = DATA) -> dict:
    """Aggregate only validated 16-model Stage34 artifacts and source hashes."""
    rows = {}
    for model in MODELS:
        tag = safe_model(model)
        path = folder / f'stage34_{tag}_intermediate_eft_running.json'
        report = read_json(path)
        if (report.get('format') != FORMAT or
                report.get('status') != 'SuccessSupportedRenormalisableIntermediateRunning' or
                report.get('model') != model or
                report.get('mu_initial_gev') != START or
                report.get('mu_final_gev') != STOP or
                report.get('physical_coupling_mapping_certified') is not False or
                report.get('new_stage29_quartics_evolved_in_eft1') is not False or
                report.get('fermion_wilson_matching_performed') is not False or
                report.get('active_wilson_operator_running_performed') is not False):
            fail('Stage34 per-model numerical report invalid: ' + model)
        src = report.get('input_sha256', {})
        if not isinstance(src, dict) or not src or any(
                sha(Path(p)) != digest for p, digest in src.items()):
            fail('Stage34 numerical input provenance stale: ' + model)
        initial, final = report['intermediate_initial_state'], report['intermediate_final_state']
        rows[model] = {'report_sha256': sha(path),
                       'new_quartics_unmatched': len(report['excluded_new_quartic_symbols']),
                       'reused_quartics_evolved': len(report['reused_quartics_evolved_in_eft1']),
                       'lambdaH_initial': initial['sm']['lambdaH'],
                       'lambdaH_final': final['sm']['lambdaH'],
                       'solver_nfev': report['solver_nfev']}
    if len(rows) != 16 or sum(v['new_quartics_unmatched'] for v in rows.values()) != 38 or \
            sum(v['reused_quartics_evolved'] for v in rows.values()) != 12:
        fail('Stage34 global registry coverage totals differ from 38 new/12 reused')
    summary = {'format': 'stage34-16-model-supported-intermediate-running-v1',
               'status': 'Success16SupportedIntermediateSegments',
               'models_checked': 16, 'mu_initial_gev': START, 'mu_final_gev': STOP,
               'new_quartics_not_in_intermediate_eft': 38,
               'reused_quartics_evolved': 12,
               'physical_coupling_mapping_certified': False,
               'wilson_coefficient_matching_performed': False,
               'models': rows}
    output = folder / 'stage34_16_model_intermediate_eft_summary.json'
    output.write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    return summary


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    group = ap.add_mutually_exclusive_group(required=True)
    group.add_argument('--model', choices=MODELS)
    group.add_argument('--all', action='store_true')
    group.add_argument('--summary', action='store_true')
    ap.add_argument('--input-dir', type=Path, default=DATA)
    args = ap.parse_args()
    try:
        if args.summary:
            data = write_summary(folder=args.input_dir)
            print(f"Stage34 {data['models_checked']}/16 intermediate EFT segments: "
                  f"{args.input_dir / 'stage34_16_model_intermediate_eft_summary.json'}")
            return
        for model in (MODELS if args.all else (args.model,)):
            data = run_model(model, folder=args.input_dir)
            print(f"Stage34 {model}: supported EFT1 running PASSED; "
                  f"beta couplings={data['beta_symbol_count']}, "
                  f"nfev={data['solver_nfev']}, "
                  f"new quartics unmatched={len(data['excluded_new_quartic_symbols'])}; "
                  'Wilson matching = NOT PERFORMED', flush=True)
        if args.all:
            write_summary(folder=args.input_dir)
            print('Stage34 combined summary: '
                  + str(args.input_dir / 'stage34_16_model_intermediate_eft_summary.json'))
    except (ValueError, KeyError, TypeError, OSError, RuntimeError,
            OverflowError, json.JSONDecodeError) as exc:
        ap.exit(2, str(exc)+'\n')


if __name__ == '__main__':
    main()
