"""Stage 35: verified scalar-threshold handoff for the 16 ordinary T3 models.

Consume Stage 34 full intermediate-EFT endpoints at mu_S=1e3 GeV and call
Numerical.running.ScalarThresholdBoundary.build_final_sm_boundary.

This is a renormalisable *projection*: it removes scalar BSM fields without
finite scalar-threshold corrections. No C5/C12 Wilson matching, Weinberg
running, or neutrino observable calculation is performed. New Stage29 quartics
remain excluded from the Stage34 EFT and are carried only as a provenance-
checked mu_F sidecar; no extrapolated mu_S value is invented.

  python -m scripts.Stage35ScalarThresholdHandoff --model A_+0
  python -m scripts.Stage35ScalarThresholdHandoff --all
  python -m scripts.Stage35ScalarThresholdHandoff --summary
"""
from __future__ import annotations

import argparse
from dataclasses import replace
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from scripts import Stage34IntermediateEFTRunning as s34

DATA = s34.DATA
ROOT = s34.ROOT
MODELS = s34.MODELS
START = s34.START
STOP = s34.STOP
FORMAT = 'stage35-ordinary-t3-scalar-threshold-handoff-v1'
SUMMARY_FORMAT = 'stage35-16-model-scalar-threshold-handoff-v1'
STATUS = 'SuccessRenormalisableScalarThresholdProjection'
SUMMARY_STATUS = 'Success16RenormalisableScalarThresholdProjections'


def fail(message: str) -> None:
    raise ValueError('Stage35 fail-closed: ' + message)


def read_json(path: Path) -> dict[str, Any]:
    return s34.read_json(path)


def sha(path: Path) -> str:
    return s34.sha(path)


def _number(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        fail(label + ': expected numeric value')
    v = float(value)
    if not math.isfinite(v):
        fail(label + ': nonfinite value')
    return v


def _close(a: Any, b: Any, label: str, *, rtol: float = 3e-10,
           atol: float = 3e-11) -> None:
    lhs = _number(a, label + '.actual')
    rhs = _number(b, label + '.recorded')
    if not math.isclose(lhs, rhs, rel_tol=rtol, abs_tol=atol):
        fail(f'{label}: state={lhs!r} record={rhs!r}')


def _compare_json(actual: Any, recorded: Any, label: str) -> None:
    """Compare all numerical snapshot fields, not just the Higgs quartic."""
    if isinstance(actual, dict):
        if not isinstance(recorded, dict) or set(actual) != set(recorded):
            fail(label + ': snapshot key inventory changed')
        for key in actual:
            _compare_json(actual[key], recorded[key], label + '.' + key)
        return
    if isinstance(actual, (int, float)) and not isinstance(actual, bool):
        _close(actual, recorded, label)
        return
    if actual != recorded:
        fail(label + ': snapshot value changed')


def _source_path(raw: str) -> Path:
    if not isinstance(raw, str) or not raw.strip():
        fail('Malformed recorded Stage34 input path')
    path = Path(raw.replace('\\', '/'))
    if not path.is_absolute():
        if len(raw) >= 2 and raw[1] == ':':
            fail('Cross-platform unresolved Windows source path: ' + raw)
        path = ROOT / path
    return path


def verify_input_hashes(report: dict[str, Any]) -> None:
    mapping = report.get('input_sha256')
    if not isinstance(mapping, dict) or len(mapping) < 5:
        fail('Stage34 input SHA-256 inventory missing')
    for raw, digest in mapping.items():
        if not isinstance(digest, str) or len(digest) != 64 or any(
                c not in '0123456789abcdef' for c in digest):
            fail('Stage34 source hash is malformed: ' + str(raw))
        path = _source_path(raw)
        if not path.is_file() or sha(path) != digest:
            fail('Stage34 source file changed or unavailable: ' + str(raw))


def decode_final_state(raw: Any):
    """Reuse Stage34's strict canonical-state JSON decoder at mu_S.

    That decoder explicitly accepts only mu_F. Check the actual scale first,
    decode an otherwise identical state with that input-scale convention, then
    replace *only* the scale metadata and revalidate the canonical state.
    No physical couplings are changed.
    """
    if not isinstance(raw, dict):
        fail('Stage34 final state must be an object')
    _close(raw.get('mu_gev'), STOP, 'EFT1 final-state scale', rtol=1e-12, atol=0)
    input_scale_copy = dict(raw)
    input_scale_copy['mu_gev'] = START
    initial_scale_state = s34.decode_intermediate_state(input_scale_copy)
    state = replace(initial_scale_state, mu_gev=STOP).validated()
    return state


def _checked_report(model: str, folder: Path) -> tuple[dict, dict, Path, Path]:
    """Validate Stage34 individual/aggregate flags, report hash and all inputs."""
    tag = s34.safe_model(model)
    report_path = folder / f'stage34_{tag}_intermediate_eft_running.json'
    summary_path = folder / 'stage34_16_model_intermediate_eft_summary.json'
    report, summary = read_json(report_path), read_json(summary_path)
    if (summary.get('format') != 'stage34-16-model-supported-intermediate-running-v1'
            or summary.get('status') != 'Success16SupportedIntermediateSegments'
            or summary.get('models_checked') != 16
            or summary.get('mu_initial_gev') != START
            or summary.get('mu_final_gev') != STOP
            or summary.get('new_quartics_not_in_intermediate_eft') != 38
            or summary.get('reused_quartics_evolved') != 12
            or summary.get('physical_coupling_mapping_certified') is not False
            or summary.get('wilson_coefficient_matching_performed') is not False
            or not isinstance(summary.get('models'), dict)
            or set(summary['models']) != set(MODELS)):
        fail('Stage34 global summary invalid')
    row = summary['models'][model]
    if (not isinstance(row, dict)
            or row.get('report_sha256') != sha(report_path)):
        fail('Stage34 numerical report SHA-256 differs from summary: ' + model)
    required = {
        'format': s34.FORMAT,
        'status': 'SuccessSupportedRenormalisableIntermediateRunning',
        'model': model,
        'mu_initial_gev': START,
        'mu_final_gev': STOP,
        'one_loop_scalar_only_rgbeta_evolution_performed': True,
        'native_quartic_class_projected': True,
        'eft1_field_content': ['SM', 'S1', 'S2'],
        'active_wilson_operator_running_performed': False,
        'fermion_wilson_matching_performed': False,
        'physical_coupling_mapping_certified': False,
        'new_stage29_quartics_evolved_in_eft1': False,
    }
    if any(report.get(k) != v for k, v in required.items()):
        fail('Stage34 individual report status/convention invalid: ' + model)
    extras = report.get('excluded_new_quartic_symbols')
    reused = report.get('reused_quartics_evolved_in_eft1')
    sidecar = report.get('new_quartic_sidecar_at_mu_f_only')
    if (not isinstance(extras, list) or not isinstance(reused, list)
            or not isinstance(sidecar, dict)
            or any(not isinstance(s, str) for s in extras + reused)
            or len(extras) != len(set(extras))
            or len(reused) != len(set(reused))
            or set(extras) != set(sidecar)
            or len(extras) != row.get('new_quartics_unmatched')
            or len(reused) != row.get('reused_quartics_evolved')
            or report.get('truncated_new_quartic_beta_system') is not bool(extras)
            or report.get('solver_nfev') != row.get('solver_nfev')
            or report.get('n_saved_points', 0) < 3):
        fail('Stage34 quartic/solver inventory mismatch: ' + model)
    for name, complex_value in sidecar.items():
        if not name.startswith('lambdaT3GeneralQ'):
            fail('Unexpected sidecar coupling name')
        s34._complex(complex_value, 'mu_F sidecar ' + name)
    _close(report.get('intermediate_initial_state', {}).get('sm', {}).get('lambdaH'),
           row.get('lambdaH_initial'), 'Stage34 initial Higgs quartic')
    _close(report.get('intermediate_final_state', {}).get('sm', {}).get('lambdaH'),
           row.get('lambdaH_final'), 'Stage34 final Higgs quartic')
    verify_input_hashes(report)
    return report, row, report_path, summary_path


def run_model(model: str, *, folder: Path = DATA, output: Path | None = None) -> dict:
    if model not in MODELS:
        fail('Unknown ordinary T3 model: ' + str(model))
    folder = Path(folder)
    report, _row, stage34_path, summary_path = _checked_report(model, folder)
    init_state = s34.decode_intermediate_state(report['intermediate_initial_state'])
    final_state = decode_final_state(report['intermediate_final_state'])
    if init_state.representation != final_state.representation:
        fail('Intermediate representation changed during Stage34 evolution')
    reused = tuple(report['reused_quartics_evolved_in_eft1'])
    for name in reused:
        if getattr(init_state, name, None) is None or getattr(final_state, name, None) is None:
            fail('Reused quartic missing from EFT1 state: ' + name)
    trajectory = report.get('trajectory')
    if (not isinstance(trajectory, list)
            or len(trajectory) != report['n_saved_points']
            or not all(isinstance(x, dict) for x in trajectory)):
        fail('Stage34 saved trajectory inventory invalid')
    scales = [_number(x.get('mu_gev'), 'trajectory mu_gev') for x in trajectory]
    if (not all(scales[i] > scales[i+1] for i in range(len(scales)-1))):
        fail('Stage34 saved trajectory scales are not descending')
    _close(scales[0], START, 'Stage34 trajectory start scale', rtol=1e-12, atol=0)
    _close(scales[-1], STOP, 'Stage34 trajectory final scale', rtol=1e-12, atol=0)
    _compare_json(s34._snapshot(init_state, reused), trajectory[0], 'EFT1 start')
    _compare_json(s34._snapshot(final_state, reused), trajectory[-1], 'EFT1 end')

    from Numerical.running.ScalarThresholdBoundary import (
        build_final_sm_boundary, scalar_threshold_masses)
    mass = scalar_threshold_masses(final_state)
    saved_masses = report.get('scalar_threshold_singular_masses_gev')
    if not isinstance(saved_masses, list) or len(saved_masses) != 2:
        fail('Stage34 scalar threshold mass inventory missing')
    for i, (actual, expected) in enumerate(zip(mass.masses_gev, saved_masses)):
        _close(actual, expected, f'Stage34 scalar mass[{i}]', rtol=3e-10, atol=1e-7)
    final_sm = build_final_sm_boundary(final_state, matching_scale_gev=STOP)
    if final_sm.mu_gev != STOP:
        fail('Scalar-only -> SM threshold projection changed matching scale')
    for name in ('gY', 'g2', 'g3', 'lambdaH', 'yu', 'yd', 'ye'):
        source = getattr(final_state.sm, name)
        target = getattr(final_sm.sm, name)
        if not np.array_equal(np.asarray(source), np.asarray(target)):
            fail('Finite/unexpected renormalisable threshold change: ' + name)

    output_report = {
        'format': FORMAT,
        'status': STATUS,
        'model': model,
        'mu_scalar_threshold_gev': STOP,
        'stage34_intermediate_running_verified': True,
        'scalar_threshold_method': 'no-finite-renormalisable-correction; remove S1/S2',
        'scalar_threshold_masses_gev': [float(x) for x in mass.masses_gev],
        'intermediate_final_state_at_mu_s': s34.s33.json_state(final_state),
        'final_sm_boundary_at_mu_s': s34.s33.json_state(final_sm),
        'new_quartics_unmatched_at_mu_f': list(report['excluded_new_quartic_symbols']),
        'new_quartic_sidecar_at_mu_f_only': report['new_quartic_sidecar_at_mu_f_only'],
        'reused_quartics_evolved_to_mu_s': list(reused),
        'new_quartics_matched_or_evolved_below_mu_f': False,
        'scalar_wilson_matching_performed': False,
        'fermion_wilson_matching_performed': False,
        'low_energy_weinberg_running_performed': False,
        'neutrino_observables_calculated': False,
        'physical_coupling_mapping_certified': False,
        'source_stage34_report': stage34_path.name,
        'source_stage34_report_sha256': sha(stage34_path),
        'source_stage34_summary_sha256': sha(summary_path),
        'limitations': [
            'Renormalisable scalar removal only; no finite scalar matching correction',
            'C5/C12 Wilson coefficient matching not performed',
            '38 new Stage29 quartic slots are excluded from intermediate EFT beta system',
            'New quartic sidecar has a value only at mu_F, not at mu_S',
            'Physical Matchete/RGBeta coupling mapping uncertified',
            'No SM+Weinberg evolution or neutrino-mass prediction',
        ],
    }
    target = output or (folder / f'stage35_{s34.safe_model(model)}_scalar_threshold_handoff.json')
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(output_report, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    return output_report


def write_summary(*, folder: Path = DATA) -> dict:
    folder = Path(folder)
    stage34_summary_path = folder / 'stage34_16_model_intermediate_eft_summary.json'
    upstream = read_json(stage34_summary_path)
    if (upstream.get('format') != 'stage34-16-model-supported-intermediate-running-v1'
            or upstream.get('status') != 'Success16SupportedIntermediateSegments'
            or set(upstream.get('models', {})) != set(MODELS)):
        fail('Stage34 global summary invalid for Stage35 aggregation')
    rows = {}
    for model in MODELS:
        tag = s34.safe_model(model)
        path = folder / f'stage35_{tag}_scalar_threshold_handoff.json'
        item = read_json(path)
        if (item.get('format') != FORMAT
                or item.get('status') != STATUS
                or item.get('model') != model
                or item.get('mu_scalar_threshold_gev') != STOP
                or item.get('stage34_intermediate_running_verified') is not True
                or item.get('scalar_wilson_matching_performed') is not False
                or item.get('fermion_wilson_matching_performed') is not False
                or item.get('new_quartics_matched_or_evolved_below_mu_f') is not False
                or item.get('low_energy_weinberg_running_performed') is not False
                or item.get('physical_coupling_mapping_certified') is not False
                or item.get('neutrino_observables_calculated') is not False):
            fail('Stage35 per-model result invalid: ' + model)
        upstream_path = folder / f'stage34_{tag}_intermediate_eft_running.json'
        if (item.get('source_stage34_report') != upstream_path.name
                or item.get('source_stage34_report_sha256') != sha(upstream_path)
                or item.get('source_stage34_summary_sha256') != sha(stage34_summary_path)
                or upstream['models'][model].get('report_sha256') != sha(upstream_path)):
            fail('Stage35 upstream SHA-256 mismatch: ' + model)
        masses = item.get('scalar_threshold_masses_gev')
        if (not isinstance(masses, list) or len(masses) != 2
                or not all(math.isfinite(float(x)) and float(x) > 0 for x in masses)):
            fail('Stage35 scalar mass inventory invalid: ' + model)
        final = item.get('final_sm_boundary_at_mu_s')
        if not isinstance(final, dict) or final.get('mu_gev') != STOP:
            fail('Stage35 final SM boundary missing: ' + model)
        rows[model] = {
            'report_sha256': sha(path),
            'scalar_threshold_masses_gev': masses,
            'lambdaH_final': final['sm']['lambdaH'],
            'new_quartics_unmatched_at_mu_f': len(item['new_quartics_unmatched_at_mu_f']),
            'reused_quartics_evolved_to_mu_s': len(item['reused_quartics_evolved_to_mu_s']),
        }
    if (len(rows) != 16
            or sum(r['new_quartics_unmatched_at_mu_f'] for r in rows.values()) != 38
            or sum(r['reused_quartics_evolved_to_mu_s'] for r in rows.values()) != 12):
        fail('Stage35 full 38/12 scalar inventory does not match Stage34')
    summary = {
        'format': SUMMARY_FORMAT, 'status': SUMMARY_STATUS,
        'models_checked': 16, 'mu_scalar_threshold_gev': STOP,
        'renormalisable_scalar_threshold_projections': 16,
        'new_quartics_unmatched_at_mu_f': 38,
        'reused_quartics_evolved_to_mu_s': 12,
        'scalar_wilson_matching_performed': False,
        'physical_coupling_mapping_certified': False,
        'source_stage34_summary_sha256': sha(stage34_summary_path),
        'models': rows,
    }
    path = folder / 'stage35_16_model_scalar_threshold_summary.json'
    path.write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--model', choices=MODELS)
    group.add_argument('--all', action='store_true')
    group.add_argument('--summary', action='store_true')
    parser.add_argument('--input-dir', type=Path, default=DATA)
    args = parser.parse_args()
    try:
        if args.summary:
            result = write_summary(folder=args.input_dir)
            print(f"Stage35 {result['models_checked']}/16 scalar handoffs: "
                  f"{args.input_dir / 'stage35_16_model_scalar_threshold_summary.json'}")
            return
        targets = MODELS if args.all else (args.model,)
        for model in targets:
            result = run_model(model, folder=args.input_dir)
            print(f"Stage35 {model}: scalar threshold projection PASSED; "
                  f"mS masses = {result['scalar_threshold_masses_gev']}; "
                  'C5/C12 matching = NOT PERFORMED', flush=True)
        if args.all:
            write_summary(folder=args.input_dir)
            print('Stage35 16/16 scalar handoffs: '
                  + str(args.input_dir / 'stage35_16_model_scalar_threshold_summary.json'))
    except (ValueError, KeyError, TypeError, OSError, RuntimeError,
            OverflowError, json.JSONDecodeError) as exc:
        parser.exit(2, str(exc) + '\n')


if __name__ == '__main__':
    main()
