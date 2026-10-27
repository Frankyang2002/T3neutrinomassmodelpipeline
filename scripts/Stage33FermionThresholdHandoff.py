"""Stage 33: reconstruct the genuine Stage 30/31 fermion-threshold UV endpoint.

The Stage 30/31 compact reports contain no complete evolved Yukawa/mass state.
Replay their *same* coupled RGBeta ODE from the hashed sources, verify EVERY
saved compact point, and project the reconstructed canonical state to the
existing scalar-only renormalisable EFT boundary at mu_F = 1e5 GeV.

This performs NO Wilson-coefficient matching or intermediate-EFT running.
New general quartics remain an explicit unprojected sidecar, not silently
identified with couplings in the existing intermediate state.

  python -m scripts.Stage33FermionThresholdHandoff --model A_+0
  python -m scripts.Stage33FermionThresholdHandoff --all
"""
from __future__ import annotations

import argparse
from dataclasses import fields, is_dataclass
from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from scripts import Stage30CoupledA0UVRunning as s30
from scripts import Stage31CoupledUVRunning as s31

ROOT = s30.ROOT
DATA = s30.DATA
MODELS = (s31.BASELINE,) + s31.MODELS
FORMAT = 'stage33-t3-fermion-threshold-handoff-v1'


def fail(message: str) -> None:
    raise ValueError('Stage33 fail-closed: ' + message)


def _must_close(actual: Any, expected: Any, label: str, *, rtol: float = 5e-8,
                atol: float = 2e-10) -> None:
    try:
        a, b = float(actual), float(expected)
    except (TypeError, ValueError, OverflowError) as exc:
        fail(f'{label}: cannot compare numbers: {exc}')
    if not (math.isfinite(a) and math.isfinite(b) and
            math.isclose(a, b, rel_tol=rtol, abs_tol=atol)):
        fail(f'{label}: replay={a!r} recorded={b!r}')


def json_state(value: Any) -> Any:
    """Lossless-for-double-precision JSON schema for canonical dataclass states.

    Matrix and scalar complex values use explicit real/imag parts; None stays
    null for absent representation-dependent couplings.
    """
    if is_dataclass(value) and not isinstance(value, type):
        return {f.name: json_state(getattr(value, f.name)) for f in fields(value)}
    if isinstance(value, np.ndarray):
        if value.ndim != 2:
            fail('Unexpected numerical-state array rank')
        arr = np.asarray(value, dtype=np.complex128)
        if not np.all(np.isfinite(arr.real)) or not np.all(np.isfinite(arr.imag)):
            fail('Nonfinite numerical-state matrix')
        return {'real': arr.real.tolist(), 'imag': arr.imag.tolist()}
    if isinstance(value, (complex, np.complexfloating)):
        z = complex(value)
        if not (math.isfinite(z.real) and math.isfinite(z.imag)):
            fail('Nonfinite complex numerical-state scalar')
        return {'real': z.real, 'imag': z.imag}
    if isinstance(value, (float, np.floating)):
        v = float(value)
        if not math.isfinite(v):
            fail('Nonfinite numerical-state real scalar')
        return v
    if isinstance(value, (int, np.integer)) and not isinstance(value, bool):
        return int(value)
    if value is None or isinstance(value, (bool, str)):
        return value
    fail('Unexpected canonical-state field type: ' + type(value).__name__)


def _recorded_path(raw: str) -> Path:
    if not isinstance(raw, str) or not raw.strip():
        fail('Malformed recorded source path')
    path = Path(raw.replace('\\', '/'))
    if not path.is_absolute():
        if len(raw) >= 2 and raw[1] == ':':
            fail('Windows absolute path on non-Windows host: ' + raw)
        path = ROOT / path
    # Provenance never needs traversal outside repository except absolute
    # source keys already explicitly recorded by Stage 30/31.
    return path


def verify_sources(record: dict) -> None:
    mapping = record.get('input_sha256')
    if not isinstance(mapping, dict) or len(mapping) < 5:
        fail('Missing original Stage30/31 source hash inventory')
    for raw_path, expected in mapping.items():
        if not isinstance(expected, str) or len(expected) != 64:
            fail('Malformed original SHA-256 record')
        actual = s30.file_hash(_recorded_path(raw_path))
        if actual != expected:
            fail('Original running input changed since integration: ' + raw_path)


def verify_saved_trajectory(record: dict, solved: Any, layout: Any,
                            names: tuple[str, ...], reused: tuple[str, ...]) -> None:
    """Verify all stored UV sample points, including optional/reused quartics."""
    from Numerical.core.StateVector import unpack_uv_state

    rows = record.get('trajectory')
    if not isinstance(rows, list) or len(rows) < 3 or len(rows) != solved.y.shape[1]:
        fail('Saved trajectory sampling differs from replay')
    if solved.y.shape[0] != layout.size + 2 * len(names):
        fail('Replayed ODE state has incorrect dimension')
    if record.get('n_saved_points') != len(rows):
        fail('Saved trajectory point count inconsistent')
    for i, row in enumerate(rows):
        if not isinstance(row, dict):
            fail('Invalid Stage30/31 compact trajectory point')
        _must_close(math.exp(float(solved.t[i])), row.get('mu_gev'), f'scale[{i}]',
                    rtol=2e-11, atol=1e-6)
        st = unpack_uv_state(solved.y[:layout.size, i], layout,
                             mu_gev=math.exp(float(solved.t[i])))
        for key in ('gY', 'g2', 'g3', 'lambdaH'):
            _must_close(getattr(st.sm, key), row.get(key), f'{key}[{i}]')
        saved_new = row.get('extra_quartics')
        if not isinstance(saved_new, dict) or set(saved_new) != set(names):
            fail('Saved general quartic set differs from replay')
        for j, name in enumerate(names):
            obj = saved_new[name]
            if not isinstance(obj, dict):
                fail('Malformed saved new quartic: ' + name)
            _must_close(solved.y[layout.size+2*j, i], obj.get('real'), name + '.real')
            _must_close(solved.y[layout.size+2*j+1, i], obj.get('imag'), name + '.imag')
        if reused:
            saved_reused = row.get('reused_quartics')
            if not isinstance(saved_reused, dict) or set(saved_reused) != set(reused):
                fail('Saved reused quartics differ from replay')
            for name in reused:
                z = complex(getattr(st, name))
                _must_close(z.real, saved_reused[name].get('real'), name + '.real')
                _must_close(z.imag, saved_reused[name].get('imag'), name + '.imag')
    if record.get('nfev') != solved.nfev:
        fail('Solver evaluation count differs from original integration')


def _verify_report(record: dict, model: str, names: tuple[str, ...],
                   reused: tuple[str, ...]) -> None:
    if record.get('model') != model or \
       record.get('status') != 'SuccessOneModelNumericalUVRunning' or \
       record.get('one_loop_beta_evolution_performed') is not True or \
       record.get('full_coupled_uv_beta_system_evolved') is not True or \
       record.get('physical_coupling_mapping_certified') is not False or \
       record.get('mu_initial_gev') != s30.START or \
       record.get('mu_final_gev') != s30.STOP or \
       record.get('source_one_loop_beta_convention') != s30.BETA_CONVENTION or \
       record.get('new_quartic_couplings_evolved') != len(names) or \
       record.get('sidecar_extra_real_dimensions') != 2 * len(names):
        fail('Original UV report header changed: ' + model)
    if model != s31.BASELINE:
        if record.get('reused_quartic_couplings_evolved') != len(reused) or \
           set(record.get('new_quartic_symbols', [])) != set(names) or \
           set(record.get('reused_quartic_symbols', [])) != set(reused):
            fail('Stage31 reused/new quartic inventory changed: ' + model)


def _verify_summary(summary: dict, model: str, original: Path) -> None:
    if summary.get('format') != 'stage31-16-model-uv-summary-v1' or \
       summary.get('status') != 'Success16ModelConditionalUVRunning' or \
       summary.get('total_models') != 16 or \
       (summary.get('new_quartic_slots_total'),
        summary.get('reused_quartic_slots_total')) != (38, 12) or \
       summary.get('physical_coupling_mapping_certified') is not False:
        fail('Stage31 combined summary is missing or inconsistent')
    key = 'baseline_stage30_sha256' if model == s31.BASELINE else None
    digest = (summary.get(key) if key else
              summary.get('stage31_sha256', {}).get(model))
    if digest != s30.file_hash(original):
        fail('Original UV trajectory hash differs from Stage31 summary: ' + model)


def _replay(model: str, folder: Path, record: dict):
    """Replay original 16-model ODE, no production UV machinery modified."""
    registry = s30.read_json(s30.REGISTRY)
    stage29 = s30.read_json(folder / 'stage29_convention_uv_handoff.json')
    s30.verify_stage29_registry_provenance(stage29)
    if model == s31.BASELINE:
        names = s30.quartic_names(registry)
        reused = ()
        boundary = s30.uv_quartics(stage29, registry)
        config = s30.read_json(s30.CONFIG)
        _st, base_vec, layout = s30.benchmark_state(config)
        export = s30.read_json(folder / 'stage30_A_p0_full_uv_betas.json')
        betas = s30.validate_export(export, names, layout.names)
    else:
        entry, names, reused = s31.inventories(registry, model)
        values = s31.stage29_values(stage29, entry, names, reused)
        conf = s31.read_json(s31.config_path(entry))
        base_vec, layout = s31.prepare_state(conf, entry, reused, values)
        boundary = {name: values[name] for name in names}
        safe = model.replace('+', 'p').replace('-', 'm')
        export = s30.read_json(folder / f'stage31_{safe}_full_uv_betas.json')
        betas = s31.validate_export(export, entry, names, reused, layout.names)

    if record.get('baseline_real_dimensions') != layout.size:
        fail('Original UV vector size differs from reconstructed layout')
    y0 = np.concatenate((base_vec, np.asarray(
        [x for name in names for x in (boundary[name].real, boundary[name].imag)],
        dtype=float)))
    rhs = lambda t, y: s30.combined_rhs(t, y, layout=layout,
                                       names=names, betas=betas)
    _ = rhs(math.log(s30.START), y0)
    solved = s30.integrate_vector(
        y0, rhs, n_points=record['n_saved_points'],
        rtol=float(record['rtol']), atol=float(record['atol']))
    verify_saved_trajectory(record, solved, layout, names, reused)
    from Numerical.core.StateVector import unpack_uv_state
    uv_final = unpack_uv_state(solved.y[:layout.size, -1], layout,
                               mu_gev=s30.STOP).validated()
    extra_final = {
        name: complex(float(solved.y[layout.size+2*j, -1]),
                      float(solved.y[layout.size+2*j+1, -1]))
        for j, name in enumerate(names)}
    return uv_final, extra_final, reused, solved


def handoff(model: str, *, folder: Path = DATA, output: Path | None = None) -> dict:
    if model not in MODELS:
        fail('Unknown ordinary T3 model: ' + model)
    safe = model.replace('+', 'p').replace('-', 'm')
    report_path = folder / (f'stage30_A_p0_coupled_uv_running.json' if model ==
                            s31.BASELINE else f'stage31_{safe}_coupled_uv_running.json')
    report = s30.read_json(report_path)
    summary_path = folder / 'stage31_16_model_uv_summary.json'
    summary = s30.read_json(summary_path)
    _verify_summary(summary, model, report_path)
    registry = s30.read_json(s30.REGISTRY)
    if model == s31.BASELINE:
        new, reused = s30.quartic_names(registry), ()
    else:
        _, new, reused = s31.inventories(registry, model)
    _verify_report(report, model, new, reused)
    verify_sources(report)
    uv_state, extra_final, reused, solved = _replay(model, folder, report)

    from Numerical.running.FermionThresholdBoundary import (
        fermion_threshold_masses, build_intermediate_scalar_boundary)
    fermion = fermion_threshold_masses(uv_state)
    intermediate = build_intermediate_scalar_boundary(
        uv_state, matching_scale_gev=s30.STOP)
    for name in reused:
        _must_close(complex(getattr(uv_state, name)).real,
                    complex(getattr(intermediate, name)).real,
                    name + ' reused threshold real')
        _must_close(complex(getattr(uv_state, name)).imag,
                    complex(getattr(intermediate, name)).imag,
                    name + ' reused threshold imag')
    payload = {
        'format': FORMAT,
        'status': 'SuccessRenormalisableFermionThresholdProjection',
        'model': model,
        'mu_fermion_threshold_gev': s30.STOP,
        'uv_running_replayed_and_trajectory_verified': True,
        'nfev_replay': int(solved.nfev),
        'fermion_threshold_method': 'no-finite-renormalisable-correction; remove F/y1/y2/MF',
        'heavy_fermion_singular_masses_gev': fermion.masses_gev.tolist(),
        'uv_threshold_state': json_state(uv_state),
        'intermediate_scalar_boundary': json_state(intermediate),
        'new_quartic_sidecar_at_threshold': {
            name: json_state(value) for name, value in extra_final.items()},
        'reused_quartics_in_intermediate_state': list(reused),
        'new_quartics_integrated_into_intermediate_eft': False,
        'wilson_coefficient_matching_performed': False,
        'intermediate_eft_running_performed': False,
        'physical_coupling_mapping_certified': False,
        'source_report': report_path.relative_to(folder).as_posix(),
        'source_report_sha256': s30.file_hash(report_path),
        'source_summary_sha256': s30.file_hash(summary_path),
        'limitations': [
            'Renormalisable threshold projection only; no C5/C12 Wilson matching',
            'New quartic sidecar has no representation in the existing scalar-only EFT state',
            'No intermediate-EFT running or neutrino-observable calculation',
            'Matchete/RGBeta quartic component mapping UNCERTIFIED',
            'Stage29 illustrative quartic boundary values, not refitted data',
        ],
    }
    if output is None:
        output = folder / f'stage33_{safe}_fermion_threshold_handoff.json'
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, allow_nan=False) + '\n',
                      encoding='utf-8')
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--model', choices=MODELS)
    group.add_argument('--all', action='store_true')
    parser.add_argument('--input-dir', type=Path, default=DATA)
    ns = parser.parse_args()
    try:
        targets = MODELS if ns.all else (ns.model,)
        result = {}
        for model in targets:
            data = handoff(model, folder=ns.input_dir)
            masses = data['heavy_fermion_singular_masses_gev']
            result[model] = {'status': data['status'],
                             'report_sha256': data['source_report_sha256'],
                             'heavy_fermion_singular_masses_gev': masses,
                             'new_sidecar_count': len(data['new_quartic_sidecar_at_threshold']),
                             'reused_count': len(data['reused_quartics_in_intermediate_state'])}
            print(f'Stage33 {model}: threshold projection PASSED; '
                  f'MF singular values = {masses}; Wilson matching = NOT PERFORMED',
                  flush=True)
        if ns.all:
            combined = {'format': 'stage33-16-model-fermion-threshold-handoff-v1',
                        'status': 'Success16RenormalisableThresholdProjections',
                        'models': result, 'models_checked': len(result),
                        'physical_coupling_mapping_certified': False,
                        'wilson_coefficient_matching_performed': False}
            out = ns.input_dir / 'stage33_16_model_fermion_threshold_summary.json'
            out.write_text(json.dumps(combined, indent=2) + '\n', encoding='utf-8')
            print('Stage33 16/16 threshold projections:', out)
    except (ValueError, KeyError, TypeError, OSError, RuntimeError,
            OverflowError, json.JSONDecodeError) as exc:
        parser.exit(2, str(exc) + '\n')


if __name__ == '__main__':
    main()
