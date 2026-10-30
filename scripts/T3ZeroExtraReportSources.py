"""Verified readers for the saved 16 x 2 zero-additional-scalar control.

New report-only module; never mutates the original physics data or GitHub.
"""
from __future__ import annotations

import json
from pathlib import Path

from scripts import T3Frozen64Sources as policy

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'output/studies/general_scalar_extended/zero_extra_all16'
MODES = ('held_zero', 'boundary_zero')
NARROW_GAMMA = tuple(round(.95 + .01*i, 2) for i in range(11))
FORMAT = 't3-all16-zero-extra-scalar-control-v1'


def run_root(run_id: str) -> Path:
    policy.require(bool(policy.RUN_ID.fullmatch(run_id)), 'invalid run ID')
    return BASE / run_id


def load_summary(run_id: str) -> dict:
    base = run_root(run_id)
    summary = policy.read_json(base / 'summary.json')
    settings = policy.read_json(base / 'control_settings.json')
    policy.require(summary.get('format') == FORMAT and summary.get('settings') == settings,
                   'zero-extra control summary/settings do not match')
    models = tuple(settings.get('models', ()))
    modes = tuple(settings.get('modes', ()))
    policy.require(models == policy.MODELS and modes == MODES,
                   'this report requires all 16 models in both saved modes')
    bench = settings.get('benchmark')
    policy.require(bench in ('optimised_low','optimised_high'),
                   'zero-extra source is not an optimised benchmark')
    rows = summary.get('rows')
    policy.require(isinstance(rows, list) and len(rows) == 32,
                   'expected 32 zero-extra results, not partial source data')
    index = {}
    for row in rows:
        k = (row.get('mode'), row.get('model'))
        policy.require(k[0] in MODES and k[1] in policy.MODELS and k not in index,
                       'duplicate or invalid zero-extra result row')
        policy.require(row.get('benchmark') == bench, 'inconsistent benchmark')
        index[k] = row
    policy.require(len(index) == 32, 'missing saved control result')
    policy.require(not summary.get('failures'), 'control campaign contains failures')
    return {'benchmark': bench, 'settings': settings, 'rows': index,
            'summary_sha256': policy.digest(base / 'summary.json'),
            'settings_sha256': policy.digest(base / 'control_settings.json')}


def sources(run_id: str, meta: dict, mode: str, model: str) -> dict:
    """Check exact named files and their internal numeric/provenance invariants."""
    policy.require(mode in MODES and model in policy.MODELS, 'invalid control cell')
    bench = meta['benchmark']
    base = run_root(run_id)
    folder = base / bench / mode
    tp = folder / policy.tag(model) / 'continuous_trajectory.json'
    np = folder / 'neutrino_v1' / (policy.tag(model) + '_neutrino.json')
    t, n = policy.read_json(tp), policy.read_json(np)
    th, nh = policy.digest(tp), policy.digest(np)
    settings = t.get('settings', {})
    policy.require(t.get('model') == model and n.get('model') == model and
                   t.get('boundary_mode') == 'FROZEN_UV' and
                   t.get('zero_extra_control_mode') == mode and
                   settings.get('zero_extra_control_mode') == mode and
                   settings.get('benchmark') == bench and
                   settings.get('frozen_uv') is True and
                   t.get('uv_ode_evaluations') == 0 and
                   t.get('scalar_endpoint_physically_admissible') is True,
                   f'{mode}/{model}: not the saved frozen zero-extra numerical calculation')
    policy.require(n.get('source_continuous_sha256') == th and
                   n.get('full_general_scalar_physical_c5_certified') is False and
                   n.get('sm_weinberg_rge_performed') is True,
                   f'{mode}/{model}: central neutrino provenance incompatible')
    policy.require(float(meta['rows'][mode, model]['lambdaT3_input_abs']) > 0,
                   'original lambdaT3 was accidentally zeroed')
    for key in ('extra_quartic_initial_magnitudes',):
        values = t.get(key, {})
        policy.require(isinstance(values, dict) and all(float(x) == 0. for x in values.values()),
                       f'{mode}/{model}: additional quartic initial value not zero')
    q = t.get('quadratic_mu2_at_fermion_threshold_gev2')
    policy.require(q is not None and float(q['real']) == 0 and float(q['imag']) == 0,
                   f'{mode}/{model}: quadratic initial boundary not zero')
    if mode == 'held_zero':
        policy.require(float(meta['rows'][mode, model]['max_extra_quartic_endpoint']) == 0. and
                       float(meta['rows'][mode, model]['abs_mu2_endpoint_gev2']) == 0.,
                       f'{mode}/{model}: held-zero restriction violated')
    return {'cell': {'model': model, 'benchmark': bench, 'mode': mode,
                     'case': 'zero_extra_' + mode, 'yukawa_mode': 'optimised',
                     'status': 'NEUTRINO_COMPLETE',
                     'trajectory': policy.relative(ROOT, tp),
                     'trajectory_sha256': th,
                     'neutrino': policy.relative(ROOT, np),
                     'neutrino_sha256': nh},
            'trajectory': t, 'neutrino': n,
            'source_root': folder}
