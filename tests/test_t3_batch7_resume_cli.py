"""Batch 7 CLI resume regressions using synthetic saved output files.

These tests do not require any Windows physics batch output. They exercise the
actual runner entrypoint and its CLI arguments, protecting reuse of B_-1.
"""
from __future__ import annotations

import json
from pathlib import Path

from scripts import RunT3PoleClosureBatch7 as runner


def _completed(root: Path, out: Path, model: str, *, running: bool = True):
    path = out / runner.tag(model) / 'source_pole_and_matching_ledger.json'
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({'model': model, 'synthetic': True}), encoding='utf-8')
    return {
        'model': model,
        'status': runner.STATUS,
        'low_energy_running_computed': running,
        'mu_low_gev': 173.0,
        'vev_gev': 246.22,
        'output_sha256': runner.hash_file(path),
        'output_path': str(path.relative_to(root)),
    }


def _manifest(out: Path, records: dict):
    out.mkdir(parents=True, exist_ok=True)
    (out/'summary.json').write_text(
        json.dumps({'schema': runner.SCHEMA, 'results': records}), encoding='utf-8')


def test_resume_skips_completed_bminus1_and_retains_record(tmp_path, monkeypatch):
    out = tmp_path / 'output' / 'parallel_t3' / 'native_pole_closure_all16_batch7_v1'
    b = _completed(tmp_path, out, 'B_-1')
    _manifest(out, {'B_-1': b})
    monkeypatch.setattr(runner, 'MODELS', ('B_-1',))

    def forbidden(*args, **kwargs):
        raise AssertionError('B_-1 should not be recomputed on --resume')

    monkeypatch.setattr(runner, 'run_one', forbidden)
    rc = runner.main(['--root', str(tmp_path), '--resume'])
    assert rc == 0
    summary = json.loads((out/'summary.json').read_text(encoding='utf-8'))
    assert summary['completed_count'] == 1
    assert summary['results']['B_-1']['output_sha256'] == b['output_sha256']


def test_resume_runs_missing_model_but_not_completed_bminus1(tmp_path, monkeypatch):
    out = tmp_path / 'output' / 'parallel_t3' / 'native_pole_closure_all16_batch7_v1'
    b = _completed(tmp_path, out, 'B_-1')
    _manifest(out, {'B_-1': b})
    monkeypatch.setattr(runner, 'MODELS', ('B_-1', 'A_-4'))
    processed = []

    def fake_run(root, output, model, *, mu_low, vev, skip_running):
        processed.append(model)
        assert (mu_low, vev, skip_running) == (173.0, 246.22, False)
        return _completed(root, output, model)

    monkeypatch.setattr(runner, 'run_one', fake_run)
    assert runner.main(['--root', str(tmp_path), '--resume']) == 0
    assert processed == ['A_-4']
    summary = json.loads((out/'summary.json').read_text(encoding='utf-8'))
    assert summary['completed_count'] == 2
    assert summary['results']['B_-1']['output_sha256'] == b['output_sha256']


def test_resume_requires_running_only_when_requested(tmp_path, monkeypatch):
    out = tmp_path / 'output' / 'parallel_t3' / 'native_pole_closure_all16_batch7_v1'
    existing = _completed(tmp_path, out, 'B_-1', running=False)
    _manifest(out, {'B_-1': existing})
    monkeypatch.setattr(runner, 'MODELS', ('B_-1',))

    def forbidden(*args, **kwargs):
        raise AssertionError('Pole-only --resume should reuse pole-only output')

    monkeypatch.setattr(runner, 'run_one', forbidden)
    assert runner.main(['--root', str(tmp_path), '--resume',
                        '--skip-low-energy-running']) == 0

    recomputed = []

    def fake_run(root, output, model, *, mu_low, vev, skip_running):
        recomputed.append(model)
        return _completed(root, output, model, running=True)

    monkeypatch.setattr(runner, 'run_one', fake_run)
    assert runner.main(['--root', str(tmp_path), '--resume']) == 0
    assert recomputed == ['B_-1']
