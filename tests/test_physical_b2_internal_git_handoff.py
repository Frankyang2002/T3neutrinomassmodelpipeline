"""Regression tests for Git provenance handoff to native Wolfram without RunProcess.

Synthetic local fixtures; successful tests do NOT certify the physical action.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'Numerical/lo_certification/physical_normalisation'))
from b2_gauge_coupling import AuditBlocked, SOURCE_BLOBS

SPEC = importlib.util.spec_from_file_location(
    'b2_native_launcher', ROOT / 'scripts/RunPhysicalB2RGBetaInternalActionNative.py'
)
assert SPEC is not None and SPEC.loader is not None
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)


def provenance() -> dict:
    return {
        'repository_head_measured': 'c' * 40,
        'all_six_source_blobs_unchanged': True,
        'validated_git_blobs': [
            {'path': rel, 'git_blob_sha': digest,
             'validated_with_git_clean_filter': True}
            for rel, digest in SOURCE_BLOBS.items()
        ],
    }


def test_full_expected_sources_checked_and_env_measured():
    hostile = {k: 'f' * 40 for k in M.ENV_KEYS}
    original = dict(hostile)
    result = M.provenance_environment(provenance(), hostile)
    assert hostile == original
    assert result['T3_B2_PROVENANCE_HEAD'] == 'c' * 40
    assert all(result[key] == SOURCE_BLOBS[rel]
               for rel, key in M.PROVENANCE_ENV.items())
    assert set(M.PROVENANCE_ENV) == set(M.ALLOWED_B2_PATHS)
    assert len(M.ENV_KEYS) == 4


def test_dirty_source_rejected():
    for change in ('head', 'hash', 'flag', 'count', 'duplicate'):
        p = provenance()
        if change == 'head':
            p['repository_head_measured'] = 'no-head'
        elif change == 'hash':
            p['validated_git_blobs'][0]['git_blob_sha'] = '0' * 40
        elif change == 'flag':
            p['validated_git_blobs'][0]['validated_with_git_clean_filter'] = False
        elif change == 'count':
            p['validated_git_blobs'].pop()
        else:
            p['validated_git_blobs'][1]['path'] = p['validated_git_blobs'][0]['path']
        with pytest.raises(AuditBlocked):
            M.provenance_environment(p, {})


def test_not_verified_source_flag_rejected():
    p = provenance()
    p['all_six_source_blobs_unchanged'] = False
    with pytest.raises(AuditBlocked):
        M.provenance_environment(p, {})


def fake_root(tmp_path: Path) -> Path:
    root = tmp_path / 'repo'
    (root / 'scripts').mkdir(parents=True)
    (root / M.WL_SCRIPT).write_text('(* test fixture only *)', encoding='utf-8')
    return root


def test_runner_uses_no_wolfram_git_and_propagates_env(monkeypatch, tmp_path):
    root = fake_root(tmp_path)
    out = root / 'native.json'
    seen = {}
    monkeypatch.setattr(M, 'source_provenance', lambda actual_root: provenance())

    def fake_run(cmd, **kwargs):
        seen.update({'cmd': cmd, **kwargs})
        out.write_text('{"format":"test"}', encoding='utf-8')
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(M.subprocess, 'run', fake_run)
    monkeypatch.setattr(M, 'read_report', lambda path: {'format': 'test'})
    monkeypatch.setattr(M, 'native_snapshot_guard',
                        lambda data, actual_root, snapshot: None)
    monkeypatch.setattr(M, 'validate_native_data',
                        lambda data: {'families_checked': 5,
                                      'native_registry_ratios_measured': 30})
    result = M.run_native(root, out, wolfram='fake-wolfram',
                          env={'T3_B2_PROVENANCE_HEAD': '0'*40})
    assert seen['cmd'] == ['fake-wolfram', '-file', str(root / M.WL_SCRIPT), str(out)]
    assert seen['cwd'] == root
    assert seen['env']['T3_B2_PROVENANCE_HEAD'] == 'c' * 40
    assert result['physical_source_certified'] is False


def test_blocked_native_or_missing_executable_cannot_be_certified(monkeypatch, tmp_path):
    root = fake_root(tmp_path)
    out = root / 'native.json'
    monkeypatch.setattr(M, 'source_provenance', lambda actual_root: provenance())
    monkeypatch.setattr(M.subprocess, 'run',
                        lambda *args, **kwargs: SimpleNamespace(returncode=7))
    with pytest.raises(AuditBlocked, match='exited 7'):
        M.run_native(root, out, env={})

    def fail(*args, **kwargs):
        raise FileNotFoundError('not found')
    monkeypatch.setattr(M.subprocess, 'run', fail)
    with pytest.raises(AuditBlocked, match='cannot launch wolframscript'):
        M.run_native(root, out, env={})


def test_existing_output_refused_before_wolfram(monkeypatch, tmp_path):
    root = fake_root(tmp_path)
    out = root / 'native.json'
    out.touch()
    monkeypatch.setattr(M, 'source_provenance', lambda *_: pytest.fail('should not run'))
    with pytest.raises(AuditBlocked, match='refusing existing'):
        M.run_native(root, out)


def test_restricted_flag_checked_before_native(monkeypatch, tmp_path):
    root = fake_root(tmp_path)
    monkeypatch.setattr(M, 'source_provenance', lambda *_: pytest.fail('should not run'))
    with pytest.raises(AuditBlocked, match='nonrestricted'):
        M.run_native(root, root / 'new.json',
                     env={'T3_EXPERIMENTAL_GENERAL_FERMIONS': '1'})


def test_native_script_has_guard_and_no_external_git_execution():
    text = (ROOT / M.WL_SCRIPT).read_text(encoding='utf-8')
    assert 'RunProcess[' not in text
    assert 'B4Git[' not in text
    for key in M.ENV_KEYS:
        assert f'Environment["{key}"]' in text
    assert 'B4SHA1Q' in text
    assert '"git_head_at_capture"->$B4HeadFromPython' in text
    assert '"git_source_blobs_at_capture"->$B4BlobsFromPython' in text
    assert '"physical_source_certified"->False' in text
    assert 'B4JSONUnexpected[$B4Probe,"$"]' in text
