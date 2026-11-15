"""Synthetic Git-fixture tests for Windows EOL-safe thesis output archiver."""
from __future__ import annotations

import csv
import importlib.util
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / 'scripts/maintenance/ArchiveThesis45UnusedOutputs.py'
spec = importlib.util.spec_from_file_location('thesis45_blob_archive', SCRIPT)
assert spec and spec.loader
arch = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = arch
spec.loader.exec_module(arch)


def git(repo: Path, *args: str) -> str:
    result = subprocess.run(['git', '-C', str(repo), *args], check=True,
                            capture_output=True, text=True)
    return result.stdout.strip()


@pytest.fixture()
def repo(tmp_path: Path, monkeypatch):
    root = tmp_path / 'checkout'
    root.mkdir()
    git(root, 'init', '-q')
    git(root, 'config', 'user.email', 'fixture@example.invalid')
    git(root, 'config', 'user.name', 'Fixture')
    git(root, 'config', 'core.autocrlf', 'true')
    (root / 'output' / 'oldrun').mkdir(parents=True)
    (root / 'audits').mkdir()
    (root / 'output' / 'oldrun' / 'data.csv').write_bytes(b'x,y\n1,2\n')
    (root / 'output' / 'oldrun' / 'run.json').write_bytes(b'{"x": 1}\n')
    git(root, 'add', '--', 'output/oldrun')
    git(root, 'commit', '-qm', 'frozen source')
    expected_count = 2
    # Verify the fixture canonical bytes are the ones stated in the manifest.
    expected_bytes = sum(archived_size(root, rel) for rel in (
        'output/oldrun/data.csv', 'output/oldrun/run.json'
    ))
    # Simulate Windows CRLF checkout. Linux Git may flag the manually
    # converted file as dirty despite equivalent canonical blob content.
    # We isolate EOL/shape verification from the independent git-status gate.
    for name in ('data.csv', 'run.json'):
        path = root / 'output' / 'oldrun' / name
        path.write_bytes(path.read_bytes().replace(b'\n', b'\r\n'))
    manifest = root / 'audits' / 'thesis45_archive_candidates.csv'
    with manifest.open('w', encoding='utf-8', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=[
            'path', 'expected_files', 'expected_bytes', 'classification', 'reason', 'audit_commit'
        ])
        writer.writeheader()
        writer.writerow({
            'path': 'output/oldrun', 'expected_files': expected_count,
            'expected_bytes': expected_bytes, 'classification': 'archive_if_unreferenced',
            'reason': 'Synthetic test', 'audit_commit': git(root, 'rev-parse', 'HEAD'),
        })
    monkeypatch.setattr(arch, 'ROOT', root)
    monkeypatch.setattr(arch, 'PLAN', manifest)
    monkeypatch.setattr(arch, 'REPORT', root / 'audits' / 'preflight.csv')
    monkeypatch.setattr(arch, 'OUTSIDE_ARCHIVES', tmp_path / 'external_archive')
    monkeypatch.setattr(arch, 'git_status', lambda path: '')
    return root


def archived_size(root: Path, path: str) -> int:
    proc = subprocess.run(['git', '-C', str(root), 'cat-file', '-s', 'HEAD:' + path],
                          check=True, capture_output=True, text=True)
    return int(proc.stdout.strip())


def inspect():
    return arch.evaluate(arch.load_candidates())[0]


def test_crlf_worktree_is_ready_when_canonical_blob_matches(repo: Path):
    assert b'\r\n' in (repo / 'output/oldrun/data.csv').read_bytes()
    row = inspect()
    assert row['status'] == 'READY'
    assert row['found_bytes'] > row['expected_bytes']
    assert row['committed_bytes'] == row['expected_bytes']


def test_content_change_is_held_even_if_git_blob_inventory_matches(repo: Path):
    path = repo / 'output/oldrun/data.csv'
    path.write_bytes(path.read_bytes().replace(b'1,2', b'9,2'))
    row = inspect()
    assert row['status'] == 'HOLD'
    assert 'local-content-differs-from-HEAD' in row['reason_codes']


def test_extra_local_file_is_held(repo: Path):
    (repo / 'output/oldrun/extra.json').write_text('{}')
    row = inspect()
    assert row['status'] == 'HOLD'
    assert 'local-file-set-differs-from-HEAD' in row['reason_codes']


def test_stale_manifest_is_held(repo: Path):
    path = repo / 'audits/thesis45_archive_candidates.csv'
    text = path.read_text()
    lines = text.splitlines()
    cells = lines[1].split(',')
    cells[2] = str(int(cells[2]) + 5)
    lines[1] = ','.join(cells)
    path.write_text('\n'.join(lines) + '\n')
    row = inspect()
    assert row['status'] == 'HOLD'
    assert 'committed-inventory-differs-from-audited-main' in row['reason_codes']


def test_references_protect_archive(repo: Path):
    (repo / 'scripts').mkdir()
    (repo / 'scripts' / 'uses_source.py').write_text('x = "output/oldrun"\n')
    row = inspect()
    assert row['status'] == 'HOLD'
    assert 'referenced-outside-own-folder' in row['reason_codes']


def test_zip_roundtrip_keeps_exact_local_bytes(repo: Path):
    row = inspect()
    assert row['status'] == 'READY'
    original = {x.name: x.read_bytes() for x in (repo / 'output/oldrun').iterdir()}
    target = arch.OUTSIDE_ARCHIVES / 'test_archive.zip'
    record, paths = arch.make_verified_zip([row], target)
    assert len(paths) == 2
    assert record['physics_source_moved'] is False
    shutil.rmtree(repo / 'output/oldrun')
    arch.restore_zip(target)
    assert {x.name: x.read_bytes() for x in (repo / 'output/oldrun').iterdir()} == original


def test_protected_thesis_path_not_manifest_eligible():
    with pytest.raises(ValueError, match='Protected'):
        arch.clean_rel('output/parallel_t3/eft1_llss_thesis_consistent_v1')
