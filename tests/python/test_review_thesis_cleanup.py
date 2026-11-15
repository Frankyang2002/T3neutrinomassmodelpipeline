"""Safety tests for local thesis cleanup review and reversible archival."""
from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / 'scripts/maintenance/ReviewThesisCleanup.py'
spec = importlib.util.spec_from_file_location('thesis_review_test_module', SCRIPT)
assert spec and spec.loader
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)


def g(root: Path, *args: str) -> bytes:
    return subprocess.run(['git', '-C', str(root), *args], capture_output=True, check=True).stdout


def setup_repo(tmp_path: Path) -> tuple[Path, Path, dict]:
    root = tmp_path / 'repo'
    root.mkdir()
    g(root, 'init', '-q')
    g(root, 'config', 'user.name', 'Fixture')
    g(root, 'config', 'user.email', 'fixture@example.invalid')
    (root/'README.md').write_text('Current readme', encoding='utf-8')
    (root/'README_BATCH1.md').write_text('Historical note\n', encoding='utf-8')
    (root/'WORKER_B_README.md').write_text('Referenced by code\n', encoding='utf-8')
    (root/'output').mkdir()
    (root/'output'/'thesis_result.json').write_text('{"neutrino_mass": 0.1}', encoding='utf8')
    g(root, 'add', '.')
    g(root, 'commit', '-qm', 'initial')
    archive = tmp_path/'T3_thesis_external_archive'
    head = g(root, 'rev-parse', 'HEAD').decode().strip()
    entries = [
        {'path': 'README.md', 'category': 'KEEP', 'decision': 'KEEP', 'bytes': 14, 'references': []},
        {'path': 'README_BATCH1.md', 'category': 'CANDIDATE', 'decision': 'BLOCKED_REFERENCE',
         'bytes': 16, 'references': ['output/cleanup_audit/plan.json']},
        {'path': 'WORKER_B_README.md', 'category': 'CANDIDATE', 'decision': 'BLOCKED_REFERENCE',
         'bytes': 19, 'references': ['scripts/main.py']},
        {'path': 'output/thesis_result.json', 'category': 'REVIEW', 'decision': 'REVIEW',
         'bytes': 22, 'references': []},
    ]
    # Keep sizes exact for copy safety.
    for x in entries:
        x['bytes'] = (root/x['path']).stat().st_size
    plan = {'schema': module.SCHEMA, 'repository': str(root), 'source_head': head,
            'archive_root': str(archive), 'entries': entries}
    plan_path = tmp_path/'plan.json'
    plan_path.write_text(json.dumps(plan), encoding='utf8')
    return root, plan_path, plan


def prepare_approved(tmp_path: Path, plan: dict, path: str, explanation: str = 'Inspected: historical note only; outside archive kept') -> Path:
    report_dir = tmp_path/'review'
    csvfile, summary = module.analyze(plan, report_dir)
    with csvfile.open(newline='', encoding='utf-8-sig') as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        if row['path'] == path:
            row['approval'] = 'ARCHIVE'
            row['justification'] = explanation
    with csvfile.open('w', newline='', encoding='utf-8-sig') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return csvfile


def test_read_only_review_has_no_approvals(tmp_path: Path):
    root, planpath, plan = setup_repo(tmp_path)
    csvfile, summary = module.analyze(plan, tmp_path/'review')
    assert summary['source_rows'] == 4
    assert summary['archive_default_approvals'] == 0
    with csvfile.open(newline='', encoding='utf-8-sig') as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 2
    assert all(not r['approval'] for r in rows)
    assert (root/'README_BATCH1.md').exists()


def test_synthetic_fixture_references_are_not_live_runtime_consumers():
    row = {'path': 'README_BATCH1.md', 'category': 'CANDIDATE',
           'references': ['tests/python/test_public_repository_cleanup.py']}
    assert module.historical_review(row)['strong_reference_count'] == 0


def test_index_and_test_refs_block_archival():
    for ref in ('docs/research/INDEX.md', 'tests/python/test_research_organization.py', 'scripts/run.py'):
        row = {'path': 'docs/research/model-and-symmetry/WORKER_OLD.md',
               'category': 'CANDIDATE', 'references': [ref]}
        assert module.historical_review(row)['suggestion'] == 'BLOCKED_LIVE_REFERENCE'


def test_thesis_data_and_core_code_never_eligible():
    for path, cat in [('README.md', 'KEEP'), ('output/studies/thesis_qualitative/a.json', 'KEEP'),
                      ('scripts/RunT3Uniform001StageRunning.py', 'REVIEW'),
                      ('Numerical/core/BetaVector.py', 'KEEP'),
                      ('docs/THESIS_SCOPE.md', 'CANDIDATE')]:
        assert module.historical_review({'path': path, 'category': cat}) is None


def test_refuse_unsupported_approvals(tmp_path: Path):
    root, plan_path, plan = setup_repo(tmp_path)
    path = prepare_approved(tmp_path, plan, 'WORKER_B_README.md')
    with pytest.raises(RuntimeError, match='Not eligible'):
        module.apply(root, plan, path, module.CONFIRM)
    assert (root/'WORKER_B_README.md').exists()


def test_refuse_wrong_git_head_and_dirty_worktree(tmp_path: Path):
    root, plan_path, plan = setup_repo(tmp_path)
    approved = prepare_approved(tmp_path, plan, 'README_BATCH1.md')
    wrong = {**plan, 'source_head': '0'*40}
    with pytest.raises(RuntimeError, match='Git HEAD differs'):
        module.apply(root, wrong, approved, module.CONFIRM)
    (root/'README.md').write_text('Edited', encoding='utf8')
    with pytest.raises(RuntimeError, match='dirty'):
        module.apply(root, plan, approved, module.CONFIRM)


def test_refuse_inadequate_justification(tmp_path: Path):
    root, plan_path, plan = setup_repo(tmp_path)
    approved = prepare_approved(tmp_path, plan, 'README_BATCH1.md', 'delete')
    with pytest.raises(RuntimeError, match='Justification'):
        module.apply(root, plan, approved, module.CONFIRM)


def test_archive_roundtrip_and_preserve_results(tmp_path: Path):
    root, plan_path, plan = setup_repo(tmp_path)
    approved = prepare_approved(tmp_path, plan, 'README_BATCH1.md')
    original = (root/'README_BATCH1.md').read_bytes()
    manifest = module.apply(root, plan, approved, module.CONFIRM)
    assert manifest is not None and manifest.exists()
    assert not (root/'README_BATCH1.md').exists()
    assert (root/'output'/'thesis_result.json').is_file()
    assert (root/'WORKER_B_README.md').is_file()
    archivefile = manifest.parent/'README_BATCH1.md'
    assert archivefile.read_bytes() == original
    module.restore(root, manifest)
    assert (root/'README_BATCH1.md').read_bytes() == original


def test_archive_refuses_changed_source_and_unknown_reference(tmp_path: Path):
    root, plan_path, plan = setup_repo(tmp_path)
    approved = prepare_approved(tmp_path, plan, 'README_BATCH1.md')
    (root/'README_BATCH1.md').write_bytes(b'altered\n')
    with pytest.raises(RuntimeError, match='dirty'):
        module.apply(root, plan, approved, module.CONFIRM)
    assert (root/'README_BATCH1.md').exists()
