"""Synthetic regression tests for the conservative one-time-file archiver."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess

import pytest

from tools.AuditOneTimeArtifacts import archive, make_plan, write_report, text_references


def git(root, *args):
    return subprocess.run(['git', '-C', str(root), *args], check=True,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)


def sandbox(tmp_path):
    root = tmp_path / 'project'
    root.mkdir()
    git(root, 'init', '-q')
    git(root, 'config', 'user.email', 'test@example.invalid')
    git(root, 'config', 'user.name', 'Test')
    (root / 'overnight_refinement.log').write_text('historical run: Accepted\n', encoding='utf-8')
    (root / 't -q --tb=line').write_text('Commands marked with * may be preceded by a number.\n', encoding='utf-8')
    (root / 'README_OLD.txt').write_text('old instructions', encoding='utf-8')
    (root / 'output').mkdir()
    (root / 'output' / 'saved_result.json').write_text('{"accepted": 32}', encoding='utf-8')
    (root / 'scripts').mkdir()
    (root / 'scripts' / 'important.py').write_text('def run(): return 42\n', encoding='utf-8')
    git(root, 'add', '-A')
    git(root, 'commit', '-qm', 'initial')
    return root


def by_name(plan, name):
    return next(row for row in plan['targets'] if row['path'] == name)


def test_unreferenced_exact_logs_are_candidates_but_scientific_outputs_not(tmp_path):
    root = sandbox(tmp_path)
    plan = make_plan(root)
    assert {r['path'] for r in plan['targets']} == {'overnight_refinement.log', 't -q --tb=line'}
    assert all(r['classification'] == 'ARCHIVE_CANDIDATE' for r in plan['targets'])
    assert any(r['path'] == 'README_OLD.txt' for r in plan['review_only'])
    assert (root / 'output/saved_result.json').exists()


def test_named_reference_blocks_archival(tmp_path):
    root = sandbox(tmp_path)
    (root / 'scripts' / 'important.py').write_text('print("overnight_refinement.log")\n')
    plan = make_plan(root)
    assert by_name(plan, 'overnight_refinement.log')['classification'] == 'REVIEW_ONLY'
    assert by_name(plan, 'overnight_refinement.log')['references'] == ['scripts/important.py']


def test_modified_tracked_file_is_blocked(tmp_path):
    root = sandbox(tmp_path)
    (root / 'overnight_refinement.log').write_text('new results not yet committed')
    plan = make_plan(root)
    assert by_name(plan, 'overnight_refinement.log')['git_state'] == 'MODIFIED_IN_GIT'
    assert by_name(plan, 'overnight_refinement.log')['classification'] == 'REVIEW_ONLY'


def test_accidental_filename_with_changed_content_is_blocked(tmp_path):
    root = sandbox(tmp_path)
    (root / 't -q --tb=line').write_text('some scientific results')
    git(root, 'add', '-A')
    git(root, 'commit', '-qm', 'update filename content')
    plan = make_plan(root)
    assert by_name(plan, 't -q --tb=line')['reason'] == 'UNEXPECTED_CONTENT_OR_SIZE'


def test_external_archive_copies_verifies_and_unlinks_only_candidates(tmp_path):
    root = sandbox(tmp_path)
    before = hashlib.sha256((root / 'overnight_refinement.log').read_bytes()).hexdigest()
    archive_dir = tmp_path / 'external_archive'
    moved = archive(root, make_plan(root), archive_dir)
    assert moved == ['overnight_refinement.log', 't -q --tb=line']
    assert not (root / 'overnight_refinement.log').exists()
    assert not (root / 't -q --tb=line').exists()
    assert hashlib.sha256((archive_dir / 'overnight_refinement.log').read_bytes()).hexdigest() == before
    manifest = json.loads((archive_dir / 'manifest.json').read_text())
    assert len(manifest['archived_files']) == 2
    assert (root / 'README_OLD.txt').exists()
    assert (root / 'scripts/important.py').exists()
    assert (root / 'output/saved_result.json').exists()


def test_archive_refuses_destination_inside_worktree_or_existing(tmp_path):
    root = sandbox(tmp_path)
    plan = make_plan(root)
    with pytest.raises(RuntimeError, match='outside'):
        archive(root, plan, root / 'archive')
    dest = tmp_path / 'existing'
    dest.mkdir()
    with pytest.raises(RuntimeError, match='already exists'):
        archive(root, plan, dest)
    assert (root / 'overnight_refinement.log').exists()


def test_cross_buffer_reference_finds_name(tmp_path):
    root = sandbox(tmp_path)
    (root / 'docs').mkdir()
    needle = 'overnight_refinement.log'
    (root / 'docs' / 'ref.md').write_bytes(b'x' * (128 * 1024 - 11) + needle.encode() + b'\n')
    refs, _, errors = text_references(root, [needle])
    assert errors == []
    assert refs[needle] == ['docs/ref.md']


def test_report_does_not_modify_source(tmp_path):
    root = sandbox(tmp_path)
    write_report(root, make_plan(root), Path('output/cleanup_audit/one_time_files_v1'))
    report = root / 'output/cleanup_audit/one_time_files_v1/summary.md'
    assert 'README_OLD.txt' in report.read_text()
    assert (root / 'overnight_refinement.log').exists()
    assert make_plan(root)['targets'][0]['classification'] == 'ARCHIVE_CANDIDATE'


def test_scan_errors_fail_closed(tmp_path, monkeypatch):
    root = sandbox(tmp_path)
    import tools.AuditOneTimeArtifacts as audit
    def fake_scan(_root, names):
        return {name: [] for name in names}, 1, ['missing: unavailable']
    monkeypatch.setattr(audit, 'text_references', fake_scan)
    plan = audit.make_plan(root)
    assert all(row['classification'] == 'REVIEW_ONLY' for row in plan['targets'])
