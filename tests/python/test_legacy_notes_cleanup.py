"""Tests for the conservative legacy root-note archival process."""
from __future__ import annotations
import json
from pathlib import Path
import subprocess
import sys

import pytest

PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT))
from tools.AuditLegacyNotes import (archive_notes, discover_notes, plan,
                                    read_sources, sha256, tracked_clean)


def git(root, *args):
    return subprocess.run(['git', *args], cwd=root, check=True,
                          capture_output=True, text=True)


@pytest.fixture
def repo(tmp_path):
    root = tmp_path / 'project'
    root.mkdir()
    git(root, 'init')
    git(root, 'config', 'user.email', 'test@example.invalid')
    git(root, 'config', 'user.name', 'Test')
    (root/'pipeline.py').write_text('print("ok")\n')
    (root/'README_NOTE.txt').write_text('Old instructions\n')
    (root/'README_REFERENCED.txt').write_text('Referenced instructions\n')
    (root/'README_THESIS_COMPARISONS.txt').write_text('Thesis provenance\n')
    (root/'docs').mkdir()
    (root/'docs'/'topic.md').write_text('See README_REFERENCED.txt for history.\n')
    (root/'output').mkdir()
    (root/'output'/'noisy.json').write_text('"README_NOTE.txt"\n')
    git(root, 'add', '.')
    git(root, 'commit', '-m', 'test-fixture')
    return root


def test_note_discovery(repo):
    assert discover_notes(repo) == ['README_NOTE.txt', 'README_REFERENCED.txt',
                                     'README_THESIS_COMPARISONS.txt']


def test_references_exclude_output_artifacts(repo):
    refs, count, errors = read_sources(repo, discover_notes(repo))
    assert not errors
    assert refs['README_NOTE.txt'] == []
    assert len(refs['README_REFERENCED.txt']) == 1


def test_classifications(repo):
    rows = {x['path']:x for x in plan(repo)['notes']}
    assert rows['README_NOTE.txt']['classification'] == 'ARCHIVE_CANDIDATE'
    assert rows['README_REFERENCED.txt']['classification'] == 'REVIEW_REFERENCED'
    assert rows['README_THESIS_COMPARISONS.txt']['classification'] == 'RETAIN_THESIS_PROVENANCE'


def test_uncommitted_notes_never_archive(repo):
    (repo/'README_NOTE.txt').write_text('Edited locally\n')
    rows = {x['path']:x for x in plan(repo)['notes']}
    assert rows['README_NOTE.txt']['classification'] == 'REVIEW_LOCAL_CHANGES'


def test_verified_archive_preserves_referenced_and_thesis(repo, tmp_path):
    archive = tmp_path/'external'
    old_sha = sha256(repo/'README_NOTE.txt')
    assert archive_notes(repo, plan(repo), archive) == 1
    assert not (repo/'README_NOTE.txt').exists()
    assert (repo/'README_REFERENCED.txt').exists()
    assert (repo/'README_THESIS_COMPARISONS.txt').exists()
    assert sha256(archive/'README_NOTE.txt') == old_sha
    manifest = json.loads((archive/'manifest.json').read_text())
    assert len(manifest['archived_files']) == 1
    assert manifest['archived_files'][0]['path'] == 'README_NOTE.txt'


def test_archive_exists_refused(repo, tmp_path):
    archive = tmp_path/'existing'; archive.mkdir()
    with pytest.raises(RuntimeError, match='already exists'):
        archive_notes(repo, plan(repo), archive)


def test_archive_inside_project_refused(repo):
    with pytest.raises(RuntimeError, match='outside project root'):
        archive_notes(repo, plan(repo), repo/'badarchive')


def test_references_block_previous_clean_note(repo):
    (repo/'docs'/'second.md').write_text('README_NOTE.txt\n')
    rows = {x['path']:x for x in plan(repo)['notes']}
    assert rows['README_NOTE.txt']['classification'] == 'REVIEW_REFERENCED'


def test_modified_archive_checksum_rejected(repo, tmp_path):
    audit = plan(repo)
    (repo/'README_NOTE.txt').write_text('Changed since inventory\n')
    with pytest.raises(RuntimeError, match='Source changed'):
        archive_notes(repo, audit, tmp_path/'external')
    assert (repo/'README_NOTE.txt').exists()
