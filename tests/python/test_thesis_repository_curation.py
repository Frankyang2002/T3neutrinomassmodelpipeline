"""Conservative archival behavior on disposable local Git repositories."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import pytest

SRC = Path(__file__).resolve().parents[2] / 'scripts/maintenance/CurateThesisRepository.py'
spec = importlib.util.spec_from_file_location('thesis_curate', SRC)
assert spec and spec.loader
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)


def git(root: Path, *args: str) -> None:
    subprocess.run(['git', '-C', str(root), *args], check=True, capture_output=True)


def fixture(tmp_path: Path) -> Path:
    root = tmp_path / 'local_checkout'
    root.mkdir()
    git(root, 'init', '-q')
    git(root, 'config', 'user.name', 'Test')
    git(root, 'config', 'user.email', 'test@example.invalid')
    content = {
        'README.md': 'Current thesis project\n',
        'README_BATCH2.md': 'An old research note\n',
        'README_BATCH7.md': 'Another historical note\n',
        'docs/PHYSICS.md': 'Physics: see README_BATCH7.md\n',
        'docs/research/old_experiment.md': 'Old experiment\n',
        'docs/README.md': 'Scientific documentation\n',
        'Reports/output/unused/plot.aux': 'temporary LaTeX aux\n',
        'Reports/output/unused/plot.tex': '\\documentclass{article}\n',
        'Reports/output/thesis_general_5_2/Figure5_2.pdf': 'thesis figure bytes\n',
        'output/studies/thesis_qualitative/reference/result.json': '{}\n',
        'output/batch1/symbolic.json': '{}\n',
        'scripts/RunCoreTests.ps1': 'Write-Host baseline\n',
        'scripts/maintenance/OpaqueWorker.py': 'print("retained")\n',
        'Numerical/parallel_t3/rough_inputs/builder.py': 'print("physics")\n',
    }
    for name, text in content.items():
        path = root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding='utf8')
    git(root, 'add', '.')
    git(root, 'commit', '-qm', 'initial')
    return root


def test_scope_not_based_on_stage_or_worker_names():
    keep, _ = module.file_kind('Numerical/parallel_t3/wilson/strict_one_loop_batch2.py', set())
    assert keep == 'KEEP'
    review, _ = module.file_kind('scripts/maintenance/OpaqueWorker.py', set())
    assert review == 'REVIEW'
    assert module.file_kind('Reports/output/thesis_general_5_2/fig.png', set())[0] == 'KEEP'
    assert module.file_kind('docs/research/old_experiment.md', set())[0] == 'CANDIDATE'


def test_audit_reference_block_and_restore(tmp_path):
    root = fixture(tmp_path)
    archive = tmp_path / 'outside'
    path = module.make_plan(root, archive)
    content = json.loads(path.read_text(encoding='utf8'))
    rows = {r['path']: r for r in content['entries']}
    assert rows['README_BATCH2.md']['decision'] == 'READY_ARCHIVE'
    assert rows['README_BATCH7.md']['decision'] == 'BLOCKED_REFERENCE'
    assert rows['Reports/output/unused/plot.aux']['decision'] == 'READY_ARCHIVE'
    assert rows['output/batch1/symbolic.json']['decision'] == 'REVIEW'
    assert rows['output/studies/thesis_qualitative/reference/result.json']['decision'] == 'KEEP'
    manifest = module.apply(root, path, module.CONFIRM)
    assert manifest.exists()
    assert not (root / 'README_BATCH2.md').exists()
    assert not (root / 'Reports/output/unused/plot.aux').exists()
    assert (root / 'README_BATCH7.md').exists()
    assert (root / 'scripts/maintenance/OpaqueWorker.py').exists()
    assert (root / 'Reports/output/thesis_general_5_2/Figure5_2.pdf').exists()
    module.restore(root, manifest)
    assert (root / 'README_BATCH2.md').read_text(encoding='utf8') == 'An old research note\n'
    assert (root / 'Reports/output/unused/plot.aux').read_text(encoding='utf8') == 'temporary LaTeX aux\n'


def test_requires_clean_checkout_before_archiving(tmp_path):
    root = fixture(tmp_path)
    plan = module.make_plan(root, tmp_path / 'outside')
    (root / 'README.md').write_text('Uncommitted edit\n', encoding='utf8')
    with pytest.raises(RuntimeError, match='Tracked checkout has local changes'):
        module.apply(root, plan, module.CONFIRM)
    assert (root / 'README_BATCH2.md').exists()


def test_incorrect_confirmation_does_nothing(tmp_path):
    root = fixture(tmp_path)
    plan = module.make_plan(root, tmp_path / 'outside')
    with pytest.raises(RuntimeError, match='confirmation'):
        module.apply(root, plan, 'no')
    assert (root / 'README_BATCH2.md').exists()


def test_archive_cannot_be_inside_checkout(tmp_path):
    root = fixture(tmp_path)
    with pytest.raises(RuntimeError, match='outside'):
        module.make_plan(root, root / 'output' / 'archive')


def test_explicit_review_can_archive_unreferenced_non_thesis_data(tmp_path):
    root = fixture(tmp_path)
    outside = tmp_path / 'outside'
    plan = module.make_plan(root, outside)
    reviews = tmp_path / 'reviewed.csv'
    reviews.write_text('path,decision,justification\n'
                       'output/batch1/symbolic.json,ARCHIVE,Checked no thesis figures or executable source consume this saved artifact\n',
                       encoding='utf8')
    manifest = module.apply(root, plan, module.REVIEW_CONFIRM, reviews)
    assert not (root / 'output/batch1/symbolic.json').exists()
    assert (root / 'output/studies/thesis_qualitative/reference/result.json').exists()
    module.restore(root, manifest)
    assert (root / 'output/batch1/symbolic.json').exists()


def test_review_cannot_override_protected_thesis_results(tmp_path):
    root = fixture(tmp_path)
    plan = module.make_plan(root, tmp_path / 'outside')
    reviews = tmp_path / 'reviewed.csv'
    reviews.write_text('path,decision,justification\n'
                       'Reports/output/thesis_general_5_2/Figure5_2.pdf,ARCHIVE,Claim it is obsolete despite exact thesis mapping\n',
                       encoding='utf8')
    with pytest.raises(RuntimeError, match='Only audited REVIEW'):
        module.apply(root, plan, module.REVIEW_CONFIRM, reviews)
