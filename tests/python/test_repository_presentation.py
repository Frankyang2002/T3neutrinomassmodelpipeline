"""Self-contained safety tests for the reversible presentation organiser."""
from __future__ import annotations

import json
from pathlib import Path
import subprocess

import pytest

from scripts.maintenance import OrganizeRepositoryPresentation as p


def git(root: Path, *args: str) -> str:
    result = subprocess.run(["git", *args], cwd=root, check=True, capture_output=True, text=True)
    return result.stdout.strip()


def fixture_repo(tmp_path: Path) -> Path:
    root = tmp_path / "T3"
    root.mkdir()
    git(root, "init", "-q")
    git(root, "config", "user.name", "T3 test")
    git(root, "config", "user.email", "t3@example.invalid")
    git(root, "config", "core.autocrlf", "false")
    (root / "scripts").mkdir()
    (root / "scripts" / "active.py").write_text("print('active calculation')\n", encoding="utf-8")
    (root / "README.md").write_text("# Scientific overview\n", encoding="utf-8")
    (root / "CLEANUP_README.md").write_bytes(b"Old packaging notes\n")
    (root / "README_REPLACEMENT.md").write_bytes(b"Alternative front page\n")
    (root / "backup_probe_20261031_215644" / "Lagrangian").mkdir(parents=True)
    (root / "backup_probe_20261031_215644" / "Lagrangian" / "old.wl").write_text("old source\n")
    git(root, "add", ".")
    git(root, "commit", "-qm", "baseline")
    return root


def paths(plan: dict, key: str) -> set[str]:
    return {x["path"] for x in plan[key]}


def test_preview_only(tmp_path):
    root = fixture_repo(tmp_path)
    plan = p.write_plan(root)
    assert "CLEANUP_README.md" in paths(plan, "approved")
    assert "backup_probe_20261031_215644" in paths(plan, "approved")
    assert (root / "CLEANUP_README.md").is_file()
    assert (root / p.PLAN_REL).is_file()
    assert (root / "scripts/active.py").is_file()


def test_referenced_file_is_held(tmp_path):
    root = fixture_repo(tmp_path)
    (root / "scripts/active.py").write_text("OLD_DOC = 'README_REPLACEMENT.md'\n")
    git(root, "add", ".")
    git(root, "commit", "-qm", "add active reference")
    plan = p.build_plan(root)
    assert "README_REPLACEMENT.md" in paths(plan, "held")
    assert "CLEANUP_README.md" in paths(plan, "approved")


def test_modified_tracked_file_is_held(tmp_path):
    root = fixture_repo(tmp_path)
    (root / "CLEANUP_README.md").write_text("Modified local notes\n")
    assert "CLEANUP_README.md" in paths(p.build_plan(root), "held")


def test_untracked_file_in_directory_is_held(tmp_path):
    root = fixture_repo(tmp_path)
    (root / "backup_probe_20261031_215644" / "untracked.tmp").write_text("keep me")
    assert "backup_probe_20261031_215644" in paths(p.build_plan(root), "held")


def test_archive_and_restore_with_sha256(tmp_path):
    root = fixture_repo(tmp_path)
    p.write_plan(root)
    archive = p.apply_plan(root, confirm=p.CONFIRM)
    assert not (root / "CLEANUP_README.md").exists()
    assert not (root / "backup_probe_20261031_215644").exists()
    assert (root / "scripts/active.py").is_file()
    manifest = json.loads((archive / "manifest.json").read_text())
    assert manifest["format"] == p.FORMAT
    count = p.restore(root, archive)
    assert count == 3
    assert (root / "CLEANUP_README.md").read_bytes() == b"Old packaging notes\n"
    assert (root / "backup_probe_20261031_215644/Lagrangian/old.wl").is_file()


def test_restore_refuses_overwrite(tmp_path):
    root = fixture_repo(tmp_path)
    p.write_plan(root)
    archive = p.apply_plan(root, confirm=p.CONFIRM)
    (root / "README_REPLACEMENT.md").write_text("New data")
    with pytest.raises(p.PresentationError, match="Refusing to overwrite"):
        p.restore(root, archive)


def test_apply_requires_reviewed_plan(tmp_path):
    root = fixture_repo(tmp_path)
    with pytest.raises(p.PresentationError, match="First run preview"):
        p.apply_plan(root, confirm=p.CONFIRM)
    p.write_plan(root)
    with pytest.raises(p.PresentationError, match="Explicit confirmation"):
        p.apply_plan(root, confirm="not confirmed")


def test_apply_rejects_stale_plan(tmp_path):
    root = fixture_repo(tmp_path)
    p.write_plan(root)
    (root / "CLEANUP_README.md").write_text("Different content")
    with pytest.raises(p.PresentationError, match="stale"):
        p.apply_plan(root, confirm=p.CONFIRM)
    assert (root / "CLEANUP_README.md").exists()


def test_archive_must_be_in_external_sibling(tmp_path):
    root = fixture_repo(tmp_path)
    with pytest.raises(p.PresentationError, match="external archive"):
        p.restore(root, tmp_path / "elsewhere")


def test_windows_crlf_clean_filter_is_allowed(tmp_path):
    root = fixture_repo(tmp_path)
    git(root, "config", "core.autocrlf", "true")
    candidate = root / "CLEANUP_README.md"
    candidate.write_bytes(candidate.read_bytes().replace(b"\n", b"\r\n"))
    assert "CLEANUP_README.md" in paths(p.build_plan(root), "approved")


def test_git_head_change_requires_new_preview(tmp_path):
    root = fixture_repo(tmp_path)
    p.write_plan(root)
    (root / "scripts/active.py").write_text("print('new function')\n")
    git(root, "add", ".")
    git(root, "commit", "-qm", "new script")
    with pytest.raises(p.PresentationError, match="stale"):
        p.apply_plan(root, confirm=p.CONFIRM)
