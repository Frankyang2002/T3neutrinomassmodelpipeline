"""Synthetic-repository safety checks for the local-only T3 bulk cleanup utility."""
from __future__ import annotations

import json
from pathlib import Path
import subprocess

import pytest

from tools import ThesisBulkCleanup as bulk


def git(root: Path, *args: str) -> None:
    subprocess.run(["git", "-C", str(root), *args], check=True,
                   capture_output=True, text=True)


@pytest.fixture
def repo(tmp_path: Path) -> Path:
    root = tmp_path / "study"
    root.mkdir()
    git(root, "init", "-q")
    git(root, "config", "user.email", "cleanup-test@example.invalid")
    git(root, "config", "user.name", "T3 Cleanup Test")
    target = root / "Reports/output/case"
    target.mkdir(parents=True)
    for suffix in (".tex", ".pdf", ".aux", ".out", ".toc", ".log", ".compile.log"):
        (target / ("report" + suffix)).write_text("test " + suffix, encoding="utf-8")
    (target / "orphan.aux").write_text("keep me", encoding="utf-8")
    (root / "output").mkdir()
    (root / "output/physics.log").write_text("precious log", encoding="utf-8")
    (root / "output/scientific.json").write_text("{}", encoding="utf-8")
    git(root, "add", ".")
    git(root, "commit", "-qm", "test baseline")
    return root


def test_selects_only_rebuildable_latex_byproducts(repo: Path) -> None:
    paths = {row["path"] for row in bulk.find_candidates(repo)}
    assert paths == {f"Reports/output/case/report{suffix}" for suffix in
                     (".aux", ".out", ".toc", ".log", ".compile.log")}


def test_preview_does_not_change_worktree(repo: Path) -> None:
    plan = bulk.preview(repo, repo.parent / "external")
    data = json.loads(plan.read_text(encoding="utf-8"))
    assert data["count"] == 5
    assert not bulk.git(repo, "status", "--porcelain", "--untracked-files=no").strip()
    assert (repo / "Reports/output/case/report.aux").exists()


def test_apply_requires_confirmation_and_fresh_head(repo: Path) -> None:
    external = repo.parent / "external"
    plan = bulk.preview(repo, external)
    with pytest.raises(RuntimeError, match="requires --confirm"):
        bulk.apply(repo, external, plan, "WRONG")
    (repo / "Reports/output/case/report.aux").write_text("changed", encoding="utf-8")
    with pytest.raises(RuntimeError, match="not clean"):
        bulk.apply(repo, external, plan, bulk.CONFIRM)


def test_apply_and_restore_verified(repo: Path) -> None:
    outside = repo.parent / "external"
    plan = bulk.preview(repo, outside)
    archived = bulk.apply(repo, outside, plan, bulk.CONFIRM)
    assert not (repo / "Reports/output/case/report.aux").exists()
    assert (repo / "Reports/output/case/report.pdf").exists()
    assert (repo / "Reports/output/case/report.tex").exists()
    assert (repo / "Reports/output/case/orphan.aux").exists()
    assert (repo / "output/scientific.json").exists()
    assert (repo / "output/physics.log").exists()
    assert (archived / "manifest.json").exists()
    bulk.restore(repo, archived)
    assert (repo / "Reports/output/case/report.aux").read_text() == "test .aux"
    assert not bulk.git(repo, "status", "--porcelain", "--untracked-files=no").strip()


def test_restore_preflights_all_conflicts(repo: Path) -> None:
    external = repo.parent / "external"
    archived = bulk.apply(repo, external, bulk.preview(repo, external), bulk.CONFIRM)
    (repo / "Reports/output/case/report.aux").write_text("conflict", encoding="utf-8")
    with pytest.raises(RuntimeError, match="destination exists"):
        bulk.restore(repo, archived)
    assert not (repo / "Reports/output/case/report.out").exists()


def test_duplicates_report_does_not_remove_files(repo: Path) -> None:
    a = repo / "output/copy_a.txt"
    b = repo / "output/copy_b.txt"
    a.write_text("same content", encoding="utf-8")
    b.write_text("same content", encoding="utf-8")
    git(repo, "add", ".")
    git(repo, "commit", "-qm", "add duplicate")
    out = bulk.duplicates(repo, repo.parent / "external")
    assert out.exists()
    assert a.exists() and b.exists()
    summary = json.loads((out.parent / "summary.json").read_text())
    assert summary["duplicate_paths_beyond_one"] >= 1
