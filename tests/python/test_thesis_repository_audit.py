"""Synthetic, isolated regression tests for thesis repository inventory/retirement."""
from __future__ import annotations

import json
from pathlib import Path
import subprocess

import pytest

from tools import ThesisRepositoryAudit as a


def run(root: Path, *args: str) -> None:
    subprocess.run(["git", *args], cwd=root, check=True, capture_output=True)


@pytest.fixture
def repo(tmp_path: Path) -> Path:
    root = tmp_path / "repo"
    root.mkdir()
    run(root, "init", "-q")
    run(root, "config", "user.name", "Test User")
    run(root, "config", "user.email", "test@example.com")
    run(root, "config", "core.autocrlf", "false")
    files = {
        "README.md": "T3 thesis project.\n",
        "pipeline.py": "from scripts import UsedByPhysics\n",
        "Lagrangian/Model.wl": "(* model *)\n",
        "Numerical/State.py": "# state\n",
        "scripts/UsedByPhysics.py": '"""Real production script."""\nx = 1\n',
        "scripts/UsedByTests.py": '"""Tested helper."""\nx = 1\n',
        "scripts/UnusedStudy.py": '"""Old numerical experiment."""\nx = 1\n',
        "scripts/UsedByResult.py": '"""Fitted source."""\nx = 1\n',
        "scripts/RunCoreTests.ps1": "Write-Output 'core'\n",
        "tools/AuditLegacyNotes.py": '"""Old root note migration."""\n',
        "tools/MoveStudyFigures.py": '"""One-time plot mover."""\n',
        "tools/RebuildReports.py": '"""Current figure reconstruction."""\n',
        "tools/AuditScriptDependencies.py": '"""Current script scanner."""\n',
        "tests/test_used.py": "from scripts import UsedByTests\n",
        "output/fits/ref.json": '{"runner":"scripts/UsedByResult.py"}',
        "docs/old_study.md": "Run tools/MoveStudyFigures.py for old figures.\n",
        "t -q --tb=line": "SUMMARY OF LESS COMMANDS\nCommands marked with * may be preceded by a number\nExit.\n",
    }
    for rel, data in files.items():
        p = root / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(data, encoding="utf-8", newline="\n")
    run(root, "add", "-A")
    run(root, "commit", "-qm", "fixture")
    return root


def row(report: dict, path: str) -> dict:
    return next(r for r in report["rows"] if r["path"] == path)


def test_role_classification_is_physics_aware() -> None:
    assert a.file_role("Lagrangian/model/x.wl") == "SCIENTIFIC_IMPLEMENTATION"
    assert a.file_role("RGE/running/foo.py") == "SCIENTIFIC_IMPLEMENTATION"
    assert a.file_role("Numerical/fitting/x.py") == "SCIENTIFIC_IMPLEMENTATION"
    assert a.file_role("output/foo.json") == "GENERATED_OR_SAVED_RESULT"
    assert a.file_role("tests/test_x.py") == "VERIFICATION"
    assert a.file_role("docs/research/eft-matching/x.md") == "RESEARCH_DOCUMENT"
    assert a.file_role("scripts/maintenance/foo.py") == "REPOSITORY_MAINTENANCE"


def test_path_escape_rejected(repo: Path) -> None:
    for value in ("../outside", "/etc/passwd", "tools/../../bad", "tools\\file.py"):
        with pytest.raises(ValueError):
            a.checked_path(repo, value)


def test_complete_file_inventory_and_references(repo: Path) -> None:
    r = a.evaluate(repo)
    assert r["coverage_complete"]
    assert r["tracked_files"] == 17
    assert len(r["files"]) == 17
    assert row(r, "scripts/UsedByPhysics.py")["status"] == "KEEP_REFERENCED"
    assert row(r, "scripts/UsedByTests.py")["status"] == "KEEP_REFERENCED"
    assert row(r, "scripts/UsedByResult.py")["status"] == "KEEP_PROVENANCE"
    assert row(r, "scripts/UnusedStudy.py")["status"] == "REVIEW_STANDALONE"
    assert row(r, "tools/AuditLegacyNotes.py")["status"] == "SUGGEST_ARCHIVE"
    assert row(r, "tools/RebuildReports.py")["status"] == "KEEP"
    assert row(r, "scripts/RunCoreTests.ps1")["status"] == "KEEP"
    assert row(r, "t -q --tb=line")["status"] == "SUGGEST_ARCHIVE"
    assert row(r, "tools/MoveStudyFigures.py")["doc_refs"] == 1
    assert "tools/MoveStudyFigures.py" not in r["suggested_selection"]
    assert {"tools/AuditLegacyNotes.py", "t -q --tb=line"} <= set(r["suggested_selection"])


def test_preview_preserves_user_edited_selection(repo: Path) -> None:
    r = a.evaluate(repo)
    out = a.emit(repo, r)
    selection = out / "reviewed_selection.txt"
    selection.write_text("# manual\n", encoding="utf-8")
    a.emit(repo, r)
    assert selection.read_text(encoding="utf-8") == "# manual\n"
    assert (out / "all_files.csv").is_file()
    assert (out / "scripts_and_tools.csv").is_file()
    assert (out / "reference_witnesses.csv").is_file()
    assert (out / "plan.json").is_file()
    assert (out / "summary.md").is_file()


def test_archive_verify_and_restore_roundtrip(repo: Path) -> None:
    r = a.evaluate(repo)
    a.emit(repo, r)
    before = {p: (repo / p).read_bytes() for p in ("tools/AuditLegacyNotes.py", "t -q --tb=line")}
    folder = a.archive(repo, r, list(before), allow_research=False, allow_documented=False)
    assert folder is not None
    assert (folder / "manifest.json").is_file()
    assert all(not (repo / p).exists() for p in before)
    manifest = json.loads((folder / "manifest.json").read_text())
    assert len(manifest["entries"]) == 2
    assert manifest["git_head"] == r["head"]
    a.restore(repo, folder)
    for p, data in before.items():
        assert (repo / p).read_bytes() == data


def test_cannot_archive_scientific_dependency(repo: Path) -> None:
    r = a.evaluate(repo)
    with pytest.raises(RuntimeError, match="not archive-eligible"):
        a.archive(repo, r, ["scripts/UsedByPhysics.py"], allow_research=True, allow_documented=True)


def test_no_auto_archiving_unreferenced_research_code(repo: Path) -> None:
    r = a.evaluate(repo)
    with pytest.raises(RuntimeError, match="needs --allow-research-scripts"):
        a.archive(repo, r, ["scripts/UnusedStudy.py"], allow_research=False, allow_documented=True)
    assert (repo / "scripts/UnusedStudy.py").exists()


def test_explicit_research_script_can_be_restored(repo: Path) -> None:
    r = a.evaluate(repo)
    folder = a.archive(repo, r, ["scripts/UnusedStudy.py"], allow_research=True, allow_documented=False)
    assert folder is not None
    a.restore(repo, folder)
    assert (repo / "scripts/UnusedStudy.py").exists()


def test_documented_path_requires_explicit_override(repo: Path) -> None:
    r = a.evaluate(repo)
    with pytest.raises(RuntimeError, match="--allow-documented"):
        a.archive(repo, r, ["tools/MoveStudyFigures.py"], allow_research=False, allow_documented=False)
    assert (repo / "tools/MoveStudyFigures.py").exists()


def test_dirty_file_held(repo: Path) -> None:
    (repo / "tools/AuditLegacyNotes.py").write_text("modified\n", encoding="utf-8")
    fresh = a.evaluate(repo)
    assert row(fresh, "tools/AuditLegacyNotes.py")["status"] == "HOLD_LOCAL_CHANGE"


def test_git_head_change_rejects_stale_plan(repo: Path) -> None:
    r = a.evaluate(repo)
    (repo / "README.md").write_text("new readme\n", encoding="utf-8")
    run(repo, "add", "README.md")
    run(repo, "commit", "-qm", "change")
    with pytest.raises(RuntimeError, match="Git HEAD changed"):
        a.archive(repo, r, ["tools/AuditLegacyNotes.py"], allow_research=False, allow_documented=False)


def test_new_untracked_reference_blocks_a_review_candidate(repo: Path) -> None:
    original = a.evaluate(repo)
    (repo / "extra_runner.py").write_text("from tools import AuditLegacyNotes\n", encoding="utf-8")
    with pytest.raises(RuntimeError, match="not archive-eligible"):
        a.archive(repo, original, ["tools/AuditLegacyNotes.py"], allow_research=False, allow_documented=False)


def test_restore_conflict_and_hash_tamper_rejected(repo: Path) -> None:
    r = a.evaluate(repo)
    folder = a.archive(repo, r, ["tools/AuditLegacyNotes.py"], allow_research=False, allow_documented=False)
    assert folder is not None
    (repo / "tools/AuditLegacyNotes.py").write_text("user data", encoding="utf-8")
    with pytest.raises(RuntimeError, match="overwrite"):
        a.restore(repo, folder)
    (repo / "tools/AuditLegacyNotes.py").unlink()
    (folder / "files/tools/AuditLegacyNotes.py").write_text("tampered", encoding="utf-8")
    with pytest.raises(RuntimeError, match="checksum"):
        a.restore(repo, folder)


def test_unreadable_oversized_text_disables_suggestions(repo: Path) -> None:
    sizes = a.git_names(repo, "ls-files", "-z")
    targets = a.targets_from_inventory(sizes)
    _, _, issues = a.reference_scan(repo, sizes, [], targets, max_text_bytes=6)
    assert any(r["problem"] == "oversized_text" for r in issues)


def test_cli_refuses_apply_without_explicit_confirmation(repo: Path) -> None:
    assert a.main(["--root", str(repo)]) == 0
    with pytest.raises(SystemExit):
        a.main(["--root", str(repo), "--apply"])


def test_sha_matches_git_head(repo: Path) -> None:
    assert a.is_git_blob_unchanged(repo, "tools/AuditLegacyNotes.py")
    assert a.is_git_blob_unchanged(repo, "t -q --tb=line")
    (repo / "tools/AuditLegacyNotes.py").write_text("drift", encoding="utf-8")
    assert not a.is_git_blob_unchanged(repo, "tools/AuditLegacyNotes.py")


def test_tool_index_is_not_counted_as_execution_dependency(repo: Path) -> None:
    index = repo / "tools/README.md"
    index.write_text("All tools: AuditLegacyNotes.py, MoveStudyFigures.py\n", encoding="utf-8")
    r = a.evaluate(repo)
    assert row(r, "tools/AuditLegacyNotes.py")["status"] == "SUGGEST_ARCHIVE"
    assert row(r, "tools/MoveStudyFigures.py")["doc_refs"] == 1


def test_known_retirement_mentions_are_not_runtime_uses_of_pager(repo: Path) -> None:
    p = repo / "tools/AuditOneTimeArtifacts.py"
    p.write_text("TARGETS = {'t -q --tb=line': 'pager capture'}\n", encoding="utf-8")
    r = a.evaluate(repo)
    assert row(r, "t -q --tb=line")["status"] == "SUGGEST_ARCHIVE"
