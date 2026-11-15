"""Pure-Python safety tests for the presentation-note archive utility."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import subprocess

MODULE = Path(__file__).resolve().parents[2] / "scripts/maintenance/ConsolidatePublicRepository.py"
_spec = importlib.util.spec_from_file_location("t3_public_archive", MODULE)
assert _spec and _spec.loader
from sys import modules
mod = importlib.util.module_from_spec(_spec)
modules[_spec.name] = mod
_spec.loader.exec_module(mod)


def git(path: Path, *args: str) -> None:
    subprocess.run(["git", "-C", str(path), *args], check=True, capture_output=True)


def fixture(tmp_path: Path) -> Path:
    root = tmp_path / "checkout"
    root.mkdir()
    git(root, "init", "-q")
    git(root, "config", "user.name", "Archive test")
    git(root, "config", "user.email", "archive-test@example.invalid")
    (root / "README.md").write_text("# Current public front page\n", encoding="utf8")
    (root / "README_BATCH1.md").write_text("Superseded note\n", encoding="utf8")
    (root / "WORKER_A_README.md").write_text("Another note\n", encoding="utf8")
    (root / "docs").mkdir()
    (root / "docs" / "REFERENCES.md").write_text("See WORKER_A_README.md\n", encoding="utf8")
    git(root, "add", ".")
    git(root, "commit", "-qm", "test")
    return root


def test_is_candidate_conservative():
    assert mod.is_candidate("README_BATCH1.md")
    assert mod.is_candidate("EFT1_STAGE_A_SOURCE_BRIDGE_STATUS.md")
    assert not mod.is_candidate("README.md")
    assert not mod.is_candidate("scripts/RunBatch.py")
    assert not mod.is_candidate("WORKER_A_MIXING_NUMERIC.json")


def test_plan_apply_restore_and_preserve_referenced_note(tmp_path):
    root = fixture(tmp_path)
    archive_base = tmp_path / "external_archive"
    plan = mod.build_plan(root, archive_base)
    data = json.loads(plan.read_text(encoding="utf8"))
    by_name = {r["path"]: r for r in data["root_note_candidates"]}
    assert by_name["README_BATCH1.md"]["eligible"]
    assert not by_name["WORKER_A_README.md"]["eligible"]
    assert by_name["WORKER_A_README.md"]["incoming_references"] == ["docs/REFERENCES.md"]
    manifest = mod.apply(root, archive_base, plan, mod.CONFIRM)
    assert manifest is not None
    assert not (root / "README_BATCH1.md").exists()
    assert (root / "WORKER_A_README.md").exists()
    assert (root / "README.md").exists()
    mod.restore(root, manifest)
    assert (root / "README_BATCH1.md").read_text(encoding="utf8") == "Superseded note\n"


def test_reject_stale_plan(tmp_path):
    root = fixture(tmp_path)
    external = tmp_path / "outside"
    plan = mod.build_plan(root, external)
    (root / "docs" / "REFERENCES.md").write_text("README_BATCH1.md and WORKER_A_README.md\n", encoding="utf8")
    try:
        mod.apply(root, external, plan, mod.CONFIRM)
    except RuntimeError as ex:
        assert "changed" in str(ex)
    else:
        assert False, "Must reject stale plans"
    assert (root / "README_BATCH1.md").is_file()
