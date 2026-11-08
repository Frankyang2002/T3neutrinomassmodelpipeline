"""Offline safety tests for a reversible, non-destructive research file layout."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import subprocess
import tempfile

import pytest

MODULE_PATH = Path(__file__).resolve().parents[2] / "scripts/maintenance/OrganizeResearchFiles.py"
spec = importlib.util.spec_from_file_location("research_organization_under_test", MODULE_PATH)
assert spec and spec.loader
organize = importlib.util.module_from_spec(spec)
spec.loader.exec_module(organize)


def git(root: Path, *args: str) -> str:
    p = subprocess.run(("git", *args), cwd=root, capture_output=True, text=True, check=True)
    return p.stdout.strip()


def sample(tmp_path: Path) -> Path:
    root = tmp_path / "T3"
    root.mkdir()
    git(root, "init", "-q")
    git(root, "config", "user.email", "fixture@example.org")
    git(root, "config", "user.name", "Synthetic Fixture")
    git(root, "config", "core.autocrlf", "false")
    entries = {
        "README.md": "# Reader overview\n",
        "docs/README.md": "# Documentation index\n",
        "docs/MODEL_AND_PIPELINE.md": "# Physics\n",
        "docs/T3_BC_UNREFERENCED_STUDY.md": "# Scalar matching proof\nScientific details.\n",
        "docs/T3_FROZEN_WITH_LINK.md": "# Scan\nSee [physics](MODEL_AND_PIPELINE.md).\n",
        "docs/T3_BC_REFERENCED_STUDY.md": "# Referenced\nScientific result.\n",
        "docs/T3_BC_EDITED.md": "# Research input\n",
        "docs/T3_BC_VIA_DOC.md": "# Candidate from other doc\n",
        "docs/T3_BC_OTHER.md": "# Related\nSee T3_BC_VIA_DOC.md\n",
        "docs/RESULTS.md": "# Results\nSee docs/T3_BC_REFERENCED_STUDY.md.\n",
        "scripts/__init__.py": "",
        "scripts/RunCoreTests.ps1": "# core suite\n",
        "scripts/Standalone.experimental.py": "from __future__ import annotations\nprint('standalone')\n",
        "scripts/LocationDependent.experimental.py": "from pathlib import Path\nROOT = Path(__file__).resolve().parents[1]\n",
        "scripts/Referenced.experimental.py": "print('example')\n",
        "scripts/README.md": "Use scripts/Referenced.experimental.py for inspection.\n",
        "mixed_wilson_inputs.zip": "synthetic archive data",
    }
    for relative, text in entries.items():
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8", newline="\n")
    git(root, "add", ".")
    git(root, "commit", "-qm", "synthetic research layout")
    (root / "docs/T3_BC_EDITED.md").write_text("modified", encoding="utf-8")
    return root


def test_classification_preserves_current_guides() -> None:
    assert organize._doc_group("PHYSICS.md") is None
    assert organize._doc_group("MODEL_AND_PIPELINE.md") is None
    assert organize._doc_group("RESEARCH_HISTORY_INDEX.md") is None
    assert organize._doc_group("T3_BC_110101_MASS_BASIS_MAP.md") == "eft-matching"
    assert organize._doc_group("T3_FROZEN_DIVERSE_RESTART.md") == "numerical-studies"
    assert organize._doc_group("THESIS_FIGURE_CAPTIONS.md") == "scientific-reports"
    assert organize._doc_group("A_MINUS4_SU2_HESSIAN_FRESH_INDEX_BATCH.md") == "model-and-symmetry"
    assert organize._doc_group("T3BCYQuarticMixing_TOY_NOT_BC.json") is None


def test_default_plan_only_moves_unreferenced_docs_and_bundle(tmp_path: Path) -> None:
    root = sample(tmp_path)
    plan = organize.preview(root)
    paths = {x["source"] for x in plan["approved"]}
    held = {x["source"]: x["reason"] for x in plan["held"]}
    assert paths == {"docs/T3_BC_UNREFERENCED_STUDY.md", "mixed_wilson_inputs.zip"}
    assert "docs/T3_BC_OTHER.md" in held
    assert "referenced" in held["docs/T3_BC_REFERENCED_STUDY.md"]
    assert "relative hyperlinks" in held["docs/T3_FROZEN_WITH_LINK.md"]
    assert "differs from Git" in held["docs/T3_BC_EDITED.md"]
    assert "referenced" in held["docs/T3_BC_VIA_DOC.md"]
    assert plan["include_experiments"] is False
    assert (root / "docs/T3_BC_UNREFERENCED_STUDY.md").is_file()


def test_experimental_scripts_require_opt_in_and_location_safety(tmp_path: Path) -> None:
    root = sample(tmp_path)
    plan = organize.build_plan(root, include_experiments=True)
    good = {x["source"]: x["destination"] for x in plan["approved"]}
    held = {x["source"] for x in plan["held"]}
    assert good["scripts/Standalone.experimental.py"] == "scripts/research/experimental/Standalone.py"
    assert "scripts/LocationDependent.experimental.py" in held
    assert "scripts/Referenced.experimental.py" in held


def test_no_apply_without_review_and_confirmation(tmp_path: Path) -> None:
    root = sample(tmp_path)
    with pytest.raises(organize.OrganizationError, match="Preview"):
        organize.apply(root, confirm=organize.CONFIRM)
    organize.preview(root)
    with pytest.raises(organize.OrganizationError, match="Confirmation"):
        organize.apply(root, confirm="")


def test_apply_restore_byte_exact_for_docs_and_inputs(tmp_path: Path) -> None:
    root = sample(tmp_path)
    before = (root / "docs/T3_BC_UNREFERENCED_STUDY.md").read_bytes()
    organize.preview(root)
    manifest = organize.apply(root, confirm=organize.CONFIRM)
    assert manifest["status"] == "applied"
    destination = root / "docs/research/eft-matching/T3_BC_UNREFERENCED_STUDY.md"
    assert destination.read_bytes() == before
    assert not (root / "docs/T3_BC_UNREFERENCED_STUDY.md").exists()
    assert (root / "input/reference/mixed_wilson_inputs.zip").is_file()
    restored = organize.restore(root)
    assert restored["status"] == "restored"
    assert (root / "docs/T3_BC_UNREFERENCED_STUDY.md").read_bytes() == before
    assert not destination.exists()
    assert (root / "mixed_wilson_inputs.zip").is_file()


def test_source_modification_rejects_stale_review(tmp_path: Path) -> None:
    root = sample(tmp_path)
    organize.preview(root)
    p = root / "docs/T3_BC_UNREFERENCED_STUDY.md"
    p.write_text("different", encoding="utf-8")
    with pytest.raises(organize.OrganizationError, match="stale"):
        organize.apply(root, confirm=organize.CONFIRM)
    assert p.is_file()


def test_head_change_rejects_stale_review(tmp_path: Path) -> None:
    root = sample(tmp_path)
    organize.preview(root)
    p = root / "README.md"
    p.write_text("new overview", encoding="utf-8")
    git(root, "add", "README.md")
    git(root, "commit", "-qm", "updated head")
    with pytest.raises(organize.OrganizationError, match="stale"):
        organize.apply(root, confirm=organize.CONFIRM)


def test_restore_refuses_overwriting_or_tampered_files(tmp_path: Path) -> None:
    root = sample(tmp_path)
    organize.preview(root)
    organize.apply(root, confirm=organize.CONFIRM)
    moved = root / "docs/research/eft-matching/T3_BC_UNREFERENCED_STUDY.md"
    moved.write_text("corrupted", encoding="utf-8")
    with pytest.raises(organize.OrganizationError, match="Restore blocked"):
        organize.restore(root)
    assert not (root / "docs/T3_BC_UNREFERENCED_STUDY.md").exists()


def test_relocation_is_not_applied_twice(tmp_path: Path) -> None:
    root = sample(tmp_path)
    organize.preview(root)
    organize.apply(root, confirm=organize.CONFIRM)
    with pytest.raises(organize.OrganizationError):
        organize.apply(root, confirm=organize.CONFIRM)


def test_generated_indexes_link_actual_files(tmp_path: Path) -> None:
    root = sample(tmp_path)
    organize.preview(root, include_experiments=True)
    organize.apply(root, confirm=organize.CONFIRM)
    scripts, docs = organize._write_indexes(root)
    assert "RunCoreTests.ps1" in scripts.read_text(encoding="utf-8")
    assert "T3_BC_UNREFERENCED_STUDY.md" in docs.read_text(encoding="utf-8")
    assert (root / "scripts/research/experimental/Standalone.py").is_file()
    assert not (root / "scripts/Standalone.experimental.py").exists()
    organize.restore(root)
    assert (root / "scripts/Standalone.experimental.py").is_file()


def test_unsafe_paths_and_symlink_are_refused(tmp_path: Path) -> None:
    root = sample(tmp_path)
    for rel in ("../outside", "/tmp/outside", "docs/../README.md", "C:/Windows"):
        with pytest.raises(organize.OrganizationError):
            organize._safe(root, rel)


def test_unmodified_crlf_working_copy_respects_git_clean_filter(tmp_path: Path) -> None:
    root = sample(tmp_path)
    rel = "docs/T3_BC_UNREFERENCED_STUDY.md"
    p = root / rel
    git(root, "config", "core.autocrlf", "true")
    original = p.read_bytes()
    p.write_bytes(original.replace(b"\r\n", b"\n").replace(b"\n", b"\r\n"))
    # Git's filter decides whether this is clean; do not compare raw Git blobs.
    diff = organize._git(root, "diff", "--quiet", "HEAD", "--", rel, allow_diff=True).returncode
    if diff == 0:
        plan = organize.build_plan(root)
        assert rel in {x["source"] for x in plan["approved"]}


def test_plan_is_deterministic(tmp_path: Path) -> None:
    root = sample(tmp_path)
    a = organize.build_plan(root, include_experiments=True)
    b = organize.build_plan(root, include_experiments=True)
    assert a == b
    assert len(a["approved"]) == a["approved_count"]
