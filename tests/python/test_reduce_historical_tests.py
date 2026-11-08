"""Tests for the reversible, Git-pinned historical-test retirement tool."""
from __future__ import annotations

import importlib.util
from pathlib import Path
import subprocess
import tempfile

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "scripts/maintenance/ReduceHistoricalTests.py"
SPEC = importlib.util.spec_from_file_location("reduce_historical_tests", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
cleanup = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(cleanup)


def git(root: Path, *args: str) -> str:
    completed = subprocess.run(["git", *args], cwd=root, capture_output=True, text=True,
                               check=True)
    return completed.stdout.strip()


def repository(base: Path, *, referenced: bool = False, import_ref: bool = False) -> Path:
    root = base / "T3"
    root.mkdir()
    git(root, "init", "-q")
    git(root, "config", "user.email", "test@example.invalid")
    git(root, "config", "user.name", "Test")
    files = {
        "tests/test_t3_frozen_persistent.py": "def test_old():\n    assert True\n",
        "tests/test_pipeline_plan.py": "def test_new():\n    assert True\n",
        "tests/python/test_t3_general_pdf_functionality.py": "def test_pdf():\n    assert True\n",
        "scripts/driver.py": (
            "# requires tests/test_t3_frozen_persistent.py\n"
            if referenced else "from tests.test_t3_frozen_persistent import test_old\n"
            if import_ref else "print('ok')\n"
        ),
    }
    for rel, content in files.items():
        path = root / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
    git(root, "add", ".")
    git(root, "commit", "-q", "-m", "fixture")
    return root


def test_name_selection_preserves_physics_and_reference_data():
    assert cleanup.test_group("tests/test_t3_frozen_persistent.py") == "historical-benchmark"
    assert cleanup.test_group("tests/python/test_t3_general_pdf_functionality.py") == "missing-pdf-writer"
    assert cleanup.test_group("tests/test_pipeline_plan.py") is None
    assert cleanup.test_group("tests/test_batch34_gauge_covariance.py") is None
    assert cleanup.test_group("tests/python/test_general_d_m2_scalar_threshold.py") is None
    assert cleanup.test_group("tests/reference/MaUVRGE.py") is None


def test_dry_run_collects_and_preserves_files():
    with tempfile.TemporaryDirectory() as tmp:
        root = repository(Path(tmp))
        plan = cleanup.build_plan(root, expected_head=None)
        assert plan["original_python_test_modules"] == 3
        assert plan["selected_by_name"] == 2
        assert len(plan["approved"]) == 2
        assert len(plan["held"]) == 0
        assert (root / "tests/test_t3_frozen_persistent.py").exists()
        assert cleanup.save_plan(root, plan).is_file()


def test_reference_blocks_archival():
    with tempfile.TemporaryDirectory() as tmp:
        root = repository(Path(tmp), referenced=True)
        plan = cleanup.build_plan(root, expected_head=None)
        approved = [x["path"] for x in plan["approved"]]
        held = {x["path"]: x for x in plan["held"]}
        assert approved == ["tests/python/test_t3_general_pdf_functionality.py"]
        old = held["tests/test_t3_frozen_persistent.py"]
        assert "scripts/driver.py" in old["references"]


def test_modified_test_is_held():
    with tempfile.TemporaryDirectory() as tmp:
        root = repository(Path(tmp))
        (root / "tests/test_t3_frozen_persistent.py").write_text("# local changes\n", encoding="utf-8")
        plan = cleanup.build_plan(root, expected_head=None)
        held = {x["path"]: x for x in plan["held"]}
        assert "different from Git HEAD" in held["tests/test_t3_frozen_persistent.py"]["hold_reason"]


def test_archive_restore_and_prevent_overwrite():
    with tempfile.TemporaryDirectory() as tmp:
        root = repository(Path(tmp))
        plan = cleanup.build_plan(root, expected_head=None)
        archive = cleanup.apply_plan(root, plan)
        assert archive is not None
        assert not (root / "tests/test_t3_frozen_persistent.py").exists()
        assert (root / "tests/test_pipeline_plan.py").is_file()
        assert (archive / "manifest.json").is_file()
        (root / "tests/test_t3_frozen_persistent.py").write_text("new", encoding="utf-8")
        with pytest.raises(FileExistsError):
            cleanup.restore_archive(root, archive)
        (root / "tests/test_t3_frozen_persistent.py").unlink()
        assert cleanup.restore_archive(root, archive) == 2
        assert (root / "tests/python/test_t3_general_pdf_functionality.py").is_file()


def test_archive_modified_after_plan_rejected():
    with tempfile.TemporaryDirectory() as tmp:
        root = repository(Path(tmp))
        plan = cleanup.build_plan(root, expected_head=None)
        (root / "tests/test_t3_frozen_persistent.py").write_text("changed", encoding="utf-8")
        with pytest.raises(RuntimeError, match="changed after preview"):
            cleanup.apply_plan(root, plan)
        assert (root / "tests/test_t3_frozen_persistent.py").exists()


def test_unknown_commit_rejected():
    with tempfile.TemporaryDirectory() as tmp:
        root = repository(Path(tmp))
        with pytest.raises(RuntimeError, match="Git HEAD changed"):
            cleanup.build_plan(root, expected_head="0" * 40)


def test_parent_git_checkout_is_not_accepted():
    with tempfile.TemporaryDirectory() as tmp:
        base = Path(tmp)
        parent = base / "parent"
        parent.mkdir()
        git(parent, "init", "-q")
        nested = parent / "T3"
        nested.mkdir()
        with pytest.raises(RuntimeError, match="Git root"):
            cleanup.checked_repository(nested, expected_head=None)


def test_module_import_without_py_suffix_blocks_archival():
    with tempfile.TemporaryDirectory() as tmp:
        root = repository(Path(tmp), import_ref=True)
        plan = cleanup.build_plan(root, expected_head=None)
        held = {x["path"]: x for x in plan["held"]}
        assert held["tests/test_t3_frozen_persistent.py"]["references"] == ["scripts/driver.py"]


def test_windows_crlf_clean_filter_is_not_misclassified_as_dirty():
    with tempfile.TemporaryDirectory() as tmp:
        root = repository(Path(tmp))
        git(root, "config", "core.autocrlf", "true")
        path = root / "tests/test_t3_frozen_persistent.py"
        path.write_bytes(path.read_bytes().replace(b"\n", b"\r\n"))
        plan = cleanup.build_plan(root, expected_head=None)
        assert "tests/test_t3_frozen_persistent.py" in {x["path"] for x in plan["approved"]}


def test_preview_binds_to_current_head_and_apply_requires_review():
    with tempfile.TemporaryDirectory() as tmp:
        root = repository(Path(tmp))
        assert cleanup.build_plan(root)["git_head"] == git(root, "rev-parse", "HEAD")
        with pytest.raises(RuntimeError, match="No reviewed plan"):
            cleanup.load_reviewed_plan(root)
        cleanup.save_plan(root, cleanup.build_plan(root))
        assert len(cleanup.load_reviewed_plan(root)["approved"]) == 2
        # A new commit invalidates the previously reviewed selection.
        (root / "README.txt").write_text("new commit", encoding="utf-8")
        git(root, "add", "README.txt")
        git(root, "commit", "-q", "-m", "later")
        with pytest.raises(RuntimeError, match="Git HEAD changed"):
            cleanup.apply_plan(root, cleanup.load_reviewed_plan(root))
        assert (root / "tests/test_t3_frozen_persistent.py").exists()


def test_exact_lf_conversion_is_only_allowed_for_managed_text():
    with tempfile.TemporaryDirectory() as tmp:
        root = repository(Path(tmp))
        rel = "tests/test_t3_frozen_persistent.py"
        git(root, "config", "core.autocrlf", "true")
        path = root / rel
        path.write_bytes(path.read_bytes().replace(b"\n", b"\r\n"))
        head_blob = cleanup.tracked_head_entries(root)[rel][1]
        assert cleanup.same_as_head_blob(root, rel, head_blob)
        path.write_bytes(path.read_bytes() + b"# substantive change\r\n")
        assert not cleanup.same_as_head_blob(root, rel, head_blob)


def test_main_dry_run_and_cli_apply_use_reviewed_plan():
    with tempfile.TemporaryDirectory() as tmp:
        root = repository(Path(tmp))
        assert cleanup.main(["--repo", str(root)]) == 0
        plan = cleanup.load_reviewed_plan(root)
        assert len(plan["approved"]) == 2
        assert cleanup.main(["--repo", str(root), "--apply",
                             "--confirm", cleanup.CONFIRM]) == 0
        assert not (root / "tests/test_t3_frozen_persistent.py").exists()


def test_crlf_fallback_handles_windows_hash_filter_mismatch(monkeypatch):
    with tempfile.TemporaryDirectory() as tmp:
        root = repository(Path(tmp))
        rel = "tests/test_t3_frozen_persistent.py"
        git(root, "config", "core.autocrlf", "true")
        path = root / rel
        path.write_bytes(path.read_bytes().replace(b"\n", b"\r\n"))
        head_blob = cleanup.tracked_head_entries(root)[rel][1]
        monkeypatch.setattr(cleanup, "git_clean_blob_sha1", lambda root, rel: "wrong")
        assert cleanup.same_as_head_blob(root, rel, head_blob)
        path.write_bytes(path.read_bytes() + b"# modified\r\n")
        assert not cleanup.same_as_head_blob(root, rel, head_blob)
