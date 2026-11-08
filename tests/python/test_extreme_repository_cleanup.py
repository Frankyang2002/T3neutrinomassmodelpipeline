"""Synthetic integration checks for the standalone reversible cleanup tool."""
from __future__ import annotations

import importlib.util
import shutil
import subprocess
from pathlib import Path
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "maintenance" / "ExtremeRepositoryCleanup.py"
spec = importlib.util.spec_from_file_location("extreme_cleanup", SCRIPT)
assert spec and spec.loader
cleanup = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cleanup)


class CleanupTests(unittest.TestCase):
    def create_root(self, where: Path) -> Path:
        root = where / "T3Repo"
        (root / "Lagrangian").mkdir(parents=True)
        (root / "docs").mkdir()
        (root / "scripts").mkdir()
        (root / "Reports" / "output").mkdir(parents=True)
        (root / "pipeline.py").write_text("# production driver\n", encoding="utf-8")
        (root / "README.md").write_text("# new overview\n", encoding="utf-8")
        (root / "README_OLD_BATCH.md").write_text("historical material\n", encoding="utf-8")
        (root / "docs" / "OLD_TECHNICAL_NOTE.md").write_text("old\n", encoding="utf-8")
        (root / "docs" / "MODEL_AND_PIPELINE.md").write_text("technical overview\n", encoding="utf-8")
        (root / "docs" / "REFERENCED_VALIDATION.md").write_text("keep it\n", encoding="utf-8")
        (root / "scripts" / "check.py").write_text("# docs/REFERENCED_VALIDATION.md\n", encoding="utf-8")
        (root / "Lagrangian" / "core.wl").write_text("T3_EXPERIMENTAL_GENERAL_FERMIONS\n", encoding="utf-8")
        (root / "Reports" / "output" / "figure.png").write_bytes(b"not-important-test-image")
        (root / "Reports" / "output" / "source.json").write_text('{"important": true}', encoding="utf-8")
        return root

    def test_dry_run_and_references(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.create_root(Path(tmp))
            plan = cleanup.build_plan(root)
            names = {a["path"] for a in plan["archive_candidates"]}
            held = {a["path"] for a in plan["held_for_references"]}
            self.assertIn("README_OLD_BATCH.md", names)
            self.assertIn("docs/OLD_TECHNICAL_NOTE.md", names)
            self.assertIn("docs/REFERENCED_VALIDATION.md", held)
            self.assertNotIn("docs/MODEL_AND_PIPELINE.md", names)
            self.assertNotIn("Reports/output/source.json", names)
            self.assertNotIn("Reports/output/figure.png", names)
            self.assertTrue(plan["general_source_flags"])
            # Planning must not remove files.
            self.assertTrue((root / "README_OLD_BATCH.md").is_file())

    def test_new_reference_blocks_retirement(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.create_root(Path(tmp))
            # Simulate a newly installed authoritative doc not yet git-added.
            (root / "docs" / "MODEL_AND_PIPELINE.md").write_text(
                "See docs/OLD_TECHNICAL_NOTE.md for the original derivation\n", encoding="utf-8"
            )
            plan = cleanup.build_plan(root)
            held = {a["path"] for a in plan["held_for_references"]}
            self.assertIn("docs/OLD_TECHNICAL_NOTE.md", held)

    def test_archive_and_restore(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.create_root(Path(tmp))
            plan = cleanup.build_plan(root)
            archive = cleanup.apply_plan(root, plan)
            self.assertIsNotNone(archive)
            self.assertFalse((root / "README_OLD_BATCH.md").exists())
            self.assertTrue((root / "docs" / "REFERENCED_VALIDATION.md").exists())
            count = cleanup.restore(root, archive)
            self.assertEqual(count, len(plan["archive_candidates"]))
            self.assertEqual((root / "README_OLD_BATCH.md").read_text(), "historical material\n")
            with self.assertRaises(FileExistsError):
                cleanup.restore(root, archive)

    def test_nested_untracked_fixture_inside_other_git_worktree(self):
        """Regression: Git discovery must not treat an ancestor repo as fixture."""
        if not shutil.which("git"):
            self.skipTest("git executable unavailable")
        with tempfile.TemporaryDirectory() as tmp:
            parent = Path(tmp) / "outer_checkout"
            parent.mkdir()
            subprocess.run(["git", "init", "-q", str(parent)], check=True)
            root = self.create_root(parent)
            self.assertFalse(cleanup.is_git_checkout_root(root))
            plan = cleanup.build_plan(root)
            archivable = {item["path"] for item in plan["archive_candidates"]}
            held = {item["path"] for item in plan["held_for_references"]}
            self.assertIn("README_OLD_BATCH.md", archivable)
            self.assertIn("docs/OLD_TECHNICAL_NOTE.md", archivable)
            self.assertIn("docs/REFERENCED_VALIDATION.md", held)
            self.assertEqual(plan["git_head"], "")

    def test_git_checkout_ignores_untracked_archive_candidates(self):
        """Only tracked documents are archive candidates in a real checkout."""
        if not shutil.which("git"):
            self.skipTest("git executable unavailable")
        with tempfile.TemporaryDirectory() as tmp:
            root = self.create_root(Path(tmp))
            subprocess.run(["git", "init", "-q", str(root)], check=True)
            subprocess.run(["git", "-C", str(root), "add", "README_OLD_BATCH.md"], check=True)
            subprocess.run([
                "git", "-C", str(root), "-c", "user.name=Test",
                "-c", "user.email=test@example.com", "commit", "-qm", "fixture"
            ], check=True)
            self.assertTrue(cleanup.is_git_checkout_root(root))
            plan = cleanup.build_plan(root)
            archivable = {item["path"] for item in plan["archive_candidates"]}
            self.assertIn("README_OLD_BATCH.md", archivable)
            self.assertNotIn("docs/OLD_TECHNICAL_NOTE.md", archivable)
            self.assertNotIn("Reports/output/figure.png", archivable)
            # New untracked code must still block retiring referenced tracked docs.
            (root / "scripts" / "new_check.py").write_text(
                '"README_OLD_BATCH.md"\n', encoding="utf-8"
            )
            plan = cleanup.build_plan(root)
            held = {item["path"] for item in plan["held_for_references"]}
            self.assertIn("README_OLD_BATCH.md", held)

    def test_modified_tracked_candidate_held(self):
        if not shutil.which("git"):
            self.skipTest("git executable unavailable")
        with tempfile.TemporaryDirectory() as tmp:
            root = self.create_root(Path(tmp))
            subprocess.run(["git", "init", "-q", str(root)], check=True)
            subprocess.run(["git", "-C", str(root), "add", "README_OLD_BATCH.md"], check=True)
            subprocess.run([
                "git", "-C", str(root), "-c", "user.name=Test",
                "-c", "user.email=test@example.com", "commit", "-qm", "fixture"
            ], check=True)
            (root / "README_OLD_BATCH.md").write_text("locally edited\n", encoding="utf-8")
            plan = cleanup.build_plan(root)
            held = {item["path"]: item for item in plan["held_for_references"]}
            self.assertIn("README_OLD_BATCH.md", held)
            self.assertEqual(held["README_OLD_BATCH.md"]["reason"], "locally modified tracked file")

    def test_unreadable_checkout_metadata_fails_closed(self):
        """A damaged real checkout must never be mistaken for a synthetic fixture."""
        if not shutil.which("git"):
            self.skipTest("git executable unavailable")
        with tempfile.TemporaryDirectory() as tmp:
            root = self.create_root(Path(tmp))
            (root / ".git").write_text("gitdir: missing-git-directory\n", encoding="utf-8")
            with self.assertRaises(RuntimeError):
                cleanup.build_plan(root)

    def test_presentations_are_opt_in(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.create_root(Path(tmp))
            plan = cleanup.build_plan(root, presentations=True)
            names = {a["path"] for a in plan["archive_candidates"]}
            self.assertIn("Reports/output/figure.png", names)
            self.assertNotIn("Reports/output/source.json", names)


if __name__ == "__main__":
    unittest.main()
