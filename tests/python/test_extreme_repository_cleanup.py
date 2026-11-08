"""Synthetic integration checks for the standalone reversible cleanup tool."""
from __future__ import annotations

import importlib.util
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

    def test_presentations_are_opt_in(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.create_root(Path(tmp))
            plan = cleanup.build_plan(root, presentations=True)
            names = {a["path"] for a in plan["archive_candidates"]}
            self.assertIn("Reports/output/figure.png", names)
            self.assertNotIn("Reports/output/source.json", names)


if __name__ == "__main__":
    unittest.main()
