import tempfile
import unittest
from pathlib import Path
from scripts.RunT3WorkerCRestrictedDMinus2 import preflight, REQUIRED

class TestRestrictedDMinus2Preflight(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.base = Path(self.temp.name)
        self.stage = self.base / "stage"
        self.profile = self.stage / "data/general_d_m2_scalar_threshold/profile_tree"
        self.profile.mkdir(parents=True)
        self.script = self.base / "matching.wl"
        self.script.write_text("(* stub for input checks only *)")
        self.out = self.base / "audits/new_run"
    def populate(self):
        for name in REQUIRED:
            (self.profile / name).write_bytes(b"fixture")
    def test_complete_inputs(self):
        self.populate()
        preflight(self.stage,self.out,self.script)
    def test_missing_input_blocks(self):
        self.populate()
        (self.profile / REQUIRED[0]).unlink()
        with self.assertRaisesRegex(ValueError,"missing"):
            preflight(self.stage,self.out,self.script)
    def test_no_overwrite(self):
        self.populate()
        self.out.mkdir(parents=True)
        with self.assertRaisesRegex(ValueError,"overwrite"):
            preflight(self.stage,self.out,self.script)
    def test_no_output_in_source(self):
        self.populate()
        with self.assertRaisesRegex(ValueError,"preserved"):
            preflight(self.stage,self.stage / "output/new",self.script)
    def test_no_missing_script(self):
        self.populate()
        self.script.unlink()
        with self.assertRaisesRegex(ValueError,"implementation"):
            preflight(self.stage,self.out,self.script)

if __name__ == '__main__': unittest.main()

class TestRestrictedDMinus2Runner(unittest.TestCase):
    def test_association_transforms_use_map(self):
        from scripts.RunT3WorkerCRestrictedDMinus2 import SCRIPT
        content = SCRIPT.read_text(encoding='utf-8')
        self.assertNotIn('AssociationMap[', content)
        for fragment in ('inputs=Map[', 'restricted=Map[', 'parts=Map[', 'sourceSha = Map['):
            self.assertIn(fragment, content)

    def test_wolfram_zero_without_report_fails(self):
        from unittest.mock import patch
        from scripts import RunT3WorkerCRestrictedDMinus2 as runner
        from types import SimpleNamespace
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            profile = root / 'data/general_d_m2_scalar_threshold/profile_tree'
            profile.mkdir(parents=True)
            for item in REQUIRED:
                (profile / item).write_bytes(b'fixture')
            out = root.parent / (root.name + '_unused_output')
            with patch.object(runner, 'DEFAULT', root), patch.object(runner, 'SCRIPT', runner.SCRIPT), \
                 patch('sys.argv', ['runner', '--output-dir', str(out)]), \
                 patch.object(runner.subprocess, 'run', return_value=SimpleNamespace(returncode=0)):
                self.assertEqual(runner.main(), 4)

    def test_wolfram_nonzero_kept(self):
        from unittest.mock import patch
        from scripts import RunT3WorkerCRestrictedDMinus2 as runner
        from types import SimpleNamespace
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            profile = root / 'data/general_d_m2_scalar_threshold/profile_tree'
            profile.mkdir(parents=True)
            for item in REQUIRED:
                (profile / item).write_bytes(b'fixture')
            out = root.parent / (root.name + '_unused_output')
            with patch.object(runner, 'DEFAULT', root), patch('sys.argv', ['runner', '--output-dir', str(out)]), \
                 patch.object(runner.subprocess, 'run', return_value=SimpleNamespace(returncode=8)):
                self.assertEqual(runner.main(), 8)
