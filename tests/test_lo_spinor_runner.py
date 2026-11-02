"""Fail-closed batch workflow and refusal to fake local repository data."""
import subprocess
import tempfile
import sys
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]


class RunnerTests(unittest.TestCase):
    def test_missing_upstream_release_blocks(self):
        """Force an absent release, independently of a real checkout's files.

        The default source release exists on the user's Windows checkout, so
        assuming it is missing makes this regression test falsely fail there.
        """
        with tempfile.TemporaryDirectory() as td:
            missing_release = Path(td) / 'intentionally_missing_release'
            self.assertFalse(missing_release.exists())
            proc = subprocess.run(
                [sys.executable, str(ROOT/'scripts/RunLOSpinorCertification.py'),
                 '--repo-root', str(ROOT),
                 '--source-release', str(missing_release), '--check-only'],
                capture_output=True, text=True,
            )
        self.assertEqual(proc.returncode, 2, proc.stdout + proc.stderr)
        self.assertIn('Missing source:', proc.stderr)
        self.assertIn('summary.json', proc.stderr)

    def test_output_restricted_to_owned_path(self):
        proc = subprocess.run([sys.executable, str(ROOT/'scripts/RunLOSpinorCertification.py'),
                               '--repo-root', str(ROOT), '--output', str(ROOT/'RGE/')],
                              capture_output=True, text=True)
        self.assertEqual(proc.returncode, 2)
        self.assertIn('Output must stay inside', proc.stderr)


if __name__ == '__main__':
    unittest.main()
