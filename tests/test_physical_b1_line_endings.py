"""B1 Windows source guard: line ending equivalence without content drift."""
from __future__ import annotations

import hashlib
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.RunPhysicalB1WeylAudit import guard_source, SOURCE_GIT_BLOBS, gitblob
from Numerical.lo_certification.physical_normalisation.b1_matchete_weyl import B1Error

REL = 'Lagrangian/model/SU2Invariants.wl'
LF = b'DefineCG[x_] := x;\nContract[a_, b_] := a + b;\n'


def independent_blob(data: bytes) -> str:
    return hashlib.sha1(b'blob ' + str(len(data)).encode('ascii') + b'\0' + data).hexdigest()


class GuardLineEndingTest(unittest.TestCase):
    def check_fixture(self, data: bytes, *, should_pass: bool):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            p = root / REL
            p.parent.mkdir(parents=True)
            p.write_bytes(data)
            with patch.dict(SOURCE_GIT_BLOBS, {REL: independent_blob(LF)}, clear=True):
                if should_pass:
                    record = guard_source(root)[REL]
                    self.assertTrue(record['github_main_blob_match'])
                    self.assertTrue(record['source_bytes_match_after_line_ending_normalization'])
                    self.assertEqual(record['git_blob_sha1'], independent_blob(LF))
                    self.assertEqual(record['working_tree_raw_git_blob_sha1'], independent_blob(data))
                    self.assertEqual(record['crlf_normalization_applied'], data != LF)
                else:
                    with self.assertRaisesRegex(B1Error, 'SOURCE DRIFT'):
                        guard_source(root)

    def test_lf_exact_source_accepted(self):
        self.check_fixture(LF, should_pass=True)

    def test_crlf_checkout_of_same_source_accepted(self):
        self.check_fixture(LF.replace(b'\n', b'\r\n'), should_pass=True)

    def test_modified_operator_rejected_even_with_crlf(self):
        self.check_fixture(LF.replace(b'a + b', b'a - b').replace(b'\n', b'\r\n'),
                           should_pass=False)

    def test_modified_operator_rejected_with_lf(self):
        self.check_fixture(LF.replace(b'x;', b'2*x;'), should_pass=False)

    def test_bare_carriage_return_rejected(self):
        self.check_fixture(LF.replace(b'\n', b'\r', 1), should_pass=False)

    def test_missing_source_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            with patch.dict(SOURCE_GIT_BLOBS, {REL: independent_blob(LF)}, clear=True):
                with self.assertRaisesRegex(B1Error, 'Missing local source'):
                    guard_source(Path(temp))

    def test_binary_blob_implementation_has_git_header(self):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp) / 'a.wl'
            p.write_bytes(LF)
            self.assertEqual(gitblob(p), independent_blob(LF))


if __name__ == '__main__':
    unittest.main()
