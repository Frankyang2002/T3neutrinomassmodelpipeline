"""Native package-source evidence is text only, not a physical phase certificate.

The archived diagnostic can be written by Windows PowerShell 5.1 as UTF-16LE.
These tests reproduce the exact missing-source-marker bug seen on Windows.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path
import tempfile
import unittest

from Numerical.lo_certification.spinor.boundary import native_evidence_probe


ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = (
    'T3_MATCHeTE_GAMMACC_DIAGNOSTIC_V2\r\n'
    'GammaCC::usage = "GammaCC is the charge conjugation matrix"\r\n'
    'PackageExport["GammaCC"]\r\n'
    'CConj[Field[f, Fermion, _, _]]\r\n'
    'DiracProduct[GammaCC]\r\n'
    'Proj[+/-1] is the chiral right/left projectors\r\n'
    'DiracProduct[GammaCC, Proj[-1]]\r\n'
)


class NativeEvidenceTests(unittest.TestCase):
    def probe(self, payload):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / 'native_probe.txt'
            p.write_bytes(payload)
            result = native_evidence_probe(p)
        return result

    def assert_found_but_not_physical(self, result):
        self.assertEqual(result['status'],
                         'PACKAGE_SOURCE_LOCATED_ONLY_NOT_NUMERICAL_MATRIX_PHASE')
        self.assertTrue(all(result['native_items'].values()))
        self.assertTrue(result['Matchete_pl_is_left_projector'])
        self.assertFalse(result['GammaCC_Phase_certified'])
        self.assertIsNone(result['GammaCC_to_master_weyl_factor'])

    def test_native_utf8(self):
        result = self.probe(EVIDENCE.encode('utf-8'))
        self.assert_found_but_not_physical(result)
        self.assertEqual(result['detected_text_encoding'], 'utf-8')

    def test_native_utf8_with_bom(self):
        result = self.probe(EVIDENCE.encode('utf-8-sig'))
        self.assert_found_but_not_physical(result)
        self.assertEqual(result['detected_text_encoding'], 'utf-8-sig')

    def test_windows_powershell_utf16_with_bom(self):
        result = self.probe(EVIDENCE.encode('utf-16'))
        self.assert_found_but_not_physical(result)
        self.assertEqual(result['detected_text_encoding'], 'utf-16')

    def test_utf16_le_and_be_without_bom(self):
        for codec in ('utf-16-le', 'utf-16-be'):
            with self.subTest(codec=codec):
                result = self.probe(EVIDENCE.encode(codec))
                self.assert_found_but_not_physical(result)
                self.assertEqual(result['detected_text_encoding'], codec)

    def test_undecodable_source_stays_uncertified(self):
        result = self.probe(b'\xff\x00\x91\x80\xff')
        self.assertEqual(result['status'], 'UNDECODABLE_PACKAGE_SOURCE')
        self.assertFalse(result['GammaCC_Phase_certified'])
        self.assertIsNone(result['GammaCC_to_master_weyl_factor'])

    def test_marker_only_does_not_become_physical_certificate(self):
        result = self.probe(b'GammaCC::usage = "charge conjugation"')
        self.assertFalse(result['GammaCC_Phase_certified'])
        self.assertIsNone(result['GammaCC_to_master_weyl_factor'])

    def test_new_evidence_path_does_not_reuse_legacy_filename(self):
        # Imports the actual current runner without calling main().
        spec = importlib.util.spec_from_file_location(
            't3_spinor_runner_for_evidence_test',
            ROOT/'scripts/RunLOSpinorCertification.py',
        )
        self.assertIsNotNone(spec)
        runner = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(runner)
        h = 'a' * 64
        self.assertEqual(runner.evidence_report_filename({'source_sha256': h}),
                         'spinor_evidence_v2_aaaaaaaaaaaaaaaa.json')
        self.assertEqual(runner.evidence_report_filename({}),
                         'spinor_evidence_v2_not_available.json')
        self.assertNotEqual(runner.evidence_report_filename({'source_sha256': h}),
                            'spinor_evidence_and_status.json')


if __name__ == '__main__':
    unittest.main()
