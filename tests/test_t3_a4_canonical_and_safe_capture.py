"""Static EFT safeguards and encoding-capture regression (not physics validation)."""
from pathlib import Path
import ast
import unittest

ROOT = Path(__file__).resolve().parents[1]
A = ROOT / 'Lagrangian/RunT3AMinus4FullGeneralPhysicalScalar.wl'
SAFE = ROOT / 'scripts/RunT3ThesisFinalOneBatchSafe.py'


class TestA4EFTSafeguards(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.a = A.read_text(encoding='utf-8')

    def test_original_mass_eft_order_promoted_reversibly(self):
        for name in ('MS1', 'MS2'):
            self.assertIn(f'Coupling[{name}, {{}}, 1]] :> Coupling[{name}, {{}}, 0]', self.a)
        self.assertIn('t3ADemoted === t3ATree', self.a)

    def test_canonicalisation_is_strictly_alpha_equivalent(self):
        self.assertIn('t3ARelabelEquivalent', self.a)
        self.assertIn('t3AAlphaKeys[t3ARelabelCandidate] === t3AAlphaKeys[t3AMatchInput]', self.a)
        self.assertIn('t3ASMAlphaEquivalent', self.a)
        self.assertIn('t3ASMPreserved', self.a)
        self.assertIn('t3ACheckSMRebuild = t3ACheckSafe[t3ASMRebuilt]', self.a)

    def test_diagnostic_projection_is_not_matched(self):
        self.assertIn('t3ACheckWithoutW', self.a)
        self.assertIn('t3ACheckWithoutWAndMix', self.a)
        self.assertIn('"projection_used_for_physical_match" -> False', self.a)
        self.assertIn('t3AMatchInput = t3ASelectedInput;', self.a)
        self.assertEqual(self.a.count('Match[t3AMatchInput, EFTOrder -> 5, LoopOrder -> 1]'), 1)

    def test_full_one_loop_wilson_and_finite_status_not_spoofed(self):
        self.assertIn('"EFT1_Wilson_RG_included" -> False', self.a)
        self.assertIn('"renormalized_physical_C5_verified" -> False', self.a)
        self.assertIn('A_PhysicalScalarSoftOnly_HardInheritanceAndWilsonRGPending', self.a)
        self.assertIn('Exit[5]', self.a)
        self.assertNotIn('NF=1', self.a)


class TestSafeSubprocessCapture(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = SAFE.read_text(encoding='utf-8')
        ast.parse(cls.source)

    def test_binary_read_avoids_python_readerthread(self):
        self.assertIn('errors="backslashreplace"', self.source)
        self.assertIn('kwargs.pop("text", False)', self.source)
        self.assertIn('result = _native_subprocess.run(*args, **kwargs)', self.source)

    def test_original_runner_is_not_replaced_or_reimplemented(self):
        self.assertIn('from scripts import RunT3ThesisFinalOneBatch as _original', self.source)
        self.assertIn('_original.main(argv)', self.source)
        self.assertIn('_original.subprocess = _BinaryCaptureAdapter', self.source)


if __name__ == '__main__':
    unittest.main()
