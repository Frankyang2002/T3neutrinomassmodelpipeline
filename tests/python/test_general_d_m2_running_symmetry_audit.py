"""Contract tests for the read-only direct-running symmetry diagnostic."""
from __future__ import annotations

import ast
import json
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
WOLFRAM = ROOT / "Lagrangian/RunGeneralDMinus2RunningSymmetryAudit.wl"
RUNNER = ROOT / "scripts/RunGeneralDMinus2RunningSymmetryAudit.py"


class RunningSymmetryContract(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.wl = WOLFRAM.read_text(encoding="utf-8")
        cls.py = RUNNER.read_text(encoding="utf-8")

    def test_wolfram_exists(self):
        self.assertTrue(WOLFRAM.is_file())

    def test_python_compiles(self):
        ast.parse(self.py)

    def test_no_matching_invocation(self):
        self.assertNotIn("Match[", self.wl.replace("Match[]", ""))

    def test_derives_p_and_q(self):
        self.assertIn("EFT1D2WilsonRunningWeinbergCoefficient[p, q]", self.wl)
        self.assertIn("EFT1D2WilsonRunningWeinbergCoefficient[q, p]", self.wl)

    def test_does_not_automatically_symmetrize(self):
        self.assertIn('"RunningNormalizationCertified" -> False', self.wl)
        self.assertIn('"FixedOrderSMEFTCertified" -> False', self.wl)
        self.assertIn('"antisymmetric_part_NOT_APPLIED"', self.wl)

    def test_exports_generator_definition(self):
        self.assertIn("DownValues[EFT1D2WilsonRunningWeinbergCoefficient]", self.wl)

    def test_exact_swap_comparison(self):
        self.assertIn("residual = Expand[runPQ - runQP]", self.wl)

    def test_simplify_bounded(self):
        self.assertIn("TimeConstrained[", self.wl)

    def test_no_overwrite_saved_efts(self):
        self.assertNotIn('"eft1_after_F_realbasis_tree.wxf"', self.wl)
        self.assertNotIn('"scalar_hard_smeft.wxf"', self.wl)

    def test_error_handling(self):
        self.assertIn("read_report", self.py)
        self.assertIn("InvalidWolframResultJSON", self.py)

    def test_preflight_uses_existing_runner(self):
        self.assertIn("preflight(matching_dir, script=SCRIPT)", self.py)

    def test_bound_sum_alpha_normalization_is_inert(self):
        self.assertIn("HoldPattern[Sum[body_, {r_Symbol, lo_, hi_}]]", self.wl)
        self.assertIn("Inactive[Sum][body /. r -> D2AlphaDummy", self.wl)
        self.assertIn("alphaResidual = Expand[alphaPQ - alphaQP]", self.wl)

    def test_all_sum_binders_must_be_converted(self):
        self.assertIn("originalSums =", self.wl)
        self.assertIn("normalizedSums =", self.wl)
        self.assertIn("originalSums =!= normalizedSums", self.wl)
        self.assertIn("D2ActiveSums[alphaPQ] =!= 0", self.wl)

    def test_does_not_rewrite_running_source(self):
        self.assertIn('"RunningFlavourSymmetryCertifiedModuloDummyRenaming" -> alphaEqual', self.wl)
        self.assertIn('"RunningNormalizationCertified" -> False', self.wl)
        self.assertIn('"FixedOrderSMEFTCertified" -> False', self.wl)

    def test_term_preview_is_bounded(self):
        self.assertIn("UpTo[10]", self.wl)


if __name__ == "__main__":
    unittest.main()
