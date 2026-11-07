"""Static and runner-contract tests for the standalone C5 pole audit."""
from __future__ import annotations
import importlib.util
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
WL = ROOT / "Lagrangian/RunGeneralDMinus2C5PoleCancellationAudit.wl"
PY = ROOT / "scripts/RunGeneralDMinus2C5PoleCancellationAudit.py"


class PoleAuditContract(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.wl = WL.read_text(encoding="utf-8")
        spec = importlib.util.spec_from_file_location("d2_pole_module", PY)
        assert spec and spec.loader
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        cls.mod = mod

    def test_saved_candidate_only(self):
        self.assertIn("profile_c5_flavour_candidate_total_provisional.wxf", self.wl)
        self.assertNotIn("Match[", self.wl.replace("Match[]", ""))

    def test_pole_is_not_set_numerically(self):
        self.assertIn("D2PoleResidue[expression_]", self.wl)
        self.assertIn("D2PoleMarker", self.wl)
        self.assertNotIn("1/10^-", self.wl)

    def test_exact_nf1(self):
        self.assertIn("If[nf =!= 1", self.wl)
        self.assertIn("D2ReduceNF1", self.wl)
        self.assertIn('D2Name[rep] === "NFlavor"', self.wl)

    def test_symbolic_reduction(self):
        self.assertIn("Cancel[Together[resNF1]]", self.wl)
        self.assertIn("ExactPoleCancellationVerifiedForNF1", self.wl)

    def test_renormalization_not_claimed(self):
        self.assertIn('"RenormalizationCertified" -> False', self.wl)

    def test_result_classification(self):
        data = {"Status": "ExactNF1PoleCancellationVerified", "NF": 1,
                "NoNewMatchPerformed": True, "RenormalizationCertified": False}
        self.assertEqual(self.mod.classify(0, False, data), data["Status"])
        self.assertEqual(self.mod.classify(1, False, data), "WolframFailed")
        self.assertEqual(self.mod.classify(0, True, data), "ExternalTimeout")
        self.assertEqual(self.mod.classify(0, False, None), "InvalidResultReport")
        self.assertEqual(self.mod.classify(0, False, {**data, "NF": 2}), "InconsistentMetadata")

    def test_read_empty_json(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "result.json"
            path.write_text("", encoding="utf-8")
            value, error = self.mod.read_result(path)
            self.assertIsNone(value)
            self.assertIn("Empty", error)

    def test_read_valid_json(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "result.json"
            path.write_text('{"Status":"OK"}', encoding="utf-8")
            value, error = self.mod.read_result(path)
            self.assertEqual(value["Status"], "OK")
            self.assertIsNone(error)


if __name__ == "__main__":
    unittest.main()
