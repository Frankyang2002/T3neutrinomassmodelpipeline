"""Contract checks for the additive read-only Weinberg field-content audit."""
import importlib.util
from pathlib import Path
import sys
import types
import unittest

ROOT = Path(__file__).resolve().parents[2]
RUNNER = ROOT / "scripts/RunGeneralDMinus2C5FlavourStructureAudit.py"
WL = ROOT / "Lagrangian/RunGeneralDMinus2C5FlavourStructureAudit.wl"


def load_runner():
    name = "scripts.RunGeneralDMinus2ScalarThreshold"
    stub = types.ModuleType(name)
    stub.preflight = lambda *a, **k: {"insertion_file": "/dummy/running.wl"}
    previous = sys.modules.get(name)
    sys.modules[name] = stub
    try:
        spec = importlib.util.spec_from_file_location("D2C5StructureAudit", RUNNER)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    finally:
        if previous is None:
            sys.modules.pop(name, None)
        else:
            sys.modules[name] = previous


class D2C5StructureAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.runner = load_runner()
        cls.source = WL.read_text(encoding="utf-8")

    @staticmethod
    def good():
        return {
            "Status": "WeinbergStructureAudited",
            "NoNewMatching": True,
            "AllCandidateTermsHaveTwoLeptonsTwoHiggs": True,
            "AllCandidateTermsHaveUnambiguousChirality": True,
            "HardC5IsLinearInHbar": True,
            "InheritedC5IsLinearInHbar": True,
            "FlavourIndexMappingCertified": False,
            "C5SymmetricNormalizationCertified": False,
            "FixedOrderSMEFTAssembled": False,
        }

    def test_success(self):
        self.assertEqual(self.runner.classify(0, False, self.good()),
                         "WeinbergStructureAudited")

    def test_nonzero_exit(self):
        self.assertEqual(self.runner.classify(20, False, self.good()), "WolframFailed")

    def test_timeout(self):
        self.assertEqual(self.runner.classify(0, True, self.good()), "ExternalTimeout")

    def test_missing_result(self):
        self.assertEqual(self.runner.classify(0, False, None), "MissingResult")

    def test_incorrect_field_content(self):
        r = self.good()
        r["AllCandidateTermsHaveTwoLeptonsTwoHiggs"] = False
        self.assertEqual(self.runner.classify(0, False, r), "InvalidStructureChecks")

    def test_extra_loop_factor(self):
        r = self.good()
        r["InheritedC5IsLinearInHbar"] = False
        self.assertEqual(self.runner.classify(0, False, r), "InvalidStructureChecks")

    def test_no_premature_certification(self):
        r = self.good()
        r["FixedOrderSMEFTAssembled"] = True
        self.assertEqual(self.runner.classify(0, False, r), "PrematureC5Certification")

    def test_only_saved_results(self):
        self.assertIn('profile_full_scalar_one_loop_hard_hbar_sector.wxf', self.source)
        self.assertIn('eft1_after_F_realbasis_loop.wxf', self.source)
        self.assertIn('profile_fixed_order_components_result.json', self.source)
        self.assertNotIn('RunT3Matching[', self.source)
        self.assertNotIn('Match[matchInput', self.source)

    def test_degree_handling(self):
        self.assertIn('D2FieldDegree[Power[x_, n_Integer?Positive], target_String]', self.source)
        self.assertIn('selfTest = {D2FieldDegree[hSample, "H"]', self.source)

    def test_original_extractor_and_full_records(self):
        self.assertIn('Get[FileNameJoin[{scriptDir, "RunMatching.wl"}]]', self.source)
        self.assertIn('ExtractWeinbergCoefficient[hard]', self.source)
        self.assertIn('ExtractWeinbergCoefficient[inheritedHbar]', self.source)
        self.assertIn('"FlavourIndexLabels" -> D2FlavourNames[t]', self.source)

    def test_comparison_to_previous_counts(self):
        self.assertIn('"InheritedPureSMWeinbergTermCount"', self.source)
        self.assertIn('"HardWeinbergTermCount"', self.source)

    def test_no_unjustified_flavour_claim(self):
        self.assertIn('"FlavourIndexMappingCertified" -> False', self.source)
        self.assertIn('"C5SymmetricNormalizationCertified" -> False', self.source)


if __name__ == "__main__":
    unittest.main()
