"""Contract tests for read-only D(-2) fixed-order components audit."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import types
import unittest

ROOT = Path(__file__).resolve().parents[2]
RUNNER = ROOT / "scripts/RunGeneralDMinus2FixedOrderComponents.py"
WL = ROOT / "Lagrangian/RunGeneralDMinus2FixedOrderComponents.wl"


def load_runner():
    name = "scripts.RunGeneralDMinus2ScalarThreshold"
    stub = types.ModuleType(name)
    stub.preflight = lambda *a, **k: {"insertion_file": "/dummy/running.wl"}
    old = sys.modules.get(name)
    sys.modules[name] = stub
    try:
        spec = importlib.util.spec_from_file_location("D2FixedOrderProbe", RUNNER)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    finally:
        if old is None:
            sys.modules.pop(name, None)
        else:
            sys.modules[name] = old


class FixedOrderComponentsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.runner = load_runner()
        cls.source = WL.read_text(encoding="utf-8")

    @staticmethod
    def good_report():
        return {
            "Status": "ComponentsAudited",
            "NoNewMatchPerformed": True,
            "TreeHasLinearHeavySource": False,
            "SavedHardSectorMatchesProcessedO1": True,
            "SavedHardSectorContainsHeavy": False,
            "HardC5LinearInHbar": True,
            "HardCandidateMatchesSavedText": True,
            "DirectRunningHasLog": True,
            "DirectRunningHasExplicitHbar": False,
            "InheritedPureSMHasZeroHbarPart": True,
            "InheritedPureSMHasHbarSquared": False,
            "HardWeinbergExtractionStatus": "Success",
            "InheritedPureSMWeinbergExtractionStatus": "NotFound",
            "FlavourAssemblyCertified": False,
            "FixedOrderSMEFTAssembled": False,
            "FullGeneralScalarThresholdCertified": False,
        }

    def test_good(self):
        self.assertEqual(self.runner.classify(0, False, self.good_report()),
                         "ComponentsAudited")

    def test_nonzero_exit(self):
        self.assertEqual(self.runner.classify(20, False, self.good_report()),
                         "WolframFailed")

    def test_timeout(self):
        self.assertEqual(self.runner.classify(0, True, self.good_report()),
                         "ExternalTimeout")

    def test_missing_report(self):
        self.assertEqual(self.runner.classify(0, False, None),
                         "InvalidComponentsReport")

    def test_pending_inherited(self):
        r = self.good_report()
        r["InheritedPureSMWeinbergExtractionStatus"] = "CoefficientPending"
        self.assertEqual(self.runner.classify(0, False, r),
                         "InheritedC5ExtractionPending")

    def test_no_false_assembly_certification(self):
        r = self.good_report()
        r["FixedOrderSMEFTAssembled"] = True
        self.assertEqual(self.runner.classify(0, False, r),
                         "InconsistentComponentsReport")

    def test_no_hbar_double_counting(self):
        r = self.good_report()
        r["DirectRunningHasExplicitHbar"] = True
        self.assertEqual(self.runner.classify(0, False, r),
                         "InconsistentComponentsReport")

    def test_only_saved_inputs(self):
        for name in ("eft1_after_F_realbasis_tree.wxf",
                     "eft1_after_F_realbasis_loop.wxf",
                     "profile_full_scalar_one_loop_processed.wxf",
                     "profile_full_scalar_one_loop_hard_hbar_sector.wxf"):
            self.assertIn(name, self.source)
        self.assertNotIn("Match[matchInput", self.source)
        self.assertNotIn("RunT3Matching[", self.source)
        self.assertIn('"FixedOrderSMEFTAssembled" -> False', self.source)

    def test_projected_inherited_separate(self):
        self.assertIn('inheritedPureSMTerms = Select[loopTerms, !D2ContainsHeavy[#] &]',
                      self.source)
        self.assertIn('ExtractWeinbergCoefficient[inheritedHbar]', self.source)
        self.assertIn('TreeHasLinearHeavySource', self.source)

    def test_project_original_extractor(self):
        self.assertIn('Get[FileNameJoin[{scriptDir, "RunMatching.wl"}]]', self.source)
        self.assertIn('ExtractWeinbergCoefficient[hardHbarSaved]', self.source)

    def test_running_external_hbar(self):
        self.assertIn('Get[insertionPath]', self.source)
        self.assertIn('EFT1D2WilsonRunningWeinbergCoefficient[p, q]', self.source)
        self.assertIn('directHasHbar = !FreeQ[directC5, hbar]', self.source)

    def test_json_read(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "x.json"
            p.write_text(json.dumps(self.good_report()), encoding="utf-8")
            self.assertEqual(self.runner.read_json(p), self.good_report())


if __name__ == "__main__":
    unittest.main()
