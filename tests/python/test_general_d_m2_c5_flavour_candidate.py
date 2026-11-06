"""Contract tests for read-only D(-2) C5 flavour candidate assembly."""
from __future__ import annotations

import ast
from pathlib import Path
import sys
import types
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "Lagrangian/RunGeneralDMinus2C5FlavourAssemblyCandidate.wl"
RUNNER = ROOT / "scripts/RunGeneralDMinus2C5FlavourAssemblyCandidate.py"


class SourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.wl = SRC.read_text(encoding="utf-8")
        cls.py = RUNNER.read_text(encoding="utf-8")

    def test_read_only_and_no_match_call(self) -> None:
        self.assertNotIn("Match[", self.wl.replace("NO Match[]", "").replace("Match,", ""))
        self.assertNotIn("RunT3Matching", self.wl)
        self.assertIn('"NoNewMatchPerformed" -> True', self.wl)

    def test_uses_saved_ef1_and_processed_hard(self) -> None:
        self.assertIn("profile_full_scalar_one_loop_hard_hbar_sector.wxf", self.wl)
        self.assertIn("eft1_after_F_realbasis_loop.wxf", self.wl)

    def test_enforces_saved_audits(self) -> None:
        self.assertIn("WeinbergStructureAudited", self.wl)
        self.assertIn("ComponentsAudited", self.wl)

    def test_checks_each_lepton_flavour_pair(self) -> None:
        self.assertIn("D2LeptonFlavourNames /@ allHol", self.wl)
        self.assertIn('Sort[pair] =!= {"d$$1", "d$$2"}', self.wl)

    def test_maps_only_flavour_indices(self) -> None:
        self.assertIn('SymbolName[fl] === "Flavor"', self.wl)
        self.assertIn('SymbolName[i] === "d$$1"', self.wl)
        self.assertIn('SymbolName[i] === "d$$2"', self.wl)

    def test_symmetric_convention_for_hard_and_inherited(self) -> None:
        self.assertIn("hardPQ + hardQP", self.wl)
        self.assertIn("inhPQ + inhQP", self.wl)

    def test_running_not_resymmetrized(self) -> None:
        self.assertIn("candidateRunning = hbar runPQ", self.wl)
        self.assertNotIn("candidateRunning = hbar (runPQ + runQP)", self.wl)

    def test_hbar_gate(self) -> None:
        self.assertIn("Coefficient[candidateTotal, hbar, 1]", self.wl)

    def test_all_provisional_flags(self) -> None:
        self.assertIn('"FixedOrderSMEFTCertified" -> False', self.wl)
        self.assertIn('"OperatorNormalizationIndependentlyCertified" -> False', self.wl)
        self.assertIn('"InheritedSMCanonicalizationCertified" -> False', self.wl)

    def test_runner_syntax(self) -> None:
        ast.parse(self.py)


class RunnerTests(unittest.TestCase):
    def _runner(self):
        fake = types.ModuleType("scripts.RunGeneralDMinus2ScalarThreshold")
        fake.preflight = lambda *a, **kw: {"insertion_file": "any"}
        with patch.dict(sys.modules, {fake.__name__: fake}):
            import importlib.util
            spec = importlib.util.spec_from_file_location("candidate_runner_test", RUNNER)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
        return module

    def _passed(self):
        return {
            "Status": "C5FlavourCandidateExported", "HardHolomorphicTermCount": 34,
            "InheritedHolomorphicTermCount": 6, "NoNewMatchPerformed": True,
            "EachHolomorphicTermHasTwoDifferentLeptonFlavourLabels": True,
            "HardSymmetricByConstruction": True, "InheritedSymmetricByConstruction": True,
            "RunningLeftUnsymmetrized": True, "HardCoefficientHasNativeHbar": True,
            "InheritedCoefficientHasNativeHbar": True, "RunningInsertedWithHbar": True,
            "TotalCandidateLinearInHbar": True, "FixedOrderC5CandidateAssembled": True,
            "OperatorNormalizationIndependentlyCertified": False,
            "InheritedSMCanonicalizationCertified": False,
            "PhysicalScalarMassRelationsImposed": False,
            "FixedOrderSMEFTCertified": False
        }

    def test_successful_classification(self) -> None:
        self.assertEqual(self._runner().classify(0, False, self._passed()),
                         "C5FlavourCandidateExported")

    def test_cannot_claim_certified(self) -> None:
        report = self._passed()
        report["FixedOrderSMEFTCertified"] = True
        self.assertEqual(self._runner().classify(0, False, report),
                         "InconsistentCandidateMetadata")

    def test_running_hbar_guard(self) -> None:
        report = self._passed()
        report["RunningInsertedWithHbar"] = False
        self.assertEqual(self._runner().classify(0, False, report),
                         "InconsistentCandidateMetadata")

    def test_nonzero_exit(self) -> None:
        self.assertEqual(self._runner().classify(9, False, self._passed()),
                         "WolframFailed")

    def test_timeout(self) -> None:
        self.assertEqual(self._runner().classify(None, True, self._passed()),
                         "ExternalTimeout")


if __name__ == "__main__":
    unittest.main()
