"""Contract tests for read-only D(-2) C5 flavour candidate assembly."""
from __future__ import annotations

import ast
from pathlib import Path
from tempfile import TemporaryDirectory
from contextlib import redirect_stdout
from io import StringIO
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

    def test_checks_each_lepton_flavour_pair_without_matching_heads(self) -> None:
        self.assertIn("D2LeptonFlavourNames /@ allHol", self.wl)
        self.assertIn("Select[Range[Length[leptonPairs]]", self.wl)
        self.assertIn("ListQ[leptonPairs[[#]]]", self.wl)
        self.assertIn('Sort[leptonPairs[[#]]] === {"d$$1", "d$$2"}', self.wl)
        self.assertNotIn("badPairs = Flatten @ Position", self.wl)

    def test_exports_complete_lepton_index_census_before_failing(self) -> None:
        self.assertIn('profile_c5_flavour_candidate_lepton_pairs.json', self.wl)
        self.assertIn('Export[pairAuditPath, D2JSONSafe[pairRecords], "RawJSON"]', self.wl)
        self.assertLess(self.wl.index('Export[pairAuditPath'),
                        self.wl.index('badPairs = Select'))
        self.assertIn('"BadTermPositions" -> badPairs', self.wl)
        self.assertIn('"LeptonPairInventoryPath" -> pairAuditPath', self.wl)

    def test_failure_json_sanitizes_wolfram_symbols(self) -> None:
        self.assertIn("D2JSONSafe[x_List] := D2JSONSafe /@ x", self.wl)
        self.assertIn('Export[reportPath, D2JSONSafe[failure], "RawJSON"]', self.wl)

    def test_does_not_overwrite_first_threshold(self) -> None:
        self.assertNotIn('Export[required["FirstLoop"]', self.wl)
        self.assertNotIn('Export[required["HardHbar"]', self.wl)

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

    def test_running_symmetry_modulo_bound_sum_dummies(self) -> None:
        self.assertIn("HoldPattern[Sum[body_, {r_Symbol, lo_, hi_}]]", self.wl)
        self.assertIn("Inactive[Sum][body /. r -> D2AlphaDummy", self.wl)
        self.assertIn("runSymmetricAlpha = TrueQ[Expand[alphaPQ - alphaQP] === 0]", self.wl)
        self.assertIn('"RunningSymmetricModuloSummationDummies" -> runSymmetricAlpha', self.wl)
        self.assertIn("sourceSums =!= convertedSums", self.wl)

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

    def test_empty_report_is_diagnosed_not_exception(self) -> None:
        runner = self._runner()
        with TemporaryDirectory() as d:
            p = Path(d) / "result.json"
            p.write_text("", encoding="utf-8")
            report, error = runner.read_report(p)
        self.assertIsNone(report)
        self.assertIn("empty", error)

    def test_malformed_report_is_diagnosed_not_exception(self) -> None:
        runner = self._runner()
        with TemporaryDirectory() as d:
            p = Path(d) / "result.json"
            p.write_text('0.e-7', encoding="utf-8")
            report, error = runner.read_report(p)
        self.assertIsNone(report)
        self.assertIn("Invalid", error)

    def test_missing_report_is_diagnosed_not_exception(self) -> None:
        runner = self._runner()
        with TemporaryDirectory() as d:
            report, error = runner.read_report(Path(d) / "missing.json")
        self.assertIsNone(report)
        self.assertIn("did not create", error)

    def test_inspect_existing_reports_a_valid_failure_json(self) -> None:
        runner = self._runner()
        import json
        with TemporaryDirectory() as d:
            profile = Path(d) / runner.PROFILE_REL
            profile.mkdir(parents=True)
            (profile / runner.REPORT).write_text(json.dumps({
                "Status": "Failed", "Reason": "Unmatched labels",
                "ExitCode": 9, "LeptonPairInventoryPath": "pairs.json",
                "BadTermPositions": [7, 13],
                "BadTermLeptonLabels": [["d$$1"], ["d$$2"]]
            }), encoding="utf-8")
            (profile / runner.LOG).write_text('Flavour audit stopped\n', encoding='utf-8')
            stream = StringIO()
            with redirect_stdout(stream):
                code = runner.inspect_existing(Path(d))
        self.assertEqual(code, 0)
        self.assertIn('BadTermPositions: [7, 13]', stream.getvalue())
        self.assertIn('LeptonPairInventoryPath: pairs.json', stream.getvalue())

    def test_inspect_existing_does_not_run_wolfram(self) -> None:
        runner = self._runner()
        with TemporaryDirectory() as d:
            profile = Path(d) / runner.PROFILE_REL
            profile.mkdir(parents=True)
            (profile / runner.REPORT).write_text('', encoding='utf-8')
            (profile / runner.LOG).write_text('Kernel failure\n', encoding='utf-8')
            stream = StringIO()
            with redirect_stdout(stream):
                code = runner.inspect_existing(Path(d))
        self.assertEqual(code, 0)
        self.assertIn('empty', stream.getvalue())
        self.assertIn('Kernel failure', stream.getvalue())


if __name__ == "__main__":
    unittest.main()
