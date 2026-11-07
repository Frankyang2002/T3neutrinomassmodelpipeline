"""Contract tests for the isolated Warsaw normalization-reference diagnostic."""
from __future__ import annotations

import importlib.util
from pathlib import Path
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "Lagrangian" / "RunGeneralDMinus2C5NormalizationReference.wl"
PYTHON = ROOT / "scripts" / "RunGeneralDMinus2C5NormalizationReference.py"


class NormalizationReferenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.wl = SCRIPT.read_text(encoding="utf-8")
        spec = importlib.util.spec_from_file_location("d2_norm_reference", PYTHON)
        assert spec is not None and spec.loader is not None
        cls.mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.mod)

    def test_source_exists(self) -> None:
        self.assertTrue(SCRIPT.is_file())
        self.assertTrue(PYTHON.is_file())

    def test_has_no_match_invocation(self) -> None:
        import re
        uncommented = re.sub(r"\(\*.*?\*\)", "", self.wl, flags=re.S)
        self.assertNotRegex(uncommented, r"\bMatch\s*\[")

    def test_loads_warsaw_reference_not_d2_uv(self) -> None:
        self.assertIn('LoadModel["SMEFT_Warsaw"]', self.wl)
        self.assertNotIn('BuildT3Lagrangian[', self.wl)

    def test_extracts_current_project_weinberg(self) -> None:
        self.assertIn('Get[FileNameJoin[{scriptDir, "RunMatching.wl"}]]', self.wl)
        self.assertIn('ExtractWeinbergCoefficient[refExpr]', self.wl)

    def test_strictly_checks_warsaw_llhh(self) -> None:
        self.assertIn('D2DegreeNorm[#, "H"] === 2', self.wl)
        self.assertIn('D2DegreeNorm[#, "l"] === 2', self.wl)
        self.assertIn('!FreeQ[#, Matchete`GammaCC]', self.wl)

    def test_respects_wilson_coupling(self) -> None:
        self.assertIn('D2UnitWilsonNorm[refOrdered]', self.wl)
        self.assertIn('WarsawReferenceCouplingNames', self.wl)

    def test_checks_extractor_linearity(self) -> None:
        self.assertIn('ExtractWeinbergCoefficient[3 refExpr]', self.wl)
        self.assertIn('scaledCoefficient - 3 refOrdered', self.wl)

    def test_has_unitary_gauge_vev_projection(self) -> None:
        self.assertIn('vD2Reference/Sqrt[2]', self.wl)
        self.assertIn('ewsbNN - vD2Reference^2/2', self.wl)

    def test_does_not_certify_unvalidated_physical_norm(self) -> None:
        self.assertIn('"FullPhysicalC5NormalizationCertified" -> False', self.wl)
        self.assertIn('"FullDMinus2SU2CGNormalizationCertified" -> False', self.wl)

    def test_json_preflight(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            d = Path(temp)
            with self.assertRaises(FileNotFoundError):
                self.mod.read_object(d / 'missing.json')

    def test_report_classification_valid(self) -> None:
        report = {
            "Status": "WarsawNormalizationReferenceAudited",
            "NoNewMatching": True,
            "ProvisionalCandidateUnchanged": True,
            "ProjectExtractorRecognizesWarsawWeinbergOperator": True,
            "ProjectExtractorLinearOnWarsawReference": True,
            "WarsawOperatorNeutrinoVEVComponentIsPlusVSquaredOverTwo": True,
            "FullPhysicalC5NormalizationCertified": False,
            "PhysicalMassRelationsImposed": False,
            "NumericalBenchmarkPerformed": False,
        }
        self.assertEqual(self.mod.classify(0, False, report),
                         'WarsawNormalizationReferenceAudited')
        report['ProjectExtractorLinearOnWarsawReference'] = False
        self.assertEqual(self.mod.classify(0, False, report),
                         'InconsistentReferenceReport')

    def test_report_classification_failure(self) -> None:
        self.assertEqual(self.mod.classify(12, False, None), 'WolframFailed')
        self.assertEqual(self.mod.classify(None, True, None), 'ExternalTimeout')
        self.assertEqual(self.mod.classify(0, False, None), 'MissingOrInvalidResult')

    def test_runner_is_read_only(self) -> None:
        self.assertIn('CANDIDATE_WXF =', PYTHON.read_text())
        self.assertIn('No matching, WXF rewriting', PYTHON.read_text())

    def test_rejects_unverified_candidate(self) -> None:
        import json
        with tempfile.TemporaryDirectory() as temp:
            d = Path(temp)
            profile = d / self.mod.PROFILE_REL
            profile.mkdir(parents=True)
            (profile / self.mod.CANDIDATE_WXF).write_bytes(b'test')
            (profile / self.mod.CANDIDATE_REPORT).write_text(json.dumps({
                "Status": "C5FlavourCandidateExported",
                "FixedOrderC5CandidateAssembled": True,
                "RunningSymmetricModuloSummationDummies": False,
                "TotalCandidateLinearInHbar": True,
            }))
            with self.assertRaises(ValueError):
                self.mod.preflight(d)


if __name__ == "__main__":
    unittest.main()
