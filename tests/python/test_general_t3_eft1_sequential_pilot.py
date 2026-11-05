"""Fail-closed general B(-3) EFT1 boundary and tensor-metadata tests."""

from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.RunGeneralT3EFT1SequentialPilot import (  # noqa: E402
    GENERAL_CG,
    GENERAL_COUPLING,
    GENERAL_INTERACTION,
    GENERAL_RGBETA_COUPLING,
    make_tensor_metadata,
    validate_general_boundary,
)


class GeneralBMinus3BoundaryTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.base = Path(self.tmp.name)
        self.matching = self.base / "matching"
        self.general = self.base / "general"
        self.good_fixture()

    def write(self, base: Path, name: str, obj: dict) -> None:
        path = base / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(obj, indent=2), encoding="utf-8")

    def read(self, base: Path, name: str) -> dict:
        return json.loads((base / name).read_text(encoding="utf-8"))

    def good_fixture(self) -> None:
        interactions = ["Yukawa1_GenericCG", "Yukawa2_GenericCG",
                        "T3ScalarMix_GenericCG", GENERAL_INTERACTION]
        self.write(self.general, "lagrangian.json", {
            "model": "B_-3", "status": "Success", "interactions": interactions,
            "quartic_expected": 1, "quartic_included": 1,
            "quadratic_expected": 0, "quadratic_included": 0,
            "extra_yukawa_expected": 0, "extra_yukawa_included": 0,
            "general_scalar_quadratics_requested": True,
            "extra_yukawas_requested": True,
        })
        self.write(self.general, "status.json", {
            "status": "SymbolicRGExportedMatchingUnverified", "errors": [],
        })
        self.write(self.general, "eft1_rgbeta.json", {
            "status": "SuccessCandidateMassBeta",
            "metadata": {"GeneralScalarModel": "B_-3", "GeneralScalarStage": "EFT1"},
            "report_betas": {GENERAL_RGBETA_COUPLING: "some beta expression"},
        })
        self.write(self.matching, "comparison_summary.json", {
            "BuildStatus": "Success", "MatchingStatus": "Success",
            "WeinbergExtractionStatus": "Success",
            "AcceptedInteractions": interactions,
            "ThresholdPlan": [["F"], ["S1", "S2"]],
        })
        self.write(self.matching, "sequential_matching_debug.json", {
            "Status": "Success", "RequestedThresholdCount": 2,
            "MatchedStageCount": 2,
        })
        self.write(self.matching, "data/eft1_after_F_rge_export.json", {
            "IntegratedFields": ["F"], "ActiveHeavyFields": ["S1", "S2"],
            "TreeDiagnostics": {"Fields": ["H", "NewScalar1", "NewScalar2", "l"]},
        })
        self.write(self.matching, "data/eft1_after_F_wilson_seed.json", {
            "TreeWilsonTermCount": 1,
            "TreeWilsonTerms": [{"Chirality": "PL", "TermInputForm": "..."}],
        })
        self.write(self.matching, "data/eft1_after_F_scalar_quartic_seed.json", {
            "ScalarQuarticTermCount": 2,
            "ScalarQuarticTerms": [
                {"CGNames": [GENERAL_CG],
                 "TermInputForm": f"Coupling[{GENERAL_COUPLING}, {{}}, 0]*CG[{GENERAL_CG}, {{{{}}}}]"},
                {"CGNames": [GENERAL_CG],
                 "TermInputForm": f"Bar[Coupling[{GENERAL_COUPLING}, {{}}, 0]]*CG[Bar[{GENERAL_CG}], {{{{}}}}]"},
            ],
            "CGRegistry": [{"Name": GENERAL_CG, "TensorInputForm": "SparseArray[...]",
                            "RepsInputForm": "{...}"}],
        })

    def test_good_general_boundary(self) -> None:
        report = validate_general_boundary(self.matching, self.general)
        self.assertTrue(report["general_quartic_cg_registry_present"])
        self.assertEqual(report["general_quartic_and_conjugate_in_eft1"], 2)
        self.assertEqual(report["tree_pl_wilson_term_count"], 1)
        self.assertFalse(report["physical_matchete_rgbeta_mapping_certified"])

    def test_restricted_uv_fails(self) -> None:
        file = "comparison_summary.json"
        data = self.read(self.matching, file)
        data["AcceptedInteractions"].remove(GENERAL_INTERACTION)
        self.write(self.matching, file, data)
        with self.assertRaisesRegex(ValueError, "different UV interaction inventory"):
            validate_general_boundary(self.matching, self.general)

    def test_wrong_threshold_order_fails(self) -> None:
        file = "comparison_summary.json"
        data = self.read(self.matching, file)
        data["ThresholdPlan"] = [["S1", "S2"], ["F"]]
        self.write(self.matching, file, data)
        with self.assertRaisesRegex(ValueError, "sequential threshold plan"):
            validate_general_boundary(self.matching, self.general)

    def test_missing_general_quartic_cg_fails(self) -> None:
        file = "data/eft1_after_F_scalar_quartic_seed.json"
        data = self.read(self.matching, file)
        data["CGRegistry"] = []
        self.write(self.matching, file, data)
        with self.assertRaisesRegex(ValueError, "CG tensor was not exported"):
            validate_general_boundary(self.matching, self.general)

    def test_fermion_must_be_removed_at_first_threshold(self) -> None:
        file = "data/eft1_after_F_rge_export.json"
        data = self.read(self.matching, file)
        data["TreeDiagnostics"]["Fields"].append("NewFermion")
        self.write(self.matching, file, data)
        with self.assertRaisesRegex(ValueError, "remove the fermion"):
            validate_general_boundary(self.matching, self.general)

    def test_general_rgbeta_quartic_must_exist(self) -> None:
        file = "eft1_rgbeta.json"
        data = self.read(self.general, file)
        data["report_betas"] = {}
        self.write(self.general, file, data)
        with self.assertRaisesRegex(ValueError, "RGBeta omitted"):
            validate_general_boundary(self.matching, self.general)

    def test_tensor_metadata_is_not_restricted_rgbeta(self) -> None:
        value = make_tensor_metadata()
        meta = value["metadata"]
        self.assertEqual(value["status"], "MetadataOnlyGeneralEFT1")
        self.assertEqual((meta["dS1"], meta["dS2"], meta["dF"]), (2, 2, 1))
        self.assertEqual((meta["YS1"], meta["YS2"]), ("-3/2", "-1/2"))
        self.assertEqual(meta["ActiveBSMFields"], ["S1", "S2"])
        self.assertNotIn("report_betas", value)


if __name__ == "__main__":
    unittest.main()
