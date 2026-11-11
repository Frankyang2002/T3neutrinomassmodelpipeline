"""Offline Worker C source-gate regressions; no Wolfram or Matchete needed."""
from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

from Numerical.joint_scalar_matching.worker_c_matchete_audit import (
    FILES, audit_matchete_scalar_profile, save_audit,
)


def make_profile(root: Path, *, y2conj: bool = True) -> None:
    documents = {
        "match_result": {
            "Status": "FullEFT1OneLoopPostprocessed",
            "OneLoopRawMatchReturned": True, "LoopOrder": 1, "EFTOrder": 5,
            "SimplificationPerformed": True, "HardC5Extracted": True,
            "WeinbergExtractionStatus": "Success", "WeinbergTermCount": 68,
            "FullFlavorHardC5Verified": False, "FullGeneralScalarThresholdCertified": False,
            "InheritedOneLoopAdded": False, "DirectRunningAdded": False,
            "PhysicalMassRelationsImposed": False,
        },
        "mass_audit": {"Status": "Passed", "SMAlphaComparisonCertified": True,
            "AllHeavyTermsPreservedExactly": True,
            "AllPhysicalScalarsRegisteredHeavy": True,
            "PureHeavyQuadraticAlphaEqualToFreeLag": True,
            "CheckLagrangianPassed": True,
            "Y1Y2WilsonTermCount": 4,
            "Y1Y2ConjWilsonTermCount": 4 if y2conj else 0},
        "fixed_order": {"HardWeinbergHolomorphicTermCount": 34,
            "HardWeinbergConjugateTermCount": 34,
            "HardC5HasY2Conj": y2conj, "FlavourAssemblyCertified": False},
        "flavour_candidate": {"OperatorNormalizationIndependentlyCertified": False},
        "pole_audit": {"ExactPoleCancellationVerifiedForNF1": True,
            "RenormalizationCertified": False, "NF": 1},
        "normalization_audit": {"FullPhysicalC5NormalizationCertified": False},
    }
    for key, doc in documents.items():
        (root / FILES[key]).write_text(json.dumps(doc), encoding="utf-8")
    symbol = "y2Conj" if y2conj else "y2"
    (root / FILES["candidate_inputform"]).write_text(
        f"hbar*Coupling[y1]*Coupling[{symbol}]*(1 + 1/ϵ)", encoding="utf-8")


class TestWorkerCMatcheteAudit(unittest.TestCase):
    def test_existing_d_minus2_like_source_blocks(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            make_profile(root)
            report = audit_matchete_scalar_profile(root)
            self.assertEqual(report["status"], "BLOCKED")
            self.assertEqual(report["candidate_inventory"]["total_weinberg_terms"], 68)
            self.assertEqual(report["candidate_inventory"]["holomorphic_hard_terms"], 34)
            self.assertTrue(report["gates"]["pole_cancellation_audited_at_reported_nf"])
            self.assertTrue(report["candidate_inventory"]["raw_hard_has_y2conj"])
            self.assertTrue(report["candidate_inventory"]["raw_hard_has_explicit_epsilon_symbol"])
            self.assertFalse(report["safe_to_add_as_independent_c5_correction"])

    def test_setting_no_conjugate_is_not_sufficient(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            make_profile(root, y2conj=False)
            report = audit_matchete_scalar_profile(root)
            self.assertTrue(report["gates"]["no_independent_conjugate_scalar_yukawa_in_source"])
            self.assertFalse(report["safe_to_add_as_independent_c5_correction"])
            self.assertFalse(report["gates"]["independent_nonbubble_remainder_certified"])

    def test_missing_profile_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            report = audit_matchete_scalar_profile(tmp)
            self.assertEqual(report["status"], "BLOCKED")
            self.assertFalse(report["gates"]["source_files_complete"])
            self.assertEqual(len(report["blocking_reasons"]), len(FILES) + len(
                [k for k, v in report["gates"].items() if not v and k in (
                    "source_files_complete", "exact_eft1_partition_passed",
                    "all_scalars_heavy_same_threshold", "raw_one_loop_match_returned",
                    "matchete_postprocessing_returned", "one_loop_weinberg_terms_present",
                    "no_independent_conjugate_scalar_yukawa_in_source",
                    "hard_candidate_has_no_unsubtracted_poles", "full_flavour_hard_c5_verified",
                    "operator_normalization_certified", "finite_renormalization_certified",
                    "physical_mass_relations_imposed", "independent_nonbubble_remainder_certified",
                    "physical_full_scalar_threshold_certified")]))

    def test_non_json_is_blocked(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            make_profile(root)
            (root / FILES["pole_audit"]).write_text("bad json", encoding="utf-8")
            report = audit_matchete_scalar_profile(root)
            self.assertFalse(report["gates"]["source_files_complete"])
            self.assertEqual(report["status"], "BLOCKED")

    def test_reports_have_sha_and_never_overwrite(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            make_profile(root)
            report = audit_matchete_scalar_profile(root)
            self.assertRegex(report["source_files"]["match_result"]["sha256"], r"^[0-9a-f]{64}$")
            saved = save_audit(report, root / "new_output")
            self.assertTrue(saved.is_file())
            self.assertFalse(json.loads(saved.read_text(encoding="utf-8"))[
                "safe_to_add_as_independent_c5_correction"])
            with self.assertRaises(FileExistsError):
                save_audit(report, root / "new_output")


if __name__ == "__main__":
    unittest.main()
