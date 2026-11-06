"""Contract tests; do not claim live Matchete execution from these tests."""
from __future__ import annotations

import unittest
from pathlib import Path

from scripts.RunGeneralDMinus2ExactBSMProbe import (
    ROOT, SCRIPT, classify_result, make_command,
)


class TestD2ExactBSMProbe(unittest.TestCase):
    def test_command_uses_separate_wolfram_file_and_exact_matching_dir(self):
        command = make_command(Path("/tmp/stage1"), Path("/tmp/insertion.wl"),
                               timeout_seconds=30, script=SCRIPT)
        self.assertEqual(command[1], "-file")
        self.assertEqual(command[-1], "30")
        self.assertTrue(command[2].endswith("RunGeneralDMinus2ExactBSMProbe.wl"))

    def test_illegal_timeout_rejected(self):
        with self.assertRaises(ValueError):
            make_command(Path("/tmp/a"), Path("/tmp/b"), timeout_seconds=0)

    @staticmethod
    def _audit():
        return {
            "Status": "Passed",
            "SMAlphaComparisonCertified": True,
            "SavedSMAlphaEquivalentToOriginalLSM": True,
            "SavedSMAlphaEquivalentToRelabelledLSM": False,
            "StandaloneSMTermsRemaining": 0,
            "AllHeavyTermsPreservedExactly": True,
            "ReconstructionResidualExactlyZero": True,
            "SMSubtractandContainsHeavyFields": False,
            "SMSubtractandExactlyEqualsSavedPureSM": True,
        }

    @staticmethod
    def _result():
        return {
            "Status": "CompleteTreeDiagnostic",
            "OneLoopScalarHardMatchingPerformed": False,
            "RemainingHeavyFields": [],
            "Stages": [{"Stage": stage, "Status": "Returned"} for stage in (
                "Match", "GreensSimplify", "EOMSimplify",
                "EvaluateLoopFunctions", "ReplaceEffectiveCouplings")],
        }

    def test_complete_tree_diagnostic_not_physical_one_loop(self):
        self.assertEqual(classify_result(exit_code=0, external_timeout=False,
                                         audit=self._audit(),
                                         result=self._result()),
                         "CompleteTreeDiagnostic")

    def test_return_from_match_is_not_enough(self):
        result = self._result()
        result["Stages"] = result["Stages"][:1]
        self.assertEqual(classify_result(exit_code=0, external_timeout=False,
                                         audit=self._audit(), result=result),
                         "MissingTreeProcessingStage")

    def test_subtraction_audit_is_mandatory(self):
        self.assertEqual(classify_result(exit_code=0, external_timeout=False,
                                         audit={"Status": "Failed"},
                                         result=self._result()),
                         "AuditMissingOrFailed")

    def test_cannot_accept_uncertified_dummy_index_comparison(self):
        audit = self._audit()
        audit["SMAlphaComparisonCertified"] = False
        self.assertEqual(classify_result(exit_code=0, external_timeout=False,
                                         audit=audit, result=self._result()),
                         "SMAlphaComparisonNotCertified")

    def test_cannot_accept_remaining_sm_terms(self):
        audit = self._audit()
        audit["StandaloneSMTermsRemaining"] = 1
        self.assertEqual(classify_result(exit_code=0, external_timeout=False,
                                         audit=audit, result=self._result()),
                         "SMSubtractionAuditIncomplete")

    def test_cannot_accept_changed_heavy_sector(self):
        audit = self._audit()
        audit["AllHeavyTermsPreservedExactly"] = False
        self.assertEqual(classify_result(exit_code=0, external_timeout=False,
                                         audit=audit, result=self._result()),
                         "SMSubtractionAuditIncomplete")

    def test_external_timeout_takes_priority(self):
        self.assertEqual(classify_result(exit_code=0, external_timeout=True,
                                         audit=self._audit(),
                                         result=self._result()),
                         "ExternalTimeout")

    def test_wl_contract_contains_exact_eft1_not_incremental_profile(self):
        wl = (ROOT / "Lagrangian/RunGeneralDMinus2ExactBSMProbe.wl")
        source = wl.read_text(encoding="utf-8")
        for fragment in (
            '"eft1_after_F_realbasis_tree.wxf"',
            'smRelabelled = Expand[RelabelIndices[uvBuild["LSM"]]]',
            "D2AlphaIndexKey[term_]",
            "sourceSM = Select[allTerms, !D2HasHeavy[#] &]",
            "smAlphaCertified = smAlphaOriginal || smAlphaRelabelled",
            "smStandalone = Total[sourceSM]",
            "exactBSM = Expand[exactTree - smStandalone]",
            "reconstructionResidual = Expand[exactBSM + smStandalone - exactTree]",
            "heavyPreserved", "remainingSM", "wilsonOverlap", "missingExtra",
            "inputSafe = transitionEqual && smAlphaCertified",
            'D2Step["Match"', 'D2Step["GreensSimplify"',
            'D2Step["EOMSimplify"', 'D2Step["EvaluateLoopFunctions"',
            'D2Step["ReplaceEffectiveCouplings"',
        ):
            self.assertIn(fragment, source)
        self.assertNotIn('profileBSMWithBoth =', source)


if __name__ == "__main__":
    unittest.main()
