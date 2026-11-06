"""Static and Python contract tests; real Wolfram execution is not available here."""
from __future__ import annotations

import unittest
from pathlib import Path

from scripts.RunGeneralDMinus2CanonicalCoreProbe import (
    ROOT, SCRIPT, classify_result, make_command,
)


class TestCanonicalCoreProbe(unittest.TestCase):
    @staticmethod
    def audit() -> dict:
        return {
            "Status": "Passed",
            "SMAlphaComparisonCertified": True,
            "SavedSMAlphaEquivalentToOriginalLSM": True,
            "SavedSMAlphaEquivalentToRelabelledLSM": True,
            "StandaloneSMTermsRemaining": 0,
            "AllHeavyTermsPreservedExactly": True,
            "ReconstructionResidualExactlyZero": True,
            "SMSubtractandContainsHeavyFields": False,
            "SMSubtractandExactlyEqualsSavedPureSM": True,
            "InputMode": "ExactHeavyPlusSavedSMKineticCore",
            "CorePlusDroppedSMReconstructsFullEFT1": True,
            "CoreHeavySectorExactlyPreserved": True,
            "UnclassifiedCouplinglessSMTerms": 0,
            "CheckLagrangianPassed": True,
            "FullPhysicalSMIncluded": False,
            "SMKineticTermCount": 6,
            "SMFermionKineticTermCount": 5,
            "SMHiggsKineticTermCount": 1,
            "SMGaugeKineticTermCount": 0,
            "ExplicitGaugeKineticPresenceRequired": False,
            "DroppedSMInteractionTermCount": 11,
        }

    @staticmethod
    def result() -> dict:
        return {
            "Status": "CompleteCanonicalCoreTreeDiagnostic",
            "InputMode": "ExactHeavyPlusSavedSMKineticCore",
            "MatchInputCheckLagrangianPassed": True,
            "SMKineticCoreRestored": True,
            "FullPhysicalSMIncluded": False,
            "OneLoopScalarHardMatchingPerformed": False,
            "RemainingHeavyFields": [],
            "Stages": [
                {"Stage": x, "Status": "Returned", "Head": "Plus"}
                for x in ("Match", "GreensSimplify", "EOMSimplify",
                          "EvaluateLoopFunctions", "ReplaceEffectiveCouplings")
            ],
        }

    def classify(self, audit=None, result=None, exit_code=0, external_timeout=False):
        return classify_result(exit_code=exit_code, external_timeout=external_timeout,
                               audit=self.audit() if audit is None else audit,
                               result=self.result() if result is None else result)

    def test_command_path_and_timeout(self):
        argv = make_command(Path('/tmp/stage'), Path('/tmp/insert.wl'),
                            timeout_seconds=90)
        self.assertTrue(argv[2].endswith('RunGeneralDMinus2CanonicalCoreProbe.wl'))
        self.assertEqual(argv[-1], '90')
        with self.assertRaises(ValueError):
            make_command(Path('/tmp/stage'), Path('/tmp/insertion'),
                         timeout_seconds=0)

    def test_complete_is_diagnostic_only(self):
        self.assertEqual(self.classify(), 'CompleteCanonicalCoreTreeDiagnostic')

    def test_zero_explicit_gauge_terms_do_not_invalidate_a_checked_core(self):
        audit = self.audit()
        self.assertEqual(audit["SMGaugeKineticTermCount"], 0)
        self.assertTrue(audit["CheckLagrangianPassed"])
        self.assertEqual(
            self.classify(audit=audit), 'CompleteCanonicalCoreTreeDiagnostic'
        )

    def test_wl_relies_on_checklagrangian_not_gauge_term_count(self):
        source = (ROOT / 'Lagrangian/RunGeneralDMinus2CanonicalCoreProbe.wl').read_text(
            encoding='utf-8')
        self.assertIn('"SMGaugeKineticTermCount" -> Length[smGaugeKinetic]', source)
        self.assertIn('"ExplicitGaugeKineticPresenceRequired" -> False', source)
        self.assertNotIn('Length[smGaugeKinetic] > 0', source)
        self.assertIn('CheckLagrangian[matchInput]', source)


    def test_elapsed_times_are_quantized_before_json_export(self):
        source = (ROOT / 'Lagrangian/RunGeneralDMinus2CanonicalCoreProbe.wl').read_text(
            encoding='utf-8')
        self.assertIn('D2ElapsedSeconds[seconds_?NumericQ]', source)
        self.assertIn('If[microseconds === 0, 0, N[microseconds/1000000]]', source)
        self.assertIn('"Seconds" -> D2ElapsedSeconds[AbsoluteTime[] - started]', source)
        self.assertIn('"CheckLagrangianSeconds" -> D2ElapsedSeconds[checkElapsed]', source)
        self.assertNotIn('"Seconds" -> AbsoluteTime[] - started', source)

    def test_failure_canonicality_not_certified(self):
        audit = self.audit()
        audit['CheckLagrangianPassed'] = False
        self.assertEqual(self.classify(audit=audit), 'CanonicalCoreNotCertified')

    def test_failure_sm_not_alpha_equivalent(self):
        audit = self.audit()
        audit['SMAlphaComparisonCertified'] = False
        self.assertEqual(self.classify(audit=audit), 'SMAlphaComparisonNotCertified')

    def test_failure_heavy_not_preserved(self):
        audit = self.audit()
        audit['CoreHeavySectorExactlyPreserved'] = False
        self.assertEqual(self.classify(audit=audit), 'CanonicalCoreNotCertified')

    def test_failure_sm_kinetic_not_found(self):
        audit = self.audit()
        audit['SMKineticTermCount'] = 0
        self.assertEqual(self.classify(audit=audit), 'CanonicalCoreNotCertified')

    def test_failure_numeric_match_return(self):
        result = self.result()
        result['Stages'][0]['Head'] = 'Integer'
        self.assertEqual(self.classify(result=result), 'InvalidEFTStageReturn')

    def test_failure_no_match_return_stage(self):
        result = self.result()
        result['Stages'][0]['Status'] = 'InvalidAtomicEFT'
        self.assertEqual(self.classify(result=result), 'FailedTreeProcessingStage')

    def test_failure_incomplete_eft_chain(self):
        result = self.result()
        result['Stages'] = result['Stages'][:1]
        self.assertEqual(self.classify(result=result), 'MissingTreeProcessingStage')

    def test_external_timeout(self):
        self.assertEqual(self.classify(external_timeout=True), 'ExternalTimeout')

    def test_failure_no_full_physical_sm_claim(self):
        result = self.result()
        result['FullPhysicalSMIncluded'] = True
        self.assertEqual(self.classify(result=result), 'InvalidDiagnosticInputContract')

    def test_wl_retains_exact_source_and_validates_match(self):
        wl = (ROOT / 'Lagrangian/RunGeneralDMinus2CanonicalCoreProbe.wl')
        source = wl.read_text(encoding='utf-8')
        for fragment in (
            '"eft1_after_F_realbasis_tree.wxf"',
            'D2AlphaIndexKey[term_]',
            'exactBSM = Expand[exactTree - smStandalone]',
            'smKineticTerms = Select[sourceSM, D2SMKineticQ]',
            'coreTree = Expand[exactBSM + Total[smKineticTerms]]',
            'coreResidual = Expand[coreTree + Total[smDroppedTerms] - exactTree]',
            'coreHeavyPreserved', 'smUnclassifiedCouplingless',
            'CheckLagrangian[matchInput]',
            'NumericQ[result] || AtomQ[result]',
            'D2ContainsFieldQ',
            'D2Step["Match"', 'D2Step["GreensSimplify"',
            'D2Step["EOMSimplify"',
        ):
            self.assertIn(fragment, source)
        self.assertNotIn('matchInput = RelabelIndices[exactBSM]', source)
        self.assertIn('"FullPhysicalSMIncluded" -> False', source)


if __name__ == '__main__':
    unittest.main()
