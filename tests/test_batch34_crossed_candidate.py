"""Physics necessary-condition tests for unverified crossed-Yukawa candidate."""
from __future__ import annotations

import unittest
import numpy as np
from scripts.CheckBatch34CrossedCandidate import (diagnose_candidate, u1_nontrivial_fixture)
from RGE.running.intermediate.GeneralEFT1CrossedWardCandidate import crossed_candidate, candidate_beta
from RGE.running.intermediate.GeneralEFT1SevenTermFlow import TERMS, master_beta


class CrossedWardCandidateTests(unittest.TestCase):
    def test_sm_and_nontrivial_u1_ward(self):
        report=diagnose_candidate()
        self.assertEqual(report["status"],"WARD_DIAGNOSTICS_PASSED_UNVERIFIED_CANDIDATE")
        for row in report["tests"]:
            self.assertLess(row["beta_ward_relative"],2e-11)
            self.assertLess(row["crossed_candidate_ward_abs"],1e-11)

    def test_nontrivial_u1_crossed_is_generated_and_gauge_invariant(self):
        c,p=u1_nontrivial_fixture()
        crossed=crossed_candidate(c,p)
        self.assertGreater(np.linalg.norm(crossed),1e-3)
        self.assertEqual(np.linalg.norm(c[0,0]),0)
        self.assertGreater(np.linalg.norm(crossed[0,0]),1e-3)

    def test_conjugate_anti_linearity(self):
        c,p=u1_nontrivial_fixture()
        c0=crossed_candidate(c,p)
        np.testing.assert_allclose(crossed_candidate(1j*c,p),-1j*c0,atol=1e-12)
        np.testing.assert_allclose(crossed_candidate(.51*c-.33j*c,p),(.51+.33j)*c0,atol=1e-12)

    def test_other_six_terms_remain_identical_to_v1_evaluator(self):
        c,p=u1_nontrivial_fixture()
        before=master_beta(c,p,corrected_crossed=True)
        after=candidate_beta(c,p)
        for name in TERMS:
            if name != "crossed_yukawa":
                np.testing.assert_array_equal(before[name],after[name])
        self.assertGreater(np.linalg.norm(after["crossed_yukawa"]-before["crossed_yukawa"]),1e-3)

    def test_report_never_certifies_physics(self):
        rep=diagnose_candidate()
        self.assertFalse(rep["physical_A_minus4_rge_certified"])
        self.assertFalse(rep["chatA_canonical_boundary_verified"])
        self.assertFalse(rep["chatB_running_general_EFT1_couplings_verified"])


if __name__=="__main__":unittest.main()
