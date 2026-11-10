"""Offline contracts only; real CG identities / Matchete require Wolfram."""
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
A = ROOT / "Lagrangian/RunT3AMinus4FullGeneralPhysicalScalar.wl"


class TestA4HermitianProbe(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.s = A.read_text(encoding="utf-8")

    def test_both_scalar_masses_remain_reversible(self):
        for name in ("MS1", "MS2"):
            self.assertIn(f"Coupling[{name}, {{}}, 1]] :> Coupling[{name}, {{}}, 0]", self.s)
        self.assertIn("t3ADemoted === t3ATree", self.s)

    def test_reconstructed_full_check_never_skipped_for_mismatch(self):
        self.assertIn("t3ACheckReconstructed = t3ACheckSafe[t3AReconstructed]", self.s)
        self.assertIn('"reconstructed_full_CheckLagrangian_DIAGNOSTIC_ONLY"', self.s)
        self.assertIn('"saved_SM_CG_orientations"', self.s)
        self.assertIn('"reconstructed_SM_CG_orientations"', self.s)

    def test_wilson_hermitian_partner_is_derived_not_assumed(self):
        self.assertIn('t3AWPRCount === 1 && t3AWPLCount === 1', self.s)
        self.assertIn('PlusHc[First[t3AWPRTerms]]', self.s)
        self.assertIn('t3AWHolomorphicGraphPreserved', self.s)
        self.assertIn('t3AWHcCheck = t3ACheckSafe[t3AWHcCandidate]', self.s)
        self.assertIn('"W_PlusHc_pair_matches_saved_tree"', self.s)

    def test_non_equivalent_exploratory_is_explicit_opt_in(self):
        self.assertIn('Environment["T3_A4_EXPLORATORY_SOURCE_MATCH"] === "1"', self.s)
        self.assertIn('t3AExploratoryEnabled &&', self.s)
        self.assertIn('t3AExploratoryQualified && !t3AFullGraphEquivalent', self.s)
        self.assertIn('t3AWHolomorphicGraphPreserved && TrueQ[t3AWHcCheck]', self.s)
        self.assertIn('"A_ExploratoryHermitianScalarSoft_NoEFTBoundaryCertificate"', self.s)
        self.assertIn('"A_exploratory_source_soft_"', self.s)
        self.assertIn('"exploratory_source_theory_used"', self.s)

    def test_one_joint_match_not_hard_double_counted(self):
        self.assertEqual(self.s.count('Match[t3AMatchInput, EFTOrder -> 5, LoopOrder -> 1]'), 1)
        self.assertIn('"EFT1_Wilson_RG_included" -> False', self.s)
        self.assertIn('"renormalized_physical_C5_verified" -> False', self.s)
        self.assertIn('Exit[5]', self.s)
        self.assertNotIn('NF=1', self.s)


if __name__ == "__main__":
    unittest.main()
