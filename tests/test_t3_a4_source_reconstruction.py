"""Offline checks only: real Matchete equivalence and matching must run locally."""
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "Lagrangian/RunT3AMinus4FullGeneralPhysicalScalar.wl"


class TestA4SourceReconstruction(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.s = SOURCE.read_text(encoding="utf-8")

    def test_canonical_index_graph_permutations(self):
        self.assertIn("t3AGraphIndexKey[term_]", self.s)
        self.assertIn("Tuples[Permutations /@ groups]", self.s)
        self.assertIn('t3APermutationLimit = 10000', self.s)
        self.assertIn('"TooManyDummyPermutations"', self.s)
        self.assertIn('"SM_sector_alpha_graph_equivalent_to_UV_LSM"', self.s)

    def test_source_rebuild_respects_actual_uv_interactions(self):
        self.assertIn('t3ABuild["LSM"] + FreeLag[NewScalar1, NewScalar2]', self.s)
        self.assertIn('t3ASourceScalarInteractions', self.s)
        self.assertIn('t3ASourceFInteractions', self.s)
        self.assertIn('t3ASavedWTerms', self.s)
        self.assertIn('Length[t3ASavedWTerms] =!= 2', self.s)
        self.assertIn('t3AScalarBase + Total[t3ASavedWTerms]', self.s)

    def test_full_graph_equality_before_matching(self):
        self.assertIn('t3AFullGraphEquivalent = t3AGraphComputable &&', self.s)
        self.assertIn('t3ASavedGraphKeys === t3ARebuiltGraphKeys', self.s)
        self.assertIn('If[t3AFullGraphEquivalent,', self.s)
        self.assertIn('t3ACheckReconstructed = t3ACheckSafe[t3AReconstructed]', self.s)
        self.assertIn('If[TrueQ[t3ACheckReconstructed]', self.s)
        self.assertIn('t3ASelectedInput = t3AReconstructed', self.s)
        self.assertIn('t3AMatchInput = t3ASelectedInput;', self.s)
        self.assertEqual(self.s.count('Match[t3AMatchInput, EFTOrder -> 5, LoopOrder -> 1]'), 1)

    def test_reconstructed_eft_is_never_used_without_certificate(self):
        self.assertIn('A_source_reconstruction_audit.json', self.s)
        self.assertIn('A_source_reconstructed_complete_EFT1_DIAGNOSTIC.wxf', self.s)
        self.assertIn('"full_EFT_reconstruction_used_for_matching" -> t3ASelectedFromSource', self.s)
        self.assertIn('If[!TrueQ[t3ACheck],', self.s)
        self.assertIn('"projection_used_for_physical_match" -> False', self.s)

    def test_no_change_to_three_fermions_or_physical_completion(self):
        self.assertIn('"fermion_flavours_in_tree_Wilson_definition" -> 3', self.s)
        self.assertIn('"EFT1_Wilson_RG_included" -> False', self.s)
        self.assertIn('"renormalized_physical_C5_verified" -> False', self.s)
        self.assertNotIn('NF=1', self.s)
        self.assertIn('Exit[5]', self.s)


if __name__ == "__main__":
    unittest.main()
