"""Source regression and strict guard tests; not Wolfram runtime tests."""
import unittest
from pathlib import Path
from scripts import RunT3BCYCanonicalReconstruction as p


class TestCanonicalReconstruction(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source=p.SCRIPT.read_text(encoding='utf-8')

    def test_source_absorbs_exact_llss(self):
        self.assertIn('composite=Expand[renorm+Total[majorana]]//RelabelIndices;',self.source)

    def test_does_not_average_nonhermitian(self):
        self.assertNotIn('Hc[scalarEFT]+scalarEFT',self.source)
        self.assertNotIn('Hc[composite]+composite',self.source)

    def test_independent_validation_always_attempted(self):
        self.assertIn('attempted=True;\n fullCheck=CheckAbort[',self.source)
        self.assertIn('independent_matchete_check_passed"->clean',self.source)

    def test_does_not_claim_eomeft_equality(self):
        self.assertIn('"accepted"->(clean && exact)',self.source)
        self.assertIn('"saved_eomeft_algebraic_equality_proven"->False',self.source)

    def test_quartics_preserved(self):
        for term in ('nativeH2["Expressions"]','nativeS2["Expressions"]',
                     'nativeExtra["Expressions"]','nativeTopology["Expression"]'):
            self.assertIn(term,self.source)

    def test_original_tree_is_preserved(self):
        self.assertIn('saved-composite',self.source)
        self.assertIn('full_saved_eft_unchanged_on_disk"->True',self.source)

    def test_two_modes(self):
        self.assertIn('"IndependentCanonicalFullEFTCheckPassedEOMEquivalenceUnproved"',self.source)
        self.assertIn('"YOnlyIndependentCanonicalRawEFT2NotCertified"',self.source)

    def test_fail_closed_fields(self):
        self.assertIn('IndependentCanonicalReassemblyFieldInventoryFailed',self.source)
        self.assertIn('IndependentCanonicalReassemblyCheckNotStable',self.source)
        self.assertIn('IndependentCanonicalCandidateNotEstablished',self.source)

    def test_heavy_only_y(self):
        self.assertIn('Mass->{Light,MSX}',self.source)
        self.assertIn('Mass->{Heavy,MSY}',self.source)
        self.assertIn('"heavy_scalar_fields"->{"NewScalarY"}',self.source)

    def test_hermiticity_gate_preserved(self):
        self.assertIn('SavedEFT1TreeHeavyScalarCheckFailed',self.source)
        self.assertIn('"Hermiticity","CanonicallyNormalized","HeavyMassBasis"',self.source)

    def test_no_new_physical_certification(self):
        self.assertIn('"physical_matching_performed"->False',self.source)
        self.assertIn('"physical_sequential_C5_computed"->False',self.source)

    def test_invalid_model_rejected(self):
        with self.assertRaises(ValueError):p.source_paths('D_-1')


if __name__=='__main__':unittest.main()
