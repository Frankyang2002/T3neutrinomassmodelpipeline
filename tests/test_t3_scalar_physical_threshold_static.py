"""Dependency-free regression checks for the repaired A(-4) scalar threshold.

These tests validate source contracts, not Matchete execution or physical C5.
"""
from pathlib import Path
import re
import unittest

ROOT = Path(__file__).resolve().parents[1]
PARENT = ROOT / "Lagrangian/RunT3FullGeneralTwoThresholdNative.wl"
A_SCRIPT = ROOT / "Lagrangian/RunT3AMinus4FullGeneralPhysicalScalar.wl"
D_SCRIPT = ROOT / "Lagrangian/RunT3DMinus2FullGeneralPhysicalScalar.wl"


class TestPhysicalScalarContinuation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.parent = PARENT.read_text(encoding="utf-8")
        cls.a = A_SCRIPT.read_text(encoding="utf-8")
        cls.d = D_SCRIPT.read_text(encoding="utf-8")

    def test_parent_does_not_repeat_d_fermion_matching(self):
        self.assertIn('Get[t3DFile]', self.parent)
        self.assertIn('Get[t3AFile]', self.parent)
        self.assertNotIn('RunT3EFT1WilsonFreshThresholdStage[', self.parent)

    def test_a_fermion_flavour_tensor_preserved(self):
        self.assertIn('t3WCount =!= 2', self.parent)
        self.assertIn('Sum[y1[p,r]*y2[q,r]/MF[r],{r,1,3}]', self.parent)
        self.assertIn('t3ATree', self.a)
        self.assertIn('Coupling[T3EFT1W12', self.parent)
        self.assertNotIn('NF=1', self.a)
        self.assertNotIn('RunSequentialT3Matching[', self.a)

    def test_a_mass_promotion_both_scalars_only(self):
        for m in ('MS1', 'MS2'):
            self.assertIn(f'Coupling[{m}, {{}}, 1]] :> Coupling[{m}, {{}}, 0]', self.a)
            self.assertIn(f'Coupling[{m}, {{}}, 0]] :> Coupling[{m}, {{}}, 1]', self.a)
        self.assertIn('t3ADemoted === t3ATree', self.a)
        self.assertIn('t3AOriginalMassCounts === t3APromotedMassCounts', self.a)
        self.assertIn('CheckLagrangian[t3AMatchInput]', self.a)

    def test_a_single_scalar_threshold_and_postprocess(self):
        self.assertIn('ReclassifyT3Fields[t3model, {"S1", "S2"}, {"S1", "S2"}]', self.a)
        self.assertEqual(self.a.count('Match[t3AMatchInput, EFTOrder -> 5, LoopOrder -> 1]'), 1)
        self.assertIn('GreensSimplify', self.a)
        self.assertIn('EOMSimplify', self.a)
        self.assertIn('EvaluateLoopFunctions', self.a)
        self.assertIn('ReplaceEffectiveCouplings', self.a)
        self.assertIn('ExtractWeinbergCoefficient[t3AOneLoop]', self.a)

    def test_no_false_c5_completion_or_hard_double_counting(self):
        self.assertIn('F_hard_inherited" -> False', self.a)
        self.assertIn('renormalized_physical_C5_verified" -> False', self.a)
        self.assertIn('A_PhysicalScalarSoftOnly_HardInheritanceAndWilsonRGPending', self.a)
        self.assertIn('D_PhysicalScalarSoftOnly_HardInheritanceAndWilsonRGPending', self.d)
        self.assertIn('Exit[5]', self.a)


if __name__ == "__main__":
    unittest.main()
