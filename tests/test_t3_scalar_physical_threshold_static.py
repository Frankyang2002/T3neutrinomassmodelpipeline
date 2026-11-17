"""Static source-contract checks for the A(-4)/D(-2) scalar continuation.

These verify the currently implemented indirection and safeguards; they do not
run Matchete or establish a physical matched C5.
"""
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
PARENT = ROOT / 'Lagrangian/RunT3FullGeneralTwoThresholdNative.wl'
A_SCRIPT = ROOT / 'Lagrangian/RunT3AMinus4FullGeneralPhysicalScalar.wl'
D_SCRIPT = ROOT / 'Lagrangian/RunT3DMinus2FullGeneralPhysicalScalar.wl'


class TestPhysicalScalarContinuation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.parent = PARENT.read_text(encoding='utf-8')
        cls.a = A_SCRIPT.read_text(encoding='utf-8')
        cls.d = D_SCRIPT.read_text(encoding='utf-8')

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
        for mass in ('MS1', 'MS2'):
            self.assertIn(
                f'Coupling[{mass}, {{}}, 1]] :> Coupling[{mass}, {{}}, 0]',
                self.a,
            )
            self.assertIn(
                f'Coupling[{mass}, {{}}, 0]] :> Coupling[{mass}, {{}}, 1]',
                self.a,
            )
        self.assertIn('t3ADemoted === t3ATree', self.a)
        self.assertIn('t3AOriginalMassCounts === t3APromotedMassCounts', self.a)
        # This stage now uses a checked wrapper and preserves the selected input.
        self.assertIn('t3ACheckOriginal = t3ACheckSafe[t3AMatchInput]', self.a)
        self.assertIn('CheckLagrangian[expr]', self.a)

    def test_a_single_scalar_threshold_and_postprocess(self):
        self.assertIn(
            'ReclassifyT3Fields[t3model, {"S1", "S2"}, {"S1", "S2"}]',
            self.a,
        )
        # The Match call is wrapped by t3AMatchAttempt and accepts a local expr.
        # Check the one-loop selected-input invocation, not an obsolete literal.
        self.assertEqual(
            self.a.count('t3AOneLoopAttempt = t3AMatchAttempt[t3AMatchInput, 1, 3600,'),
            1,
        )
        self.assertIn('Match[expr, EFTOrder -> 5, LoopOrder -> loop]', self.a)
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


if __name__ == '__main__':
    unittest.main()
