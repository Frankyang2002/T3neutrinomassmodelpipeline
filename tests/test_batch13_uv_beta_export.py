"""Static Batch 13 restricted-Yukawa UV/EFT1 export checks.

Run: python -m unittest discover -s tests -p test_batch13_uv_beta_export.py -v
These checks do not evaluate RGBeta or certify coupling normalisations.
"""
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / "RGE/running/rgbeta/RunGeneralT3QuadraticRGBeta.wl"


class TestBatch13UVExport(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.runner = RUNNER.read_text(encoding="utf-8")

    def test_conjugate_yukawas_disabled_before_build(self):
        flag = 'SetEnvironment["T3_EXPERIMENTAL_GENERAL_FERMIONS" -> "0"]'
        self.assertIn(flag, self.runner)
        self.assertNotIn('SetEnvironment["T3_EXPERIMENTAL_GENERAL_FERMIONS" -> "1"]', self.runner)
        self.assertLess(self.runner.index(flag), self.runner.index('T3RGBetaBuild[dS1,dS2,dF,alpha]'))

    def test_restricted_uv_yukawas_checked(self):
        self.assertIn('restrictedYukawaNames = {"y1", "y2"}', self.runner)
        self.assertIn('Restricted topology Yukawa missing from RGBeta registry', self.runner)
        self.assertIn('Restricted topology Yukawa beta missing from UV export', self.runner)
        self.assertIn('"YukawaSector"->"RestrictedTopologyY1Y2"', self.runner)

    def test_conjugate_yukawa_exclusion_checked(self):
        self.assertIn('excludedYukawaNames = {"y1Conj", "y2Conj"}', self.runner)
        self.assertIn('Conjugate-scalar Yukawa unexpectedly registered', self.runner)
        self.assertIn('Conjugate-scalar Yukawa beta unexpectedly exported', self.runner)
        self.assertIn('"IndependentConjugateScalarYukawasIncluded"->False', self.runner)

    def test_eft1_removes_heavy_fermions(self):
        self.assertIn('Heavy-fermion parameter leaked into EFT1 beta export', self.runner)
        self.assertIn('T3RGBetaBuildEFT1[dS1,dS2,dF,alpha]', self.runner)

    def test_unresolved_beta_diagnostics(self):
        self.assertIn('Unresolved RGBeta beta keys: ', self.runner)
        self.assertIn('"bad_beta_expressions"', self.runner)
        self.assertIn('"quartic_projection_keys"', self.runner)
        self.assertIn('T3QFail["Unresolved beta expression"]', self.runner)
        self.assertIn('Block[{Monitor = (#1 &)},', self.runner)

    def test_quartic_projection_and_quadratics_preserved(self):
        self.assertIn('BetaTerm[Quartic,1]', self.runner)
        self.assertIn('Block[{Monitor = (#1 &)}, BetaTerm[Quartic,1]]', self.runner)
        self.assertIn('Quartic projection diagnostic: head = ', self.runner)
        self.assertIn('base = Join[base,quartics]', self.runner)
        for name in ('muT3GeneralS1Sq', 'muT3GeneralS2Sq', 'muT3GeneralS1S2'):
            self.assertIn(name, self.runner)
        self.assertIn('"physical_mapping_certified"->False', self.runner)


if __name__ == "__main__":
    unittest.main()
