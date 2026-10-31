import unittest
from Numerical.experimental.T3BCMajoranaActionBatch import run,one_trial
class TestMajoranaAction(unittest.TestCase):
    def test_many_generic_complex_units(self):
        r=run(64)
        self.assertFalse(r['physical_C5_certified'])
        self.assertLess(max(r['max_covariance_residuals'].values()),2e-12)
        self.assertGreater(r['negative_control_min_residuals']['wrong_y1_matrix_relative'],1e-3)
    def test_independent_seed(self):
        r=one_trial(987612)
        self.assertLess(r['majorana_mass'],2e-12)
        self.assertGreater(r['wrong_y1_matrix_relative'],1e-3)
    def test_native_flags_cannot_certify(self):
        self.assertFalse(run(2)['matchete_rgbeta_absolute_dictionary_certified'])
if __name__=='__main__':unittest.main()
