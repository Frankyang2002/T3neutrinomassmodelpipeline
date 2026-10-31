import unittest
from Numerical.experimental.T3BCMajoranaChiralRotation import check
class TestMajoranaChiral(unittest.TestCase):
    def test_exact_matrix_identity(self):
        for seed in (1,137,992):
            r=check(seed)
            self.assertLess(r['inverse_residual'],1e-12)
            self.assertLess(r['charge_conjugation_flip_residual'],1e-12)
            self.assertFalse(r['native_matchete_dictionary_certified'])
if __name__=='__main__':unittest.main()
