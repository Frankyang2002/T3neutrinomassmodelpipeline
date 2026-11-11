import unittest
import numpy as np
from scripts.RunT3Worker1ScalarRobustness import real_fit,observables,compare_cases
class Worker1Tests(unittest.TestCase):
    def test_opposite_sign(self):
        x=np.arange(1,13,dtype=float).reshape(3,4)
        v=real_fit(x,-x)
        self.assertAlmostEqual(v['alpha'],-1)
        self.assertAlmostEqual(v['residual'],0)
        self.assertAlmostEqual(v['raw_relative'],2)
    def test_nonproportional(self):
        v=real_fit([1,0],[0,1]);self.assertAlmostEqual(v['alpha'],0);self.assertAlmostEqual(v['residual'],1)
    def test_majorana_singular_values(self):
        m=np.diag([3.,1.,2.]).astype(complex)
        self.assertEqual(observables(m)['takagi_singular_values'],[1.,2.,3.])
    def test_invalid(self):
        with self.assertRaises(ValueError):real_fit([1,2],[1])
if __name__=='__main__':unittest.main()
