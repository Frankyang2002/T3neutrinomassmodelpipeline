import unittest
from types import SimpleNamespace
import numpy as np
from Numerical.running.T3ScalarRunningSourceBridge import (
 source_quadratic_tensor,mass_affine,build_from_cached_context,compare_masses)

class WorkerBSourceTests(unittest.TestCase):
    def test_a4_mass_blocks(self):
        m=source_quadratic_tensor('A_-4',{'mS1Sq':4.,'mS2Sq':9.})
        self.assertTrue(np.allclose(np.diag(m)[4:],[4.,4.]+[9.]*6))
    def test_d2_complex_bilinear(self):
        st={'mS1Sq':25.,'mS2Sq':16.,'muT3GeneralS2Sq':{'real':3.,'imag':4.}}
        m=source_quadratic_tensor('D_-2',st)
        np.testing.assert_allclose(m[10:,10:],[[10.,-8.],[-8.,22.]])
        np.testing.assert_allclose(np.linalg.eigvalsh(m[10:,10:]),[6.,26.])
    def test_d2_missing_mu_fails(self):
        with self.assertRaises(ValueError):source_quadratic_tensor('D_-2',{'mS1Sq':1,'mS2Sq':2})
    def test_affine_mass_at_changed_scale(self):
        ref={'mS1Sq':25.,'mS2Sq':16.,'muT3GeneralS2Sq':{'real':3.,'imag':4.}}
        nxt={'mS1Sq':17.,'mS2Sq':12.,'muT3GeneralS2Sq':{'real':-2.,'imag':1.}}
        np.testing.assert_allclose(mass_affine('D_-2',ref).evaluate(nxt),source_quadratic_tensor('D_-2',nxt))
    def test_context_bridge_requires_hhss(self):
        ref={'mS1Sq':25.,'mS2Sq':16.,'q':0.2}
        ctx=SimpleNamespace(model='A_-4',reference=ref,reference_tensor=np.zeros((12,)*4),
                            columns=np.zeros((1,)+(12,)*4),directions=(('q','real'),))
        with self.assertRaises(ValueError): build_from_cached_context('A_-4',1e3,ref,ctx,hhss_source=None)
        bg=build_from_cached_context('A_-4',1e3,ref,ctx,hhss_source=np.zeros((4,4,8,8)))
        np.testing.assert_allclose(bg.eigenvalues_gev2,[16.]*6+[25.]*2)
    def test_comparison(self):
        m=np.diag([0]*4+[4]*2+[9]*6)
        self.assertEqual(compare_masses(m,[9]*6+[4]*2)['maximum_absolute_difference_gev2'],0.)

if __name__=='__main__':unittest.main()
