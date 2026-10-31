"""Tests for the grouped heavy-scalar bubble, independent of native running."""
from __future__ import annotations
import unittest
import numpy as np
from Numerical.experimental.T3BCGroupedFiniteScalarBubble import group_bubble,rge_log_check,LOOP
from Numerical.experimental.T3BCScalarBubbleThreshold import b0bar_zero


def synthetic():
    C=np.zeros((3,3,12,12,6),complex)
    q=np.zeros((12,12,12,12),float)
    for i in range(3):
        for j in range(3):
            for f in range(6):
                C[i,j,4,4,f]=(i+j+f+1)*1e-9
                C[i,j,8,8,f]=(i+j+f+2)*1e-9
    from itertools import permutations
    for key,val in [((0,0,4,4),0.17),((0,0,8,8),0.08),((0,0,4,8),0.04)]:
        for k in set(permutations(key)):q[k]=val
    return C,q

class TestBubble(unittest.TestCase):
    def test_zero_sources(self):
        c,q=synthetic()
        d,_=group_bubble(c,np.zeros_like(q),(1e6,1.4e6),1e3)
        self.assertEqual(np.linalg.norm(d),0.)
    def test_shape_and_nonzero(self):
        c,q=synthetic()
        d,x=group_bubble(c,q,(1e6,1.4e6),1e3)
        self.assertEqual(d.shape,(3,3,4,4,6))
        self.assertGreater(np.linalg.norm(d),0)
        self.assertFalse(x['physical_finite_matching_certified'])
    def test_scalar_pair_symmetry(self):
        c,q=synthetic();d,_=group_bubble(c,q,(1e6,1.4e6),1e3)
        np.testing.assert_allclose(d,d.swapaxes(2,3),atol=1e-23)
    def test_flavour_symmetry(self):
        c,q=synthetic();d,_=group_bubble(c,q,(1e6,1.4e6),1e3)
        np.testing.assert_allclose(d,d.swapaxes(0,1),atol=1e-23)
    def test_log_beta_cancellation(self):
        c,q=synthetic();x=rge_log_check(c,q,(1e6,1.4e6),1e3)
        self.assertLess(x['relative_log_derivative_residual'],1e-8)
    def test_frozen_scale_shifts(self):
        c,q=synthetic()
        a,_=group_bubble(c,q,(1e6,1.4e6),5e2)
        b,_=group_bubble(c,q,(1e6,1.4e6),2e3)
        self.assertGreater(np.linalg.norm(a-b),0.)
    def test_zero_scale_rejected(self):
        c,q=synthetic()
        with self.assertRaises(ValueError):group_bubble(c,q,(1e6,1.4e6),0.)
    def test_tachyon_rejected(self):
        c,q=synthetic()
        with self.assertRaises(ValueError):group_bubble(c,q,(-1.,1.4e6),1e3)
    def test_asymmetric_source_rejected(self):
        c,q=synthetic();q[0,1,4,4]+=3.
        with self.assertRaises(ValueError):group_bubble(c,q,(1e6,1.4e6),1e3)
    def test_b0_equal_mass(self):
        self.assertAlmostEqual(b0bar_zero(1e6,1e6,1e6),0.)
    def test_no_second_eft_in_runner(self):
        from pathlib import Path
        runner=Path(__file__).resolve().parents[1]/'scripts/RunT3BCMatchingMegaBatch.py'
        s=runner.read_text(encoding='utf8')
        self.assertIn("'one_intermediate_EFT_only':True",s)
        self.assertIn("'physical_C5_certified':False",s)

if __name__=='__main__':unittest.main()
