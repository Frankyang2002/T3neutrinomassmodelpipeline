"""Standalone tensor and loop-integral tests (no saved Windows study needed)."""
from __future__ import annotations
import importlib.util
from pathlib import Path
import unittest
import numpy as np

SOURCE=Path(__file__).resolve().parents[1]/'Numerical/experimental/T3BCScalarBubbleThreshold.py'
spec=importlib.util.spec_from_file_location('T3BCScalarBubbleThresholdTest',SOURCE)
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)

class TestScalarBubble(unittest.TestCase):
    def test_equal(self):
        self.assertAlmostEqual(mod.b0bar_zero(4.,4.,4.),0.,places=12)
    def test_massless(self):
        self.assertAlmostEqual(mod.b0bar_zero(4.,0.,4.),1.,places=12)
    def test_symmetry(self):
        self.assertAlmostEqual(mod.b0bar_zero(4.,3.,2.),mod.b0bar_zero(3.,4.,2.),places=12)
    def test_scale_derivative(self):
        h=1e-5; a,b=4.,2.75
        up=mod.b0bar_zero(a,b,np.exp(2*h))
        down=mod.b0bar_zero(a,b,np.exp(-2*h))
        self.assertAlmostEqual((up-down)/(2*h),2.,places=8)
    def test_equal_limit(self):
        val=mod.b0bar_zero(4.+1e-9,4.-1e-9,4.)
        self.assertAlmostEqual(val,0.,places=8)
    def test_unequal_number(self):
        r=.76
        expected=1+r/(1-r)*np.log(r)
        self.assertAlmostEqual(mod.b0bar_zero(1.,r,1.),expected,places=12)
    def test_invalid(self):
        with self.assertRaises(ValueError):mod.b0bar_zero(-1.,2.,1.)
    def toy(self):
        rng=np.random.default_rng(33)
        a=rng.normal(size=(3,3,12,12,6))+1j*rng.normal(size=(3,3,12,12,6))
        a=.25*(a+a.swapaxes(0,1)+a.swapaxes(2,3)+a.swapaxes(0,1).swapaxes(2,3))
        q=rng.normal(size=(12,)*4)
        q=sum(q.transpose(p) for p in __import__('itertools').permutations(range(4)))/24.
        return a,q
    def test_complete_partition(self):
        c,q=self.toy();result=mod.partition_beta(c,q)
        expected=np.einsum('abcd,ijcdf->ijabf',q[:8,:8],c,optimize=True)
        light=np.einsum('abcd,ijcdf->ijabf',q[:8,:8,:8,:8],c[:,:,:8,:8,:],optimize=True)
        self.assertTrue(np.allclose(sum(result.values()),expected-light,atol=1e-12))
    def test_finite_scheme_difference(self):
        c,q=self.toy();full,fd=mod.bubble_boundary(c,q,1.2e3,1.03e6,mode='full-mass')
        hard,hd=mod.bubble_boundary(c,q,1.2e3,1.03e6,mode='hard-massless-X')
        self.assertLess(fd['B0_YX'],hd['B0_YX'])
        self.assertGreater(np.linalg.norm(full-hard),0)
    def test_at_x_quartic_only(self):
        self.assertAlmostEqual(mod.x_equal_mass_bubble_at_threshold(1035.),0.,places=12)
    def test_log_coefficient_beta(self):
        c,q=self.toy();y=1200.;x2=1000.**2;h=1e-4
        b=mod.partition_beta(c,q)
        # Differentiation holding the masses fixed: d delta/dln(mu)=-beta_heavy/LOOP.
        # At variable mu, one must vary subtraction scale only, not Y crossing.
        def at_scale(s):
            b0xy=mod.b0bar_zero(y*y,x2,s*s)
            b0hy=mod.b0bar_zero(y*y,0,s*s)
            b0yy=mod.b0bar_zero(y*y,y*y,s*s)
            return -(b0xy*b['XY']+b0hy*b['HY']+b0yy*b['YY'])/(2*mod.LOOP)
        d=(at_scale(y*np.exp(h))-at_scale(y*np.exp(-h)))/(2*h)
        ref=-sum(b.values())/mod.LOOP
        self.assertTrue(np.allclose(d,ref,atol=1e-11,rtol=1e-9))
    def test_zero_wilson(self):
        c,q=self.toy();delta,_=mod.bubble_boundary(c*0,q,1200.,1e6,mode='full-mass')
        self.assertEqual(np.linalg.norm(delta),0)
    def test_symmetric_output(self):
        c,q=self.toy();d,_=mod.bubble_boundary(c,q,1200.,1e6,mode='full-mass')
        self.assertTrue(np.allclose(d,d.swapaxes(0,1)))
        self.assertTrue(np.allclose(d,d.swapaxes(2,3)))

if __name__=='__main__':unittest.main()
