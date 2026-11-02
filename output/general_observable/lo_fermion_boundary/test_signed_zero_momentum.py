"""Offline regression for signed zero-momentum heavy exchange diagnostic."""
from __future__ import annotations

from fractions import Fraction as Q
import unittest

from signed_zero_momentum_check import G, bil, integrate_out, amplitude_benchmarks


class ExactGrassmannTests(unittest.TestCase):
    def test_anticommutativity_and_nilpotency(self):
        a,b = G.var('a'), G.var('b')
        self.assertEqual(a*a, G())
        self.assertEqual(a*b, -(b*a))

    def test_left_derivative_includes_sign(self):
        a,b = G.var('a'), G.var('b')
        self.assertEqual((a*b).left_derivative('a'),b)
        self.assertEqual((a*b).left_derivative('b'),-a)

    def test_symmetric_weyl_bilinear(self):
        a=[G.var('a1'),G.var('a2')]
        b=[G.var('b1'),G.var('b2')]
        self.assertEqual(bil(a,b),bil(b,a))
        self.assertEqual(bil(a,a),2*a[0]*a[1])
        self.assertEqual(bil(a,b).c[('a1','b2')], Q(1))
        self.assertEqual(bil(a,b).c[('a2','b1')], Q(-1))

    def test_independent_majorana_uv_saddle(self):
        h=[G.var('h1'),G.var('h2')]
        a=[G.var('a1'),G.var('a2')]
        uv=-Q(13)*bil(h,h)/2+Q(5)*bil(h,a)
        effective, _=integrate_out(uv,['h1','h2'])
        self.assertEqual(effective,Q(25,26)*bil(a,a))
        self.assertNotEqual(effective,-Q(25,26)*bil(a,a))

    def test_independent_dirac_uv_saddle(self):
        u=[G.var('u1'),G.var('u2')]
        v=[G.var('v1'),G.var('v2')]
        a=[G.var('a1'),G.var('a2')]
        b=[G.var('b1'),G.var('b2')]
        uv=-Q(13)*bil(u,v)+Q(5)*bil(u,a)-Q(3)*bil(v,b)
        effective,_=integrate_out(uv,['u1','u2','v1','v2'])
        self.assertEqual(effective,-Q(15,13)*bil(a,b))

    def test_full_mixed_channel_and_symmetry_factors(self):
        result=amplitude_benchmarks()
        self.assertTrue(result['zero_momentum_Grassmann_saddle_exact'])
        self.assertEqual(result['majorana_same_channel_symmetry_factor'],'1/2')
        self.assertEqual(result['majorana_cross_channel_symmetry_factor'],'1')
        self.assertEqual(result['dirac_cross_channel_factor'],'1')

    def test_diagnostic_rejects_physical_certification(self):
        self.assertIs(amplitude_benchmarks()['fully_physical_matchete_to_master_normalization_certified'],False)


if __name__=='__main__':
    unittest.main(verbosity=2)
