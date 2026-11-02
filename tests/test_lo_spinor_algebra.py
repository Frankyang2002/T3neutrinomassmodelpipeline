"""Exact, offline spinor and Grassmann checks. Synthetic, not Matchete-certified."""
import unittest
from Numerical.lo_certification.spinor.spinors import G, convention_probe, exchange_probe, wedge


class AlgebraTests(unittest.TestCase):
    def test_4_component_C_and_C_inverse(self):
        r = convention_probe()
        self.assertEqual(r['transpose_C_PL_over_weyl_epsilon'], '+1')
        self.assertEqual(r['transpose_Cinv_PL_over_weyl_epsilon'], '-1')
        self.assertIsNone(r['matchete_gamma_cc_to_master_weyl_factor'])

    def test_Grassmann_exchange(self):
        r = exchange_probe()
        self.assertTrue(r['exact_grassmann_saddle_signs_verified'])
        self.assertEqual(r['majorana_distinct_cross_channels'], 6)
        self.assertIsNone(r['gamma_cc_to_master_weyl_factor'])

    def test_symmetric_Weyl_bilinear(self):
        a = G.var('a'); b = G.var('b'); c = G.var('c'); d = G.var('d')
        self.assertEqual(wedge((a, b), (c, d)), wedge((c, d), (a, b)))
        self.assertEqual(a*a, G.scalar(0))
        self.assertEqual(a*b, -(b*a))

    def test_majorana_half_vs_cross_two_flavour(self):
        a = (G.var('a1'), G.var('a2'))
        b = (G.var('b1'), G.var('b2'))
        j = (2*a[0] + 3*b[0], 2*a[1] + 3*b[1])
        lhs = wedge(j, j)
        rhs = 4*wedge(a, a) + 12*wedge(a, b) + 9*wedge(b, b)
        self.assertEqual(lhs, rhs)

    def test_real_matrices_anticonjugation(self):
        r = convention_probe()
        self.assertFalse(r['matchete_GammaCC_identified_with_this_C'])


if __name__ == '__main__':
    unittest.main()
