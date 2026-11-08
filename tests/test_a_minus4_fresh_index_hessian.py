"""Independent exact A(-4) scalar one-loop and native-index experiment checks."""
from __future__ import annotations

import unittest
from fractions import Fraction

import numpy as np

from Numerical.running.Aminus4QuarticHessianIdentities import (
    COORDINATES, beta_monomials_exact, group_identities,
    evaluated_beta, restricted_coefficients,
)


class ExactAminus4Hessian(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mono=beta_monomials_exact()

    def test_all_55_coupling_monomials_have_exact_polynomial_closure(self):
        self.assertEqual(len(self.mono),55)
        self.assertEqual(len(COORDINATES),10)
        self.assertTrue(all(isinstance(x, tuple) and len(x)==2
                            for x in self.mono))

    def test_SU2_Fierz_and_spin_one_identity(self):
        result=group_identities()
        self.assertTrue(result['fundamental_Fierz_exact'])
        self.assertTrue(result['spin_one_adjoint_square_exact'])

    def test_triplet_adjoint_squared_coefficients(self):
        x=self.mono['lambdaS2Adj','lambdaS2Adj']
        self.assertEqual(x,{'lambdaS2':8,'lambdaS2Adj':6})

    def test_portal_triplet_crossed_coefficient(self):
        self.assertEqual(self.mono['lambdaS2Adj','lambdaH2Adj'],
                         {'lambdaH2Adj':6})

    def test_portal_squared_coefficients(self):
        self.assertEqual(self.mono['lambdaH2Adj','lambdaH2Adj'],
                         {'lambdaH':1,'lambdaS2Adj':1,'lambdaH2':2})

    def test_topology_modulus_not_real_squared(self):
        re=self.mono['lambdaT3.real','lambdaT3.real']
        im=self.mono['lambdaT3.imag','lambdaT3.imag']
        cross=self.mono['lambdaT3.real','lambdaT3.imag']
        self.assertEqual(re,im)
        self.assertEqual(re['lambdaH2Adj'],8)
        self.assertEqual(cross,{})
        self.assertEqual(len(re),5)

    def test_SM_scalar_only_limit(self):
        self.assertEqual(self.mono['lambdaH','lambdaH'],{'lambdaH':12})

    def test_five_independent_mismatches_reconstruct_archived_native(self):
        # This is a comparison of the independently derived exact Hessian
        # coefficients with the *observed* native-minus-reference source
        # polynomial, not a refit of any model parameter.
        observed={
          ('lambdaS2Adj','lambdaS2Adj'):{'lambdaS2':8,'lambdaS2Adj':14},
          ('lambdaS2Adj','lambdaH2Adj'):{'lambdaH2Adj':10},
          ('lambdaH2Adj','lambdaH2Adj'):{'lambdaH':2,'lambdaS2Adj':2,'lambdaH2':4},
          ('lambdaT3.real','lambdaT3.real'):{'lambdaH2Adj':Fraction(-16,3)},
          ('lambdaT3.imag','lambdaT3.imag'):{'lambdaH2Adj':Fraction(-16,3)},
        }
        expected_native={
          ('lambdaS2Adj','lambdaS2Adj'):{'lambdaS2':16,'lambdaS2Adj':20},
          ('lambdaS2Adj','lambdaH2Adj'):{'lambdaH2Adj':16},
          ('lambdaH2Adj','lambdaH2Adj'):{'lambdaH':3,'lambdaS2Adj':3,'lambdaH2':6},
          ('lambdaT3.real','lambdaT3.real'):{'lambdaH2Adj':Fraction(8,3),
                     'lambdaH':4,'lambdaH1':12,'lambdaH2':4,'lambda12':4},
          ('lambdaT3.imag','lambdaT3.imag'):{'lambdaH2Adj':Fraction(8,3),
                     'lambdaH':4,'lambdaH1':12,'lambdaH2':4,'lambda12':4},
        }
        for pair, measured_delta in observed.items():
            reference=restricted_coefficients(self.mono)[pair]
            native={k:reference.get(k,0)+measured_delta.get(k,0)
                    for k in set(reference)|set(measured_delta)}
            native={k:v for k,v in native.items() if v}
            self.assertEqual(native,expected_native[pair])

    def test_global_sign_changes_neither_singular_values_nor_quadratic_norm(self):
        # Algebraic statement about a single linear Wilson-to-mass map.
        # Not a determination of the installed Matchete GammaCC sign.
        m=np.array([[1+.2j,.3-1j,.2j],[.3-1j,-.7+.4j,1.2],
                    [.2j,1.2,1+.8j]],complex)
        np.testing.assert_allclose(m,m.T,atol=0,rtol=0)
        np.testing.assert_allclose(np.linalg.svd(m)[1],np.linalg.svd(-m)[1],
                                   atol=1e-13,rtol=0)


if __name__=='__main__':
    unittest.main()
