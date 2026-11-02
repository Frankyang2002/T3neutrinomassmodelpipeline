"""Unit tests for OFF-SHELL native portal identity and hermitian symmetrization."""
from __future__ import annotations
import unittest
import numpy as np

from scripts.ResolveT3DirectUVPortalIdentity import (
    fit_native_identity, hermitian_portal, identity_tensor,
)


class DirectUVPortalIdentityTests(unittest.TestCase):
    def test_singlet_native_identity(self):
        t=identity_tensor(1)
        x=fit_native_identity([t])
        self.assertTrue(x['full_tensor_identity_in_real_native_span'])
        self.assertAlmostEqual(x['native_coupling_coefficients_for_unit_identity'][0],1.)

    def test_two_native_invariants_need_not_equal_identity_one_by_one(self):
        dim=2
        t=identity_tensor(dim)
        exchange=np.zeros_like(t)
        for i in range(2):
            for j in range(2):
                for a in range(2):
                    for b in range(2):
                        exchange[i,j,a,b]=int(i==b and j==a)
        result=fit_native_identity([exchange,t])
        self.assertTrue(result['full_tensor_identity_in_real_native_span'])
        self.assertAlmostEqual(result['native_coupling_coefficients_for_unit_identity'][0],0.)
        self.assertAlmostEqual(result['native_coupling_coefficients_for_unit_identity'][1],1.)
        self.assertAlmostEqual(result['max_abs_four_index_residual'],0.,places=12)

    def test_linear_combination_of_nonidentity_terms(self):
        t=identity_tensor(3)
        x=np.zeros_like(t)
        x[0,0,0,0]=1
        x[1,1,2,2]=-1
        r=fit_native_identity([t+x,x])
        self.assertTrue(r['full_tensor_identity_in_real_native_span'])
        np.testing.assert_allclose(r['native_coupling_coefficients_for_unit_identity'],[1.,-1.],atol=1e-12)

    def test_nonspanning_basis_blocks(self):
        x=np.zeros_like(identity_tensor(3))
        x[0,0,0,0]=1
        r=fit_native_identity([x])
        self.assertFalse(r['full_tensor_identity_in_real_native_span'])
        self.assertGreater(r['max_abs_four_index_residual'],0.5)

    def test_real_coefficients_only_and_plus_hc(self):
        t=identity_tensor(2)
        phase=1j*t
        self.assertAlmostEqual(np.max(abs(hermitian_portal(phase))),0.)
        r=fit_native_identity([phase])
        self.assertFalse(r['full_tensor_identity_in_real_native_span'])
        self.assertEqual(r['native_invariant_linear_rank'],0)

    def test_bad_shape(self):
        with self.assertRaisesRegex(ValueError, 'shape'):
            hermitian_portal(np.eye(3))

if __name__=='__main__':
    unittest.main()
