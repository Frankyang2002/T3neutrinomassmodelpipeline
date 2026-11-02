"""Standalone sign-sensitivity tests: no physical convention certification."""
import copy
import unittest

import numpy as np

from scripts.RunT3DirectOneLoopPortalSignSensitivity import make_variant, complex_matrix


class PortalSignSensitivityTests(unittest.TestCase):
    @staticmethod
    def seed():
        matrix=np.array([[10.0-1.,0.2],[0.2,30.0-3.]],complex)
        source={
            'model':'B_-1',
            'scalar_A_gev2':{'real':matrix.real.tolist(),'imag':matrix.imag.tolist()},
            'provenance':{'validated_uv_component_conventions':False},
        }
        frozen={
            'representation':{'d_s1':1,'d_s2':1},
            'base_state':{'ordinary':{'lambdaH1':0.02,'lambdaH2':0.06,
                                      'mS1Sq':10.0,'mS2Sq':30.0}},
            'numerical':{'vev_gev':10.0}}
        identity={'model':'B_-1','scalars':[
            {'scalar':'S1','dimension':1,'full_tensor_identity_in_real_native_span':True,
             'max_abs_four_index_residual':0.,'native_invariant_count':1,
             'native_coupling_coefficients_for_unit_identity':[1.]},
            {'scalar':'S2','dimension':1,'full_tensor_identity_in_real_native_span':True,
             'max_abs_four_index_residual':0.,'native_invariant_count':1,
             'native_coupling_coefficients_for_unit_identity':[1.]}]}
        return source,frozen,identity

    def test_only_diagonal_shifted_and_source_is_unchanged(self):
        source,frozen,identity=self.seed()
        before=copy.deepcopy(source)
        result=make_variant(source,frozen,identity)
        self.assertEqual(source,before)
        np.testing.assert_allclose(complex_matrix(result['scalar_A_gev2']),
                                   np.array([[11.,0.2],[0.2,33.]]))
        self.assertFalse(result['provenance']['validated_uv_component_conventions'])
        self.assertTrue(result['provenance']['portal_sign_sensitivity_only'])

    def test_fails_closed_if_baseline_diagonal_changed(self):
        source,frozen,identity=self.seed()
        source['scalar_A_gev2']['real'][0][0] += 1.
        with self.assertRaisesRegex(ValueError,'baseline'):
            make_variant(source,frozen,identity)

    def test_fails_closed_without_full_tensor_pass(self):
        source,frozen,identity=self.seed()
        identity['scalars'][1]['full_tensor_identity_in_real_native_span']=False
        with self.assertRaisesRegex(ValueError,'identity'):
            make_variant(source,frozen,identity)

    def test_fails_closed_if_certified_input(self):
        source,frozen,identity=self.seed()
        source['provenance']['validated_uv_component_conventions']=True
        with self.assertRaisesRegex(ValueError,'downgrade'):
            make_variant(source,frozen,identity)

if __name__=='__main__':
    unittest.main()
