"""Deterministic mass-Hessian and provenance tests; no Wolfram dependency."""
import unittest
import numpy as np
from Numerical.experimental.Batch16ScalarMassHessians import (
    audit_mass_spectrum,real_mass_hessian,analytic_mass_spectrum,SPECIAL,DIMS)
from scripts.RunBatch16ScalarMassAll16 import MODELS


def make_input(model,mu=2e5+3e5j,m1=9e6,m2=16e6):
    kind=SPECIAL.get(model,'none')
    names={'self_s1':'muT3GeneralS1Sq',
           'self_s2':'muT3GeneralS2Sq','mixed':'muT3GeneralS1S2'}
    key=names.get(kind)
    return {'status':'CandidateScalarThresholdHandoffAudited',
            'physical_neutrino_mass_prediction':False,'model':model,
            'mu_f_gev':1e5,'mu_s_gev':1e3,
            'scalar_threshold':{
                'diagonal_mass_squared_gev2':{'mS1Sq':m1,'mS2Sq':m2},
                'additional_quadratic_symbols':[key] if key else [],
                'additional_quadratic_values':{key:{'real':mu.real,'imag':mu.imag}} if key else {},
            }}


class MassHessianTests(unittest.TestCase):
    def test_all_16_coverage_and_eight_special(self):
        self.assertEqual(len(MODELS),16)
        self.assertEqual(len(SPECIAL),8)
        for name in MODELS:
            self.assertIn(name[0],DIMS)
            self.assertEqual(audit_mass_spectrum(make_input(name))['model'],name)

    def test_self_singlet_complex_splitting(self):
        r=audit_mass_spectrum(make_input('A_+0',mu=3e5+4e5j))
        spectrum=r['conventions']['rgbeta_potential_plus']['mass_squared_eigenvalues_gev2']
        self.assertTrue(np.allclose(spectrum,[8e6,10e6]+[16e6]*6))

    def test_self_triplet_degeneracy(self):
        r=audit_mass_spectrum(make_input('D_+0',mu=3e5+4e5j))
        spectrum=r['conventions']['rgbeta_potential_plus']['mass_squared_eigenvalues_gev2']
        self.assertTrue(np.allclose(spectrum,[8e6]*3+[10e6]*3+[16e6]*2))

    def test_mixed_epsilon_spectrum_fourfold(self):
        r=audit_mass_spectrum(make_input('B_-1'))
        s=r['conventions']['rgbeta_potential_plus']['mass_squared_eigenvalues_gev2']
        self.assertEqual(len(s),8)
        self.assertTrue(np.allclose(s[:4],s[0]))
        self.assertTrue(np.allclose(s[4:],s[-1]))
        self.assertLess(s[0],9e6)
        self.assertGreater(s[-1],16e6)

    def test_convention_flip_changes_hessian_not_spectrum(self):
        r=audit_mass_spectrum(make_input('C_-1'))
        a=r['conventions']['rgbeta_potential_plus']
        b=r['conventions']['matchete_lint_plus']
        self.assertTrue(np.allclose(a['mass_squared_eigenvalues_gev2'],b['mass_squared_eigenvalues_gev2']))
        self.assertFalse(np.allclose(a['real_mass_hessian_gev2'],b['real_mass_hessian_gev2']))

    def test_zero_mu_recovers_proxy(self):
        for model in MODELS:
            r=audit_mass_spectrum(make_input(model,mu=0j))
            self.assertEqual(r['maximum_shift_from_diagonal_proxy_gev2'],0.)

    def test_detects_tachyon_from_bilinear(self):
        r=audit_mass_spectrum(make_input('A_+0',mu=6e6+0j))
        self.assertFalse(r['conventions']['rgbeta_potential_plus']['tachyon_free_candidate'])
        self.assertFalse(r['conventions']['rgbeta_potential_plus']['eigenmass_threshold_hierarchy_candidate'])

    def test_reject_missing_quadratic(self):
        data=make_input('E_-2')
        data['scalar_threshold']['additional_quadratic_symbols']=[]
        with self.assertRaisesRegex(ValueError,'quadratic inventory mismatch'):
            audit_mass_spectrum(data)

    def test_reject_spurious_quadratic(self):
        data=make_input('A_-4')
        data['scalar_threshold']['additional_quadratic_symbols']=['muT3GeneralS1Sq']
        with self.assertRaisesRegex(ValueError,'quadratic inventory mismatch'):
            audit_mass_spectrum(data)

    def test_no_physical_c5_claim(self):
        r=audit_mass_spectrum(make_input('B_-1'))
        self.assertFalse(r['physical_neutrino_mass_prediction'])
        self.assertFalse(r['matchete_rgbeta_mass_parameter_mapping_certified'])
        self.assertFalse(r['mass_eigenvectors_authoritative_for_wilson_matching'])

if __name__=='__main__':
    unittest.main()
