"""Offline tests of parameter, mass-parity and inventory guards."""
import unittest
import numpy as np
from Numerical.experimental.Batch17All16Background import (
    scalar_couplings_at_mu_s, verify_quadratic_input, certify_candidate_mass_parity,
)


def enc(x):
    z=complex(x)
    return {'real':z.real,'imag':z.imag}


def endpoint(model='B_-1'):
    return {'model':model,'status':'CandidateFThresholdAndEFT1Integrated',
            'final':{'mS1Sq':enc(5e6),'mS2Sq':enc(8e6),
                     'muT3GeneralS1S2':enc(0.4e6+0.2e6j),
                     'lambdaT3':enc(.2+.1j),
                     'yu':{'real':np.eye(3).tolist(),'imag':np.zeros((3,3)).tolist()}}}


class Tests(unittest.TestCase):
    def test_decode_complex_quartic(self):
        z=scalar_couplings_at_mu_s(endpoint(),'B_-1')
        self.assertEqual(z['lambdaT3'],enc(.2+.1j))

    def test_preserve_mixing(self):
        z=scalar_couplings_at_mu_s(endpoint(),'B_-1')
        self.assertEqual(z['muT3GeneralS1S2'],enc(.4e6+.2e6j))

    def test_skip_yukawa_matrices(self):
        self.assertNotIn('yu',scalar_couplings_at_mu_s(endpoint(),'B_-1'))

    def test_reject_wrong_model(self):
        with self.assertRaisesRegex(ValueError,'incorrect Batch14'):
            scalar_couplings_at_mu_s(endpoint(),'C_-1')

    def test_reject_bad_quadratic(self):
        z=scalar_couplings_at_mu_s(endpoint(),'B_-1')
        audit={'quadratic_kind':'mixed','quadratic_parameter':'muT3GeneralS1S2',
               'quadratic_value':enc(.3e6+.2e6j)}
        with self.assertRaisesRegex(ValueError,'disagree'):
            verify_quadratic_input(audit,z,'B_-1')

    def test_accept_matching_quadratic(self):
        z=scalar_couplings_at_mu_s(endpoint(),'B_-1')
        audit={'quadratic_kind':'mixed','quadratic_parameter':'muT3GeneralS1S2',
               'quadratic_value':enc(.4e6+.2e6j)}
        verify_quadratic_input(audit,z,'B_-1')

    def test_fail_absent_quadratic(self):
        z=scalar_couplings_at_mu_s(endpoint(),'B_-1')
        with self.assertRaisesRegex(ValueError,'absent quadratic'):
            verify_quadratic_input({'quadratic_kind':'mixed','quadratic_parameter':'muMISSING'},z,'B_-1')

    def test_positive_spectrum_parity(self):
        h=np.diag([5e6]*4+[8e6]*4)
        audit={'conventions':{'rgbeta_potential_plus':{
            'mass_squared_eigenvalues_gev2':np.linalg.eigvalsh(h).tolist()}}}
        self.assertLess(certify_candidate_mass_parity(audit,h)['spectrum_relative_max_difference'],1e-12)

    def test_reject_old_spectrum(self):
        h=np.diag([5e6]*4+[8e6]*4)
        audit={'conventions':{'rgbeta_potential_plus':{
            'mass_squared_eigenvalues_gev2':[4e6]*4+[8e6]*4}}}
        with self.assertRaisesRegex(ValueError,'mass-spectrum mismatch'):
            certify_candidate_mass_parity(audit,h)

    def test_reject_bad_real_mass(self):
        ep=endpoint();ep['final']['mS1Sq']=enc(-2e6)
        with self.assertRaisesRegex(ValueError,'real and positive'):
            scalar_couplings_at_mu_s(ep,'B_-1')


if __name__=='__main__':unittest.main()
