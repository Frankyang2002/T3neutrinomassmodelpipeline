"""Offline regression checks of Batch17 A(-4) new endpoint scalar masses."""
import unittest
import numpy as np
from Numerical.experimental.Batch17RebuildBackground import (
    scalar_couplings_at_mu_s, expected_mass_at_mu_s, check_native_mass_matrix,
)


def enc(x):
    z=complex(x)
    return {'real':z.real,'imag':z.imag}


def fixture(m1=4e6,m2=9e6):
    endpoint={
        'model':'A_-4','status':'CandidateFThresholdAndEFT1Integrated',
        'final':{'mS1Sq':enc(m1),'mS2Sq':enc(m2),
                 'lambdaH1':enc(0.2),'lambdaT3':enc(.02+.04j),
                 'yu':{'real':np.eye(3).tolist(),'imag':np.zeros((3,3)).tolist()}},
    }
    h=np.diag([m1]*2+[m2]*6)
    mass={
        'model':'A_-4','status':'CandidateNativeQuadraticHessianComputed',
        'quadratic_kind':'none',
        'conventions':{'rgbeta_potential_plus':{
            'tachyon_free_candidate':True,
            'mass_squared_eigenvalues_gev2':np.linalg.eigvalsh(h).tolist(),
            'real_mass_hessian_gev2':h.tolist()}},
    }
    return endpoint,mass


class Batch17RebuildTests(unittest.TestCase):
    def test_complete_scalar_decoding(self):
        ep,_=fixture()
        vals=scalar_couplings_at_mu_s(ep)
        self.assertEqual(vals['mS1Sq'],4e6)
        self.assertEqual(vals['lambdaT3'],{'real':.02,'imag':.04})
        self.assertNotIn('yu',vals)

    def test_source_basis_mass_matrix(self):
        ep,ma=fixture()
        h=check_native_mass_matrix(ma,scalar_couplings_at_mu_s(ep))
        self.assertEqual(h.shape,(8,8))
        self.assertAlmostEqual(h[0,0],4e6)
        self.assertAlmostEqual(h[7,7],9e6)

    def test_six_triplet_real_components(self):
        ep,ma=fixture()
        h=expected_mass_at_mu_s(ma)
        self.assertEqual(np.count_nonzero(np.diag(h)==9e6),6)

    def test_new_endpoint_mass_mismatch_rejected(self):
        ep,ma=fixture()
        ep['final']['mS2Sq']=enc(7e6)
        with self.assertRaisesRegex(ValueError,'actual Batch14 scalar masses disagree'):
            check_native_mass_matrix(ma,scalar_couplings_at_mu_s(ep))

    def test_negative_mass_rejected(self):
        ep,ma=fixture();ep['final']['mS1Sq']=enc(-1.)
        with self.assertRaisesRegex(ValueError,'Invalid muS scalar running masses'):
            check_native_mass_matrix(ma,scalar_couplings_at_mu_s(ep))

    def test_complex_mass_rejected_in_matrix_parity(self):
        ep,ma=fixture();ep['final']['mS1Sq']=enc(4e6+1j)
        with self.assertRaises((ValueError,TypeError)):
            check_native_mass_matrix(ma,scalar_couplings_at_mu_s(ep))

    def test_nonpositive_hessian_rejected(self):
        ep,ma=fixture()
        ma['conventions']['rgbeta_potential_plus']['tachyon_free_candidate']=False
        with self.assertRaisesRegex(ValueError,'No positive'):
            expected_mass_at_mu_s(ma)

    def test_mixed_model_not_reused(self):
        ep,ma=fixture(); ma['quadratic_kind']='mixed'
        with self.assertRaisesRegex(ValueError,'unexpectedly contains quadratic'):
            expected_mass_at_mu_s(ma)

    def test_wrong_endpoint_model_rejected(self):
        ep,ma=fixture(); ep['model']='A_+0'
        with self.assertRaisesRegex(ValueError,'Expected the successful'):
            scalar_couplings_at_mu_s(ep)

    def test_absent_quartic_not_invented(self):
        ep,_=fixture();vals=scalar_couplings_at_mu_s(ep)
        self.assertNotIn('lambdaH2',vals)


if __name__=='__main__':unittest.main()
