"""Local NumPy-only regression tests; no Wolfram/Matchete or saved RG output required."""
import unittest
import numpy as np

from Numerical.experimental import T3BCScalarYThresholdProjection as p


def complex_U():
    A=np.array([[1.1e6,2.0e5+1.8e5j],[2.0e5-1.8e5j,1.4e6]],dtype=complex)
    _,U=np.linalg.eigh(A)
    for k in range(2):
        q=np.argmax(abs(U[:,k]))
        U[:,k]*=np.conj(U[q,k]/abs(U[q,k]))
    assert np.max(abs(U.imag))>0.01
    return U


def random_symmetric_tensor(seed=4):
    r=np.random.default_rng(seed)
    C=(r.normal(size=(3,3,12,12,6))+1j*r.normal(size=(3,3,12,12,6)))
    C=(C+C.swapaxes(0,1))/2
    C=(C+C.swapaxes(2,3))/2
    return C


class ExactYProjectionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.U=complex_U()
        cls.T=p.mixing_real_matrix(cls.U)

    def test_realification_is_orthogonal(self):
        np.testing.assert_allclose(self.T.T@self.T,np.eye(12),atol=1e-13)
        np.testing.assert_allclose(self.T@self.T.T,np.eye(12),atol=1e-13)

    def test_original_complex_field_identity_including_epsilon(self):
        self.assertLess(p.check_field_identity(self.U,self.T),1e-12)
        rng=np.random.default_rng(73)
        h,x,y=[rng.normal(size=2)+1j*rng.normal(size=2) for _ in range(3)]
        vec=p.real_field_vector(h,x,y)
        true_s2=-p.EPS@np.conj(self.U[1,0]*x+self.U[1,1]*y)
        got=(self.T@vec)[8:12]
        want=p.real_field_vector(np.zeros(2),np.zeros(2),true_s2)[8:12]
        np.testing.assert_allclose(got,want,atol=1e-12)

    def test_negative_control_wrong_epsilon_sign(self):
        wrong=self.T.copy()
        wrong[8:12,:]*=-1
        with self.assertRaisesRegex(ValueError,'epsilon/conjugation'):
            p.check_field_identity(self.U,wrong)

    def test_negative_control_drop_complex_mixing_phase(self):
        # Replacing complex U with real(U) is not unitary and cannot preserve fields.
        with self.assertRaisesRegex(ValueError,'unitary'):
            p.mixing_real_matrix(self.U.real)
        # Artificially changing the imaginary component rows spoils the exact identity.
        wrong=self.T.copy();wrong[5,5]+=0.05
        with self.assertRaisesRegex(ValueError,'epsilon/conjugation'):
            p.check_field_identity(self.U,wrong)

    def test_negative_control_scalar_order_swap(self):
        wrong=self.T.copy()
        wrong[[4,5,6,7,8,9,10,11],:]=self.T[[8,9,10,11,4,5,6,7],:]
        with self.assertRaisesRegex(ValueError,'epsilon/conjugation'):
            p.check_field_identity(self.U,wrong)

    def test_su2_and_hypercharge_covariance(self):
        self.assertLess(p.check_gauge_covariance(self.T),1e-12)
        wrong=self.T.copy()
        wrong[8:10,:]*=-1  # not an SU(2)-singlet sign redefinition
        with self.assertRaisesRegex(ValueError,'covariance'):
            p.check_gauge_covariance(wrong)

    def test_wilson_polynomial_transport_all_flavor_pairs(self):
        C=random_symmetric_tensor()
        Cprime=p.rotate_wilson(C,self.T)
        self.assertLess(p.check_tensor_identity(C,Cprime,self.T),1e-12)
        for flavor in range(6):
            self.assertGreater(np.linalg.norm(Cprime[:,:,:,:,flavor]),0)
        np.testing.assert_allclose(Cprime[:,:,:4,:4,:],C[:,:,:4,:4,:],atol=1e-12)

    def test_light_projection_exact_at_Y_zero(self):
        C=random_symmetric_tensor()
        rotated=p.rotate_wilson(C,self.T)
        light=p.light_projection(rotated)
        self.assertEqual(light.shape,(3,3,8,8,6))
        self.assertLess(p.check_light_identity(rotated,light),1e-12)
        self.assertLess(max(p.symmetry_residual(light).values()),1e-12)
        self.assertEqual(set(p.classify_survivors(light)),{'HH','HX','XX'})

    def test_negative_control_wrong_wilson_index_order(self):
        C=random_symmetric_tensor()
        rotated=p.rotate_wilson(C,self.T)
        wrong=np.einsum('ap,ijabf,bq->ijpqf',self.T.T,C,self.T.T,optimize=True)
        with self.assertRaisesRegex(ValueError,'polynomial identity'):
            p.check_tensor_identity(C,wrong,self.T)
        self.assertGreater(np.linalg.norm(rotated-wrong),1e-3)

    def test_zero_y_is_exact_restriction_not_full_matching(self):
        C=random_symmetric_tensor()
        rotated=p.rotate_wilson(C,self.T)
        light=p.light_projection(rotated)
        self.assertGreater(np.linalg.norm(rotated[:,:,8:,:,:]),0)
        self.assertLess(p.check_light_identity(rotated,light),1e-12)
        # Nothing in this module computes a finite threshold correction.
        self.assertFalse(hasattr(p,'finite_scalar_threshold_matching'))

    def test_invalid_shapes_and_nans_fail_closed(self):
        C=random_symmetric_tensor()
        with self.assertRaises(ValueError):p.rotate_wilson(C[:,:,:8,:,:],self.T)
        with self.assertRaises(ValueError):p.mixing_real_matrix(np.eye(3))
        bad=C.copy();bad[0,0,0,0,0]=np.nan
        with self.assertRaises(ValueError):p.rotate_wilson(bad,self.T)
        with self.assertRaises(ValueError):p.light_projection(C[:,:,:8,:8,:])


if __name__=='__main__':unittest.main()
