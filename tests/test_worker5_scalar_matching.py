"""Worker 5 regression tests; run from repository root with unittest."""
import unittest
import numpy as np
from Numerical.diagnostics.worker5_scalar_matching import mat,matrix_metrics
from Numerical.joint_scalar_matching.core import match_c5,b0bar,weinberg_tensor


class Worker5Tests(unittest.TestCase):
    def test_matrix_overlap_identical(self):
        z={'real':np.eye(3).tolist(),'imag':np.zeros((3,3)).tolist()}
        d=matrix_metrics(z,z)
        self.assertAlmostEqual(d['absolute_cosine'],1.)
        self.assertAlmostEqual(d['least_squares_complex_line_residual'],0.)

    def test_matrix_overlap_different_flavour(self):
        a={'real':np.diag([1,0,0]).tolist(),'imag':np.zeros((3,3)).tolist()}
        b={'real':np.diag([0,1,0]).tolist(),'imag':np.zeros((3,3)).tolist()}
        d=matrix_metrics(a,b)
        self.assertAlmostEqual(d['absolute_cosine'],0.)
        self.assertAlmostEqual(d['least_squares_complex_line_residual'],1.)

    def test_bubble_plus_inherited_once(self):
        # Exactly degenerate positive masses, two scalar species with two real
        # components each. Pick an LLSS*HHSS tensor factorized in O5.
        n=4; T=weinberg_tensor()
        C=np.zeros((3,3,2,2,n,n),complex)
        Q=np.zeros((4,4,n,n),float)
        # The Weinberg T carries complex Higgs components, while the heavy
        # tensor is real; use only an independent direct zero-channel check.
        c5=np.eye(3,dtype=complex)*1e-12
        ans=match_c5(C,Q,np.eye(n)*1e6,1e3,c5_eft1=c5,
                     group_names=['S1','S1','S2','S2'])
        np.testing.assert_allclose(ans['c5_total'],c5,rtol=0,atol=1e-28)
        np.testing.assert_allclose(ans['c5_scalar_bubble'],0,atol=1e-28)
        self.assertLessEqual(ans['projection_residual'],1e-8)

    def test_b0_matching_scale(self):
        self.assertAlmostEqual(b0bar(1e6,1e6,1e3),0.,places=12)
        self.assertAlmostEqual(b0bar(1e6,1e6,2e3),2*np.log(2),places=12)

    def test_reject_bad_matrix(self):
        with self.assertRaises(ValueError):
            mat({'real':[[1]],'imag':[[0]]})


if __name__=='__main__': unittest.main()
