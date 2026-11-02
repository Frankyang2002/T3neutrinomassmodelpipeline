"""Independent analytic/artificial fixtures; none are T3 numerical predictions."""
from __future__ import annotations
import math
import unittest
import numpy as np
from Numerical.direct_one_loop.core import (takagi,scalar_real_basis,b0_finite,evaluate,LOOP)


class DirectOneLoopCoreTests(unittest.TestCase):
    def test_takagi_arbitrary_complex_and_degenerate(self):
        rng=np.random.default_rng(171)
        q,_=np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))
        for s in ([1.,1.,3.],[1.,1.,1.],[0.,2.,3.],[1e-8,1.,2.]):
            m=q.conj()@np.diag(s)@q.conj().T
            U,masses,error=takagi(m)
            np.testing.assert_allclose(masses,sorted(s),rtol=2e-8,atol=1e-12)
            np.testing.assert_allclose(U.T@m@U,np.diag(masses),atol=2e-12,rtol=1e-10)
            self.assertLess(error,1e-9)

    def test_b0_degenerate(self):
        x=1e6
        for d in [0.,1e-12,1e-8,1e-5,1e-3]:
            actual=b0_finite(x,x*(1+d),x)
            self.assertAlmostEqual(actual,-d/2+d*d/6,delta=max(1e-15,2*d**3))
        self.assertAlmostEqual(b0_finite(1e6,1e6,1e6),0.)

    def test_ma_split_scalars_exact_and_lnv_limit(self):
        # L = -h nu N phi + h.c., phi=(R+i I)/sqrt2.
        # With V= m0^2|phi|^2 + (delta phi^2+h.c.)/2, mR2=m0^2+delta,
        # mI2=m0^2-delta. This is the exact 1/(32pi^2) convention for h.
        h=.014+0.009j; mf=1000.; m0sq=1.e6; delta=1000.
        t={'model':'synthetic_Ma_unit_test','fermion_mass_gev':[[mf]],
           'scalar_A_gev2':[[m0sq]],'scalar_B_gev2':[[delta]],
           'renormalization_scale_gev':mf,
           'vertices':[{'name':'y','field':np.array([[[h]],[[0]],[[0]]]).tolist()}]}
        r=evaluate(t)
        mr=m0sq+delta;mi=m0sq-delta
        expected=(h*h)*mf*LOOP*(b0_finite(mf**2,mr,mf**2)-b0_finite(mf**2,mi,mf**2))/2
        actual=complex(r['total']['complex_symmetric_matrix_gev']['real'][0][0],
                       r['total']['complex_symmetric_matrix_gev']['imag'][0][0])
        self.assertAlmostEqual(abs(actual-expected),0,delta=1e-14)
        self.assertTrue(r['diagnostics']['uv_finite'])
        t['scalar_B_gev2']=[[0.]]
        zero=evaluate(t)
        self.assertLess(zero['diagnostics']['total_frobenius_gev'],1e-17)

    def test_scalar_complex_mixing_and_reconstruction(self):
        A=np.array([[2.e6,1e4+3e3j],[1e4-3e3j,3.e6]],complex)
        B=np.array([[1500+700j,200j],[200j,0]],complex)
        C,eigen,error=scalar_real_basis(A,B)
        self.assertEqual(C.shape,(2,4))
        self.assertGreater(min(eigen),0.)
        self.assertLess(error,1e-9)

    def test_vertex_pair_sum_and_nontrivial_dirac_chain(self):
        # Two Weyl fermions have a Dirac mass; self-masses vanish.
        # Source tensor has phi1 at first vertex, phi2* at second.
        d=1000; z=np.zeros((3,2,2),complex)
        y=z.copy();y[0,0,0]=.08;z[1,1,1]=.12
        A=np.array([[1.e6,3000],[3000,1.6e6]])
        inp={'model':'synthetic_vectorlike','fermion_mass_gev':[[0,d],[d,0]],
           'scalar_A_gev2':A.tolist(),'scalar_B_gev2':[[0,0],[0,0]],
           'vertices':[{'name':'first','field':y.tolist()},
                       {'name':'second','conjugate':z.tolist()}]}
        report=evaluate(inp)
        full=report['total']['complex_symmetric_matrix_gev']
        np.testing.assert_allclose(full['real'],np.asarray(full['real']).T)
        self.assertTrue(report['diagnostics']['uv_finite'])
        self.assertEqual(len(report['vertex_pair_matrices_gev']),3)
        # No scalar mixing => no lepton-number-violating fermion loop.
        inp['scalar_A_gev2']=[[1e6,0],[0,1.6e6]]
        zreport=evaluate(inp)
        self.assertLess(zreport['diagnostics']['total_frobenius_gev'],1e-15)

    def test_refuse_tachyonic_and_nonsymmetric_mass(self):
        with self.assertRaises(ValueError):
            scalar_real_basis(np.diag([-1.]),np.zeros((1,1)))
        with self.assertRaises(ValueError):
            takagi(np.array([[1.,1.],[0.,1.]]))


if __name__=='__main__':unittest.main()
