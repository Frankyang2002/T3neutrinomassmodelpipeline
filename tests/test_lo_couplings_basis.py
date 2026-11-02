"""Exact real-scalar transformations, phase-sensitive mass Hessians, and factorials."""
import sys
from pathlib import Path
import unittest
import numpy as np
import sympy as sp
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from Numerical.lo_certification.couplings.basis import (
    orthogonal_matrix, mass2_matrix, transport_mass2, transport_hhss,
    transport_llss, relative_basis_rotation, self_bilinear_hessian,
    checked_eigenbasis)


class TestRealTensorGeometry(unittest.TestCase):
    def setUp(self):
        rng=np.random.default_rng(1209)
        a=rng.normal(size=(8,8)); self.rf=np.linalg.qr(a)[0]
        a=rng.normal(size=(8,8)); self.rs=np.linalg.qr(a)[0]
        self.mass=np.diag(np.arange(1,9)*10.0)
        a=rng.normal(size=(4,4,8,8))
        self.quartic=(a+a.swapaxes(0,1)+a.swapaxes(2,3)+a.transpose(1,0,3,2))/4
        self.c=rng.normal(size=(3,3,2,2,8,8))+1j*rng.normal(size=(3,3,2,2,8,8))

    def test_relative_rotation_at_distinct_scales(self):
        rel=relative_basis_rotation(self.rf,self.rs)
        np.testing.assert_allclose(self.rf@rel,self.rs,rtol=1e-12,atol=1e-12)
        self.assertGreater(np.linalg.norm(rel-np.eye(8)),0.3)

    def test_mass_covariance_and_eigenspectrum(self):
        f=transport_mass2(self.mass,self.rf)
        evals,_=checked_eigenbasis(f)
        np.testing.assert_allclose(np.sort(evals),np.diag(self.mass),rtol=1e-12)
        np.testing.assert_allclose(transport_mass2(f,self.rf.T),self.mass,atol=1e-12)

    def test_quartic_and_wilson_contraction_invariance(self):
        q=transport_hhss(self.quartic,self.rf)
        c=transport_llss(self.c,self.rf)
        before=np.einsum('fgijab,ABab->fgijAB',self.c,self.quartic,optimize=True)
        after=np.einsum('fgijab,ABab->fgijAB',c,q,optimize=True)
        np.testing.assert_allclose(before,after,rtol=1e-12,atol=1e-11)
        np.testing.assert_allclose(transport_hhss(q,self.rf.T),self.quartic,atol=1e-11)

    def test_independent_quartics_add_not_overwrite(self):
        ch=np.stack((self.quartic,self.quartic*0.01),axis=0)
        rotated=transport_hhss(ch,self.rf)
        np.testing.assert_allclose(rotated.sum(0),transport_hhss(ch.sum(0),self.rf),atol=1e-11)
        self.assertGreater(np.linalg.norm(rotated[1]),0.)

    def test_self_mass_complex_phase_and_no_false_diagonalization(self):
        m=self_bilinear_hessian(3,100.,3.+4.j)
        np.testing.assert_allclose(np.sort(np.linalg.eigvalsh(m)),[90.]*3+[110.]*3)
        self.assertNotEqual(float(m[0,3]),0.)
        np.testing.assert_allclose(m[0,3],-8.)

    def test_mass_block_nonpositive_is_rejected(self):
        with self.assertRaisesRegex(ValueError,'nonpositive'):
            mass2_matrix(np.diag([1.,-1.]),2)

    def test_orthogonality_is_not_assumed(self):
        with self.assertRaisesRegex(ValueError,'orthogonal'):
            orthogonal_matrix(np.ones((3,3)),3)

    def test_real_component_fourth_derivative_uses_factorial_4(self):
        h,s,g=sp.symbols('h s g',real=True)
        # V= g h^2 s^2/4; d^4V/dh^2 ds^2=g, not g/4 or 6g.
        v=g*h*h*s*s/sp.Integer(4)
        self.assertEqual(sp.diff(v,h,2,s,2),g)
        # For four distinct real legs, 4! permutations give no extra factor.
        x=sp.symbols('x:4',real=True)
        quart=sp.Symbol('lam')*sp.prod(x)
        self.assertEqual(sp.diff(quart,*x),sp.Symbol('lam'))


if __name__=='__main__':
    unittest.main()
