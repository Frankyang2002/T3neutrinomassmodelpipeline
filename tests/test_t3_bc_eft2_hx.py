"""Standalone EFT2 tests. Run from repository root with unittest discover."""
import sys
import types
import unittest
from unittest.mock import patch
import numpy as np

from Numerical.experimental import T3BCEFT2HX as e


def raw_sm():
    return types.SimpleNamespace(gY=0.36,g2=0.64,g3=1.1,
        ye=np.zeros((3,3),complex),yu=np.zeros((3,3),complex),
        yd=np.zeros((3,3),complex))


def dummy_beta(**kwargs):
    y=np.zeros((3,3),complex)
    return (41/6*kwargs['gY']**3,-19/6*kwargs['g2_value']**3,
            -7*kwargs['g3']**3, 12*kwargs['lambdaH_value']**2,
            y,y,y,np.zeros((3,3),complex))


def stub_modules():
    rge=types.ModuleType('RGE');running=types.ModuleType('RGE.running')
    wb=types.ModuleType('RGE.running.weinberg')
    child=types.ModuleType('RGE.running.weinberg.WeinbergRGE')
    child.sm_weinberg_beta=dummy_beta
    return {x.__name__:x for x in (rge,running,wb,child)}


class EFT2AlgebraTest(unittest.TestCase):
    def test_six_quartic_tensors_linearly_independent(self):
        B=e.template_quartics().reshape(6,-1)
        self.assertEqual(np.linalg.matrix_rank(B),6)

    def test_full_tensor_roundtrip_complex_lambda5(self):
        coefficients=np.array([.27,.19,.11,-.08,.004,-.003])
        extracted,residual=e.fit_quartics(e.quartic_tensor(coefficients))
        np.testing.assert_allclose(extracted,coefficients,atol=2e-14)
        self.assertLess(residual,1e-13)

    def test_exact_scalar_polynomial_normalization(self):
        rng=np.random.default_rng(27)
        p=rng.normal(size=8)
        h=(p[:4:2]+1j*p[1:4:2])/np.sqrt(2)
        x=(p[4::2]+1j*p[5::2])/np.sqrt(2)
        eta=np.array([x[1].conjugate(),-x[0].conjugate()])
        lam=np.array([.27,.19,.11,-.08,.004,-.003])
        z=np.vdot(h,eta)
        expected=(.5*lam[0]*np.vdot(h,h).real**2
            +.5*lam[1]*np.vdot(eta,eta).real**2
            +lam[2]*np.vdot(h,h).real*np.vdot(eta,eta).real
            +lam[3]*abs(z)**2+np.real((lam[4]+1j*lam[5])*z*z))
        q=e.quartic_tensor(lam)
        result=np.einsum('abcd,a,b,c,d',q,p,p,p,p)/24
        self.assertAlmostEqual(expected,result,places=12)

    def test_wrong_lambda5_imag_sign_rejected(self):
        lam=np.array([.27,.19,.11,-.08,.004,-.003])
        bad=e.quartic_tensor(lam)
        wrong=e.quartic_tensor(lam*np.array([1,1,1,1,1,-1]))
        self.assertGreater(np.linalg.norm(bad-wrong),1e-3)

    def test_forbidden_z2_odd_operator_detected(self):
        q=e.quartic_tensor(np.array([.25,.19,.1,-.05,.006,.002]))
        q[0,0,0,4]+=1e-2
        with self.assertRaisesRegex(ValueError,'not the inert-doublet'):
            e.fit_quartics(q)

    def test_exchange_symmetry_of_template(self):
        q=e.quartic_tensor(np.array([.25,.19,.1,-.05,.006,.002]))
        self.assertLess(np.linalg.norm(q-q.transpose((2,1,0,3))),1e-14)
        self.assertLess(np.linalg.norm(q-q.transpose((0,1,3,2))),1e-14)

    def test_inert_potential_is_gauge_covariant(self):
        from Numerical.experimental.T3BCScalarYThresholdProjection import real_doublet
        rng=np.random.default_rng(5)
        v=rng.normal(size=8)
        q=e.quartic_tensor(np.array([.25,.19,.1,-.05,.006,.002]))
        quartic=lambda u: np.einsum('abcd,a,b,c,d',q,u,u,u,u)/24
        for _ in range(7):
            a=rng.normal(size=4);a/=np.linalg.norm(a)
            g=np.array([[a[0]+1j*a[1],a[2]+1j*a[3]],
                        [-a[2]+1j*a[3],a[0]-1j*a[1]]])
            phase=rng.uniform(-np.pi,np.pi)
            U=np.zeros((8,8))
            U[:4,:4]=real_doublet(g*np.exp(1j*phase/2))
            U[4:,4:]=real_doublet(g*np.exp(-1j*phase/2))
            self.assertAlmostEqual(quartic(U@v),quartic(v),places=11)

    def test_invalid_extra_quartic_tensor(self):
        q=e.quartic_tensor(np.ones(6));bad=q.copy();bad[0,1,2,3]=np.nan
        with self.assertRaises(ValueError):e.fit_quartics(bad)

    def test_sm_running_state_roundtrip(self):
        s=raw_sm()
        q=np.array([.25,.19,.1,-.05,.006,.002])
        v=e.pack_running(s,q,1.06e6,-7200.)
        result=e.unpack_running(v)
        self.assertEqual(v.shape,(65,))
        np.testing.assert_array_equal(result['quartics'],q)
        self.assertEqual(result['mx2'],1.06e6)

    def test_inert_gauge_beta_coefficients(self):
        s=raw_sm();q=np.array([.25,.19,.1,-.05,.006,.002])
        v=e.pack_running(s,q,1.06e6,-7200.)
        with patch.dict(sys.modules,stub_modules()):
            beta=e.beta_running(np.log(1150),v)*e.LOOP
        self.assertAlmostEqual(beta[0],(41/6+1/6)*s.gY**3,places=12)
        self.assertAlmostEqual(beta[1],(-19/6+1/6)*s.g2**3,places=12)
        self.assertAlmostEqual(beta[57],12*.25**2+4*.1**2+4*.1*(-.05)
                  +2*(-.05)**2+2*(.006**2+.002**2),places=12)

    def test_complex_lambda5_zero_remains_zero(self):
        s=raw_sm();q=np.array([.25,.19,.1,-.05,0.,0.])
        v=e.pack_running(s,q,1.06e6,-7200.)
        with patch.dict(sys.modules,stub_modules()):
            beta=e.beta_running(np.log(1150),v)
        self.assertEqual(beta[61],0)
        self.assertEqual(beta[62],0)

    def test_self_consistent_X_threshold_not_full_EFT1(self):
        s=raw_sm();q=np.array([.25,.19,.1,-.05,.006,.002])
        v=e.pack_running(s,q,1.06e6,-7200.)
        with patch.dict(sys.modules,stub_modules()):
            mux,yx,sol=e.solve_threshold(1187.5,v)
        self.assertGreater(mux,900)
        self.assertLess(mux,1100)
        self.assertLess(abs(np.sqrt(yx[63])/mux-1),1e-10)

    def test_Wilson_zero_coupling_reduction(self):
        rng=np.random.default_rng(3)
        c=rng.normal(size=(3,3,8,8,6))+1j*rng.normal(size=(3,3,8,8,6))
        c=(c+c.swapaxes(0,1))/2;c=(c+c.swapaxes(2,3))/2
        z=np.zeros
        k=dict(y=z((3,3,8),complex),quartic=z((8,8,8,8)),
            scalarg=z((8,8,8,8)),mixed=z((3,3,8,8),complex),
            gammas=z((8,8),complex),gammaf=z((3,3),complex))
        self.assertLess(np.linalg.norm(e.rhs_master_eft2(c,k)),1e-12)
        self.assertEqual(e.rhs_master_eft2(c,k).shape,c.shape)

    def test_eft2_wilson_rhs_preserves_exchange_symmetries(self):
        rng=np.random.default_rng(5)
        c=rng.normal(size=(3,3,8,8,6))
        c=(c+c.swapaxes(0,1))/2;c=(c+c.swapaxes(2,3))/2
        k=dict(y=np.zeros((3,3,8),complex),quartic=e.quartic_tensor(np.ones(6)),
            scalarg=np.zeros((8,8,8,8)),mixed=np.zeros((3,3,8,8),complex),
            gammas=np.zeros((8,8),complex),gammaf=np.zeros((3,3),complex))
        beta=e.rhs_master_eft2(c,k)
        self.assertLess(np.linalg.norm(beta-beta.swapaxes(0,1)),1e-10)
        self.assertLess(np.linalg.norm(beta-beta.swapaxes(2,3)),1e-10)


if __name__=='__main__':unittest.main()
