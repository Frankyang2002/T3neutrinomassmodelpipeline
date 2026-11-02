"""Source B-1 lambdaT3 numerical-pilot tests; no physical certification."""
from __future__ import annotations
import unittest
import importlib.util
import math
import sys
import types
from pathlib import Path
from unittest.mock import patch
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
SPEC=importlib.util.spec_from_file_location('b1_source_pilot',ROOT/'RGE/running/intermediate/B1PreliminaryNumerics.py')
b1=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(b1)

class RealBasisTests(unittest.TestCase):
    def test_identity_still_contains_epsilon_and_conjugation(self):
        R=b1.real_lift_b1(np.eye(2))
        # X0 real = S10 real; Y0 real maps -eps conjugate Y -> S2_1R
        self.assertAlmostEqual(R[0,0],1.)
        self.assertAlmostEqual(R[6,4],1.)
        # Y0 imaginary maps to -S2_1 imaginary (complex conjugation)
        self.assertAlmostEqual(R[7,5],-1.)
        self.assertAlmostEqual(np.linalg.norm(R.T@R-np.eye(8)),0.,places=13)

    def test_complex_mixed_doublet_rotation_is_orthogonal(self):
        t=.35; phase=.41
        U=np.array([[np.cos(t), np.sin(t)*np.exp(1j*phase)],
                    [-np.sin(t)*np.exp(-1j*phase),np.cos(t)]])
        R=b1.real_lift_b1(U)
        self.assertTrue(np.allclose(R.T@R,np.eye(8),atol=1e-14))
        x=np.random.default_rng(1).normal(size=8)
        z=R@x
        X=np.array([x[0]+1j*x[1],x[2]+1j*x[3]])
        Y=np.array([x[4]+1j*x[5],x[6]+1j*x[7]])
        eps=np.array([[0,1],[-1,0]])
        s1=U[0,0]*X+U[0,1]*Y
        s2=-eps@(U[1,0]*X+U[1,1]*Y).conj()
        expected=np.ravel(np.stack([np.r_[s1,s2].real,np.r_[s1,s2].imag],axis=1))
        np.testing.assert_allclose(z,expected,atol=1e-14)

    def test_reject_nonunitary(self):
        with self.assertRaises(b1.PilotError): b1.real_lift_b1(2*np.eye(2))

    def test_mass_eigenmodes_after_muF_to_muS_rotation(self):
        U=np.array([[np.cos(.23),np.sin(.23)],[-np.sin(.23),np.cos(.23)]])
        mass=b1.mass_original_from_eigenbasis(U,[1.7e6,2.8e6])
        Rf=b1.real_lift_b1(np.eye(2))
        inF=Rf.T@mass@Rf
        np.testing.assert_allclose(np.sort(np.linalg.eigvalsh(inF)),
                                   np.sort(np.repeat([1.7e6,2.8e6],4)),rtol=1e-13)
        self.assertGreater(np.linalg.norm(inF-np.diag(np.diag(inF))),1000.)

    def test_quartic_tensor_transport_and_channel_separation(self):
        rng=np.random.default_rng(4)
        q=rng.normal(size=(2,4,4,8,8))
        q=(q+q.swapaxes(1,2))/2
        q=(q+q.swapaxes(3,4))/2
        U=np.array([[np.cos(.15),np.sin(.15)],[-np.sin(.15),np.cos(.15)]])
        R=b1.real_lift_b1(U)
        got=b1.q_in_fermion_basis(q,R)
        self.assertEqual(got.shape,(2,4,4,8,8))
        np.testing.assert_allclose(got.sum(0),np.einsum('ABab,ai,bj->ABij',q.sum(0),R,R))
        self.assertAlmostEqual(np.linalg.norm(q),np.linalg.norm(got),places=10)

    def test_quartic_rejects_asymmetric_source(self):
        q=np.zeros((1,4,4,8,8));q[0,0,2,0,1]=1
        with self.assertRaises(b1.PilotError):b1.q_in_fermion_basis(q,np.eye(8))

    def test_complex_symmetric_serialization(self):
        a=np.eye(3,dtype=complex)*(1+2j)
        encoded=b1.cjson(a)
        np.testing.assert_allclose(np.asarray(encoded['real'])+1j*np.asarray(encoded['imag']),a)


class MinimalCoreMockTests(unittest.TestCase):
    def core(self, residual=0):
        mod=types.ModuleType('Numerical.joint_scalar_matching.core')
        def check_tensors(c,q,n):
            if c.shape != (3,3,2,2,n,n) or q.shape != (4,4,n,n):
                raise ValueError('wrong shapes')
            return c,q
        def project_weinberg(amp,strict=False):
            return (4*amp[:,:,0,0,2,2],residual)
        def bubble(c,q,mass2,mu,reference_scale=None):
            beta=np.einsum('fgijab,ABab->fgijAB',c,q)
            m2=float(mass2[0,0]); b0=-math.log(m2/(mu*mu))
            return {'amplitude':-beta*b0/(2*16*math.pi**2)}
        mod.check_tensors=check_tensors
        mod.project_weinberg=project_weinberg
        mod.scalar_bubble_matching=bubble
        return mod

    def fixtures(self):
        C=np.zeros((3,3,2,2,8,8),complex)
        C[0,0,0,0,0,0]=2e-11+1e-11j
        q=np.zeros((4,4,8,8),float);q[2,2,0,0]=.13
        mass=np.eye(8)*1.7e6
        return C,q,mass

    def test_fixed_order_matching_scale_cancellation(self):
        C,q,mass=self.fixtures()
        with patch.dict(sys.modules,{'Numerical.joint_scalar_matching.core':self.core()}):
            result=b1.evaluate_response(C,q,mass,1e5,1e3)
        self.assertLess(result['matching_scale_cancellation_relative'],1e-11)
        self.assertTrue(result['full_SU2_Weinberg_projection_pass'])
        self.assertEqual(len(result['neutrino_mass_singular_values_ev_per_abs_NGamma']),3)
        self.assertGreater(max(result['neutrino_mass_singular_values_ev_per_abs_NGamma']),0)
        self.assertFalse(result['physical_C5_certified'])
        self.assertEqual(result['joint_to_master_Weinberg_coefficient_factor'],2)
        self.assertAlmostEqual(result['beta_C5_master_16pi2_per_N_gev_inverse']['real'][0][0],
                               8*C[0,0,0,0,0,0].real*q[2,2,0,0])

    def test_nonweinberg_output_never_claims_mass(self):
        C,q,mass=self.fixtures()
        with patch.dict(sys.modules,{'Numerical.joint_scalar_matching.core':self.core(.2)}):
            result=b1.evaluate_response(C,q,mass,1e5,1e3)
        self.assertFalse(result['full_SU2_Weinberg_projection_pass'])
        self.assertIsNone(result['neutrino_mass_singular_values_ev_per_abs_NGamma'])

    def test_scale_inversion_rejected(self):
        C,q,mass=self.fixtures()
        with patch.dict(sys.modules,{'Numerical.joint_scalar_matching.core':self.core()}):
            with self.assertRaises(b1.PilotError):b1.evaluate_response(C,q,mass,1e3,1e5)


if __name__=='__main__':unittest.main()
