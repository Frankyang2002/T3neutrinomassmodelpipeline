"""Source-independent math tests, plus native source contract tests when installed."""
from __future__ import annotations
import importlib
import json
import math
from pathlib import Path
import sys
import types
import unittest
import numpy as np

# Allow package standalone test outside the complete user's repository.
ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
try:
    import Numerical.experimental.GeneralT3MixedScalarMassBasis
except ModuleNotFoundError:
    # Offline ONLY: no claim of live/native Matchete scalar-basis verification.
    m=types.ModuleType('Numerical.experimental.GeneralT3MixedScalarMassBasis')
    class MockBasis:
        def __init__(self,model,m1,m2,mu):
            self.model=model;self.mu_sq=mu;self.m_s1_sq=m1;self.m_s2_sq=m2
            self.eigenvalues_sq,self.U=np.linalg.eigh(np.array([[m1,mu],[mu.conjugate(),m2]],complex))
            self.tachyonic=bool(self.eigenvalues_sq[0]<0);self.degenerate=False
        @classmethod
        def build(cls,*args):return cls(*args)
        def verify(self):
            np.testing.assert_allclose(self.U.conj().T@self.U,np.eye(2),atol=1e-12)
        def info(self):return {'masses_squared_gev2':self.eigenvalues_sq.tolist()}
    m.MixedDoubletMassBasis=MockBasis
    m.quadratic_mass_tensor_transport=lambda b: {'mass_tensor_diagonalization_residual_gev2':0.}
    m.verify_tensor_transport=lambda *args,**kwargs:0.
    sys.modules['Numerical.experimental.GeneralT3MixedScalarMassBasis']=m

from Numerical.experimental.T3BCRunningOneScaleTriangle import (
   scalar,matrix,packed,takagi_majorana,assemble_ordered_c5,
   uv_boundary,evaluate_one_scale,MATCH_QUARTICS
)


def cmat(A):
    A=np.asarray(A,complex)
    return {'real':A.real.tolist(),'imag':A.imag.tolist()}


def config(model='B_-1',mu=complex(9500,2200)):
    mf=np.diag([100000.,110000.,120000.]).astype(complex)
    y1=np.array([[1e-4,2e-4,3e-4],[-2e-4,3e-4,1e-4],[1e-4,-1e-4,2e-4]],complex)
    y2=1j*np.array([[1e-4,1e-4,2e-4],[2e-4,1e-4,-1e-4],[0,3e-4,1e-4]],complex)
    state={'representation':{'d_f':1 if model=='B_-1' else 3,'alpha':-1},
           'MF':cmat(mf),'y1':cmat(y1),'y2':cmat(y2),
           'mS1Sq':1.35e6,'mS2Sq':1.02e6,
           'lambdaT3':{'real':.1,'imag':0},
           'lambdaHHdagS2S2':{'real':-.01,'imag':0.002},
           'lambdaHHdagS1barS1bar':{'real':-.02,'imag':-.003}}
    return {'model':model,'status':'CandidateCoupledUVAndEFT1FiveRealMassesIntegrated',
            'solver_success':True,'physical_matchete_rgbeta_mapping_certified':False,
            'physical_neutrino_mass_prediction':False,
            'uv_range_gev':[1e7,1e5],
            'uv_threshold':{'mu_gev':1e5,'state':state,
                'muT3GeneralS1S2_gev2':{'real':mu.real,'imag':mu.imag}}}


class ShapeContractTests(unittest.TestCase):
    def test_complex_packing(self):
        x=np.diag([1+1j,2-2j,3j]);np.testing.assert_array_equal(matrix(cmat(x),'fixture'),x)
        self.assertEqual(packed(x),cmat(x))
        self.assertEqual(scalar({'real':3,'imag':-2},'test'),3-2j)
    def test_bad_matrix_rejected(self):
        with self.assertRaises(ValueError): matrix({'real':[[1]],'imag':[[0]]},'bad')
        with self.assertRaises(ValueError): matrix({'real':[[float('nan')]*3]*3,'imag':[[0]*3]*3},'bad')


class MajoranaTests(unittest.TestCase):
    def test_takagi_diagonal_and_complex_basis(self):
        rng=np.random.default_rng(42)
        a=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
        q,_=np.linalg.qr(a)
        masses=np.array([100000.,110000.,120000.])
        m=q.conj()@np.diag(masses)@q.conj().T
        got,U,metrics=takagi_majorana(m)
        np.testing.assert_allclose(got,masses,rtol=1e-9)
        np.testing.assert_allclose(U.T@m@U,np.diag(masses),rtol=2e-10,atol=1e-5)
        self.assertLess(metrics['takagi_residual'],1e-10)
    def test_majorana_not_symmetric_fails(self):
        m=np.eye(3,dtype=complex);m[0,1]=1
        with self.assertRaises(ValueError):takagi_majorana(m)
    def test_degenerate_heavy_basis_fails(self):
        with self.assertRaises(ValueError):takagi_majorana(np.diag([1e5,1e5,1.2e5]))


class CoefficientTests(unittest.TestCase):
    @staticmethod
    def simple_coefficient(mf,mx,my,U,Q):
        return 0.2/mf,0.4/mf,0.8/mf
    def test_full_three_heavy_and_complex_flavour_symmetry(self):
        ms=np.array([1e5,1.1e5,1.2e5])
        y1=np.diag([1,2+1j,3]);y2=np.diag([1j,2,1-2j])
        K=np.array([[.4,.4],[.4,1.6]],complex)
        want=sum(np.column_stack((y1[:,i].conj(),y2[:,i].conj()))@(K/ms[i])@
                 np.column_stack((y1[:,i].conj(),y2[:,i].conj())).T for i in range(3))
        got,pieces=assemble_ordered_c5(ms,y1,y2,1000,1200,np.eye(2),
                   {k:0j for k in MATCH_QUARTICS},self.simple_coefficient)
        np.testing.assert_allclose(got,want,rtol=1e-12,atol=1e-17)
        self.assertEqual(len(pieces),3)
        np.testing.assert_allclose(got,got.T)
    def test_degenerate_scalar_source_order_fails(self):
        y=np.eye(3,dtype=complex)
        with self.assertRaises(ValueError):
            assemble_ordered_c5(np.array([1000.,1100.,1200.]),y,y,1000,1200,np.eye(2),
                                {k:0j for k in MATCH_QUARTICS},self.simple_coefficient)
    def test_quartic_coverage_fails(self):
        with self.assertRaises(ValueError):
            assemble_ordered_c5(np.array([100000.,110000.,120000.]),np.eye(3),np.eye(3),
                                1000,1200,np.eye(2),{'lambdaT3':.1},self.simple_coefficient)
    def test_boundary_same_uv_scale_and_sign_mapping(self):
        result=uv_boundary('B_-1',config())
        self.assertAlmostEqual(result['mu_gev'],1e5)
        self.assertAlmostEqual(result['quartics']['lambdaT3'],.1)
        self.assertAlmostEqual(result['quartics']['lambdaT3ExtraQ020020I1'],.01-.002j)
        self.assertAlmostEqual(result['quartics']['lambdaT3ExtraQ020200I1'],.02+.003j)
        self.assertFalse(result['physical_coupling_mapping_certified'])
        self.assertFalse(result['physical_heavy_spinor_phase_mapping_certified'])
    def test_running_candidate_not_physical_c5(self):
        b=uv_boundary('C_-1',config('C_-1'))
        result=evaluate_one_scale(b,self.simple_coefficient,vev_gev=246.)
        self.assertFalse(result['physical_C5_certified'])
        self.assertFalse(result['heavy_spinor_phase_mapping_physically_certified'])
        self.assertIn('alternative_same_V_spinor_rotation_relative_C5_difference',result)
        self.assertEqual(np.asarray(result['running_heavy_eigenbasis_y1_3x3']['real']).shape,(3,3))
        self.assertFalse(result['general_extra_mixed_mass_source_matching_certified'])
        self.assertEqual(len(result['reference_mnu_singular_ev_descending_NOT_PHYSICAL']),3)
        self.assertTrue(all(x>=0 for x in result['reference_mnu_singular_ev_descending_NOT_PHYSICAL']))
    def test_scale_mismatch_fails(self):
        c=config();c['uv_threshold']['mu_gev']=1000
        with self.assertRaises(ValueError):uv_boundary('B_-1',c)
    def test_certification_flags_not_promoted(self):
        c=config();c['physical_matchete_rgbeta_mapping_certified']=True
        with self.assertRaises(ValueError):uv_boundary('B_-1',c)


class NativePresenceTests(unittest.TestCase):
    @unittest.skipUnless((ROOT/'scripts/DirectVsSequentialOneLoop.experimental.py').exists(),
                         'Requires user local complete repository')
    def test_real_direct_matchete_source_files_present(self):
        from scripts.RunT3BCRunningOneScaleTriangle import (
            import_source,archived_source,DIRECT,V18,V22,read)
        d=import_source(DIRECT,'t3_native_v22_test')
        a=import_source(V18,'t3_native_v18_test')
        source=read(V22)
        for model in ('B_-1','C_-1'):
            _,audit=archived_source(model,d,a,source)
            self.assertLess(audit['archived_frozen_C5_parity_relative'],2e-8)

if __name__=='__main__':unittest.main()
