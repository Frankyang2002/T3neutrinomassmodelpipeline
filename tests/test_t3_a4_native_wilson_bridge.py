"""Offline CG/basis tests (no Matchete/Mathematica and no fake physical C5)."""
from __future__ import annotations

import copy
import math
import sys
import types
import unittest
from unittest.mock import patch
from tempfile import TemporaryDirectory
from pathlib import Path

import numpy as np

# Real project imports its authoritative Aminus4EFT1Couplings. This minimal
# fallback exists only for running these tests in a source-file-only environment.
try:
    from Numerical.running.Aminus4EFT1Couplings import su2_generators
except ModuleNotFoundError:
    fake_su2=types.ModuleType('Numerical.running.Aminus4EFT1Couplings')
    def su2_generators(dim):
        spin=(dim-1)/2; ms=[spin-i for i in range(dim)]
        plus=np.zeros((dim,dim),complex)
        for col,m in enumerate(ms):
            if m+1 in ms:
                plus[ms.index(m+1),col]=math.sqrt((spin-m)*(spin+m+1))
        minus=plus.conj().T
        return ((plus+minus)/2,(plus-minus)/(2j),np.diag(ms).astype(complex))
    fake_su2.su2_generators=su2_generators
    sys.modules['Numerical.running.Aminus4EFT1Couplings']=fake_su2

from Numerical.running.T3A4NativeWilsonBridge import (
    NativeWilsonBridgeError, _SOURCE_CG, derive_triplet_map,
    formal_tree_master, source_cg_from_seed, source_validation,
    spherical_cg, verify_su2_singlet_cg, three_flavour_wbar,
)
from Numerical.parallel_t3.background.general_batch1_api import encode_native


def seed() -> dict:
    terms = [
        {'Chirality':'PR','TermInputForm':
         'Field[NewScalar1, Scalar, {}, {}]*Bar[Field[NewScalar2, Scalar, {}, {}]]*'
         'CG[Bar[T3MixCG],{}]*Coupling[y1,{},0]*Coupling[y2,{},0]*'
         'Coupling[MF,{},0]*GammaCC'},
        {'Chirality':'PL','TermInputForm':
         'Bar[Field[NewScalar1, Scalar, {}, {}]]*Field[NewScalar2, Scalar, {}, {}]*'
         'CG[T3MixCG,{}]*Coupling[y1,{},0]*Coupling[y2,{},0]*'
         'Coupling[MF,{},0]*GammaCC'}
    ]
    return {'IntegratedFields':['F'],'ActiveHeavyFields':['S1','S2'],
            'TreeWilsonTermCount':2,'TreeWilsonTerms':terms,
            'CGRegistry':[{'Name':'T3MixCG',
              'RepsInputForm':'{Bar[SU2L[fund]], Bar[SU2L[fund]], T3BSMd3}',
              'TensorInputForm':_SOURCE_CG}]}


def first_boundary(M=None) -> dict:
    if M is None:
        M = np.diag([100000.,110000.,120000.]).astype(complex)
    y1 = np.array([[1+1j,2,3],[0.2,-0.3j,0.4],[1,1j,2]],complex)*1e-4
    y2 = np.array([[2j,-1,3],[0.5,2,1j],[1,-0.3,4]],complex)*1e-4
    return {'model':'A_-4','all_heavy_fermions_to_integrate_out_together':True,
            'scale_GeV':1e5,
            'native_uv_rgbeta_couplings_at_muF':{
                'MF':encode_native(M),'y1':encode_native(y1),
                'y2':encode_native(y2)}}


class NativeCGBridgeTests(unittest.TestCase):
    def test_exact_saved_cartesian_tensor(self):
        c=source_cg_from_seed(seed())
        self.assertEqual(c.shape,(2,2,3))
        self.assertEqual(np.count_nonzero(c),6)
        self.assertAlmostEqual(c[0,0,0].real,3**(-.25),14)
        self.assertAlmostEqual(c[1,1,2].imag,3**(-.25),14)
        self.assertTrue(np.allclose(c,c.swapaxes(0,1)))

    def test_triplet_intertwiner_derived_unitary(self):
        a=derive_triplet_map(seed())
        U=a.cartesian_from_spherical
        target=np.array([[-1/math.sqrt(2),0,1/math.sqrt(2)],
                         [0,1,0],
                         [1j/math.sqrt(2),0,1j/math.sqrt(2)]],complex)
        np.testing.assert_allclose(U,target,atol=5e-15)
        np.testing.assert_allclose(U.conj().T@U,np.eye(3),atol=5e-15)
        self.assertAlmostEqual(a.spherical_cg_scale,math.sqrt(2)/3**.25,14)
        self.assertLess(a.max_reconstruction_residual,1e-14)

    def test_exact_source_s2_contraction_in_spherical_coordinates(self):
        c=source_cg_from_seed(seed())
        m=derive_triplet_map(seed())
        z=np.array([2+3j, -1+4j, 0.4-0.2j])
        # C_mat[ij,k]*S2_mat[k] = scale*C_sph[ij,s]*S2_sph[s].
        lhs=np.einsum('ijk,k->ij',c,m.cartesian_from_spherical.conj()@z)
        rhs=m.spherical_cg_scale*np.einsum('ijs,s->ij',spherical_cg(),z)
        np.testing.assert_allclose(lhs,rhs,atol=1e-14)

    def test_spherical_Wilson_cg_is_su2_invariant(self):
        self.assertLess(verify_su2_singlet_cg(spherical_cg()),1e-14)

    def test_real_scalar_map_orthogonal(self):
        m=derive_triplet_map(seed())
        R=m.real_s2_cartesian_from_spherical
        np.testing.assert_allclose(R@R.T,np.eye(6),atol=5e-15)
        z=np.array([0.5+2j,-3+0.1j,1-4j]); x=np.array([u for v in z for u in (v.real,v.imag)])
        actual=R@x
        expected=m.cartesian_from_spherical.conj()@z
        np.testing.assert_allclose(actual,np.array([u for v in expected for u in (v.real,v.imag)]))

    def test_reject_changed_source_sparse_tensor(self):
        bad=seed()
        bad['CGRegistry'][0]['TensorInputForm']=_SOURCE_CG.replace('I/3^(1/4)', '-I/3^(1/4)',1)
        with self.assertRaises(NativeWilsonBridgeError):
            derive_triplet_map(bad)

    def test_reject_wrong_CG_representation(self):
        bad=seed()
        bad['CGRegistry'][0]['RepsInputForm']='{SU2L[fund],SU2L[fund],T3BSMd3}'
        with self.assertRaises(NativeWilsonBridgeError):
            derive_triplet_map(bad)

    def test_reject_missing_conjugate_Wilson(self):
        bad=seed();bad['TreeWilsonTerms'][1]['Chirality']='PR'
        with self.assertRaises(NativeWilsonBridgeError):
            derive_triplet_map(bad)

    def test_real_mass_three_flavour_Wbar(self):
        f=first_boundary()
        W,report=three_flavour_wbar(f)
        from Numerical.parallel_t3.background.general_batch1_api import decode_native
        d=f['native_uv_rgbeta_couplings_at_muF']
        y1=np.asarray(decode_native(d['y1']),complex)
        y2=np.asarray(decode_native(d['y2']),complex)
        expected=(y1.conj()/np.array([1e5,1.1e5,1.2e5]))@y2.conj().T
        np.testing.assert_allclose(W,expected,atol=1e-28)
        self.assertEqual(report['heavy_fermions'],3)

    def test_reject_nondiagonal_MF_without_rotation(self):
        M=np.diag([1e5,1.1e5,1.2e5]).astype(complex);M[0,1]=.2
        with self.assertRaises(NativeWilsonBridgeError):
            three_flavour_wbar(first_boundary(M))

    def test_reject_old_or_exploratory_first_threshold(self):
        tree='GammaCC GammaCC Coupling[T3EFT1W12, Coupling[T3EFT1W12,'
        status={'model':'A_-4','first_threshold_reused':False,
                'A_joint_scalar_match_performed':True,
                'A_physical_scalar_soft_Weinberg_terms':10,
                'exploratory_source_theory_used':False}
        lag={'CheckLagrangian_original':'True','accepted_canonicalisation':'ExactSavedTree'}
        source_validation(seed(),tree,status,lag)
        status['exploratory_source_theory_used']=True
        with self.assertRaises(NativeWilsonBridgeError):
            source_validation(seed(),tree,status,lag)

    def test_existing_B1_RG_adapter_without_reintegration(self):
        # Interface-level test: use solved, labelled Batch1 arrays instead of
        # starting a new renormalizable RGBeta calculation.
        from scripts.RunT3A4WilsonMixing import build_running_context
        from Numerical.running.Aminus4EFT1Couplings import su2_generators as su2
        from types import SimpleNamespace
        names=('g2','g3','gY','lambda12','lambdaH','lambdaH1','lambdaH2',
               'lambdaH2Adj','lambdaS1','lambdaS2','lambdaS2Adj','lambdaT3',
               'mS1Sq','mS2Sq','yd','ye','yu')
        manifest={'model':'A_-4','coupling_order':{'eft1':list(names)}}
        fake_native=types.ModuleType('Numerical.parallel_t3.background.general_batch1_api')
        fake_native.NativeBatch1Trajectory=SimpleNamespace(load=lambda _:SimpleNamespace(manifest=manifest))
        basis=types.ModuleType('Numerical.running.Aminus4EFT1Couplings')
        basis.BETA_NAMES=frozenset(names)
        basis.F_LABELS=tuple('f'+str(i) for i in range(45))
        basis.S_LABELS=tuple('s'+str(i) for i in range(12))
        basis.su2_generators=su2
        def point(mu,values):
            return SimpleNamespace(Y=np.zeros((45,45,12)),
              lambda4=np.zeros((12,12,12,12)),
              couplings={k:float(values[k].real) for k in ('g2','gY','g3')})
        basis.tensor_point=point
        def generators():
            return {k:{'fermion':(np.zeros((45,45)),),
                       'scalar':(np.zeros((12,12)),)} for k in ('SU2','U1','SU3')}
        basis.gauge_generators=generators
        parent=types.ModuleType('Numerical.running.T3EFT1WilsonODETransport')
        parent.Aminus4Trajectory=lambda *args:SimpleNamespace(mu=args[1],
              t=args[0],y=args[2],quartic=args[3],
              fermion_labels=args[4],scalar_labels=args[5],
              couplings=args[6],gauge=args[7])
        with TemporaryDirectory() as d:
            folder=Path(d)
            mu=np.array([1e5,1e4,1e3]);arrays={'eft1__mu_gev':mu}
            for key in names:
                if key in ('yu','ye','yd'):
                    arrays['eft1__'+key]=np.zeros((3,3,3),complex)
                else:
                    arrays['eft1__'+key]=np.array([.5,.6,.7],complex)
            np.savez(folder/'trajectories.npz',**arrays)
            with patch.dict(sys.modules,{
                'Numerical.parallel_t3.background.general_batch1_api':fake_native,
                'Numerical.running.Aminus4EFT1Couplings':basis,
                'Numerical.running.T3EFT1WilsonODETransport':parent}):
                tr=build_running_context(folder,1e5,1e3)
            np.testing.assert_array_equal(tr.mu,mu[::-1])
            self.assertEqual(tr.y.shape,(3,45,45,12))
            self.assertEqual(tr.quartic.shape,(3,12,12,12,12))
            self.assertEqual(len(tr.couplings),3)

    def test_full_tree_master_orbits_and_polynomial(self):
        names=tuple(f'L[{p},{i}]' for p in (1,2,3) for i in (1,2))
        fake=types.ModuleType('Numerical.running.Aminus4EFT1Couplings')
        fake.F_LABELS = names + tuple(f'x{k}' for k in range(39))
        fake.FI = {x:i for i,x in enumerate(fake.F_LABELS)}
        fake.su2_generators = su2_generators
        fake.S_LABELS=tuple(f'H[{n},{ri}]' for n in (1,2) for ri in ('R','I')) + (
            'S1[1,R]','S1[1,I]') + tuple(f'S2[{m},{ri}]'
                for m in (1,2,3) for ri in ('R','I'))
        w,_=three_flavour_wbar(first_boundary())
        with patch.dict(sys.modules,{'Numerical.running.Aminus4EFT1Couplings':fake}):
            C=formal_tree_master(w,derive_triplet_map(seed()))
        self.assertEqual(C.shape,(45,45,12,12))
        self.assertTrue(np.allclose(C,C.swapaxes(0,1),atol=1e-26))
        self.assertTrue(np.allclose(C,C.swapaxes(2,3),atol=1e-26))
        self.assertTrue(np.allclose(C[:6,:6,:4,:4],0))
        # The two symmetrizations reproduce the original source polynomial.
        psi=np.array([1+.2j,-2+.4j,.5+2j,3-1j,.4j,2-.7j] + [0]*39)
        phi=np.array([0,0,0,0, 0.5,-0.2, 0.4,0.6, -0.3,0.7, 0.9,-0.8])
        pol=np.einsum('ijab,i,j,a,b->',C,psi,psi,phi,phi,optimize=True)
        zS1=(phi[4]-1j*phi[5])/math.sqrt(2)
        zs2=(phi[np.array([6,8,10])]+1j*phi[np.array([7,9,11])])/math.sqrt(2)
        cg=spherical_cg()
        expected=0j
        for p in range(3):
            for q in range(3):
                for i in range(2):
                    for j in range(2):
                        for k in range(3):
                            expected+=(w[p,q]*derive_triplet_map(seed()).spherical_cg_scale *
                                cg[i,j,k]*psi[2*p+i]*psi[2*q+j]*zS1*zs2[k])
        self.assertAlmostEqual(abs(pol-expected),0,places=23)


if __name__=='__main__':
    unittest.main()
