"""Special scalar v6 tests: real rotations, masses, provenance and C5 path."""
from __future__ import annotations

from dataclasses import replace
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import sympy as sp

from Numerical.joint_scalar_matching.source_quartic_v4 import QuarticChannels
from Numerical.joint_scalar_matching.fermion_boundary_bridge import SourceBoundary
from Numerical.joint_scalar_matching.numerical_benchmark_v5 import RunningSnapshot,evaluate_ordinary_response
from Numerical.joint_scalar_matching.special_mass_v6 import (
    SPECIAL,MU_PARAMETER,DIMS,expected_special_labels,kind_for_model,
    build_special_basis,build_special_background,_mixed_rotation,_self_rotation,
    _native_self_mass_block)
from Numerical.joint_scalar_matching.large_batch import load_background,parse_matrix
from Numerical.experimental.T3SelfBilinearMassBasis import SelfBilinearMassBasis
from Numerical.experimental.GeneralT3MixedScalarMassBasis import MixedDoubletMassBasis
from test_t3_joint_scalar_integrated_batch import clean_synthetic_source
from scripts.RunT3JointScalarAll16V6 import numerical_metrics


def fake_special(model,root,mu=2300.0+120j):
    source,Q,_=clean_synthetic_source(model,root/model)
    state={'lambdaT3':{'real':.3,'imag':0},
           MU_PARAMETER[kind_for_model(model)]:{'real':mu.real,'imag':mu.imag}}
    snapshot=RunningSnapshot(model,root/'fake-running.json','c'*64,1e5,1e3,
                            1.1e6,1.9e6,state,{'illustrative':True})
    stub=replace(source,labels=expected_special_labels(model))
    basis=build_special_basis(model,snapshot,stub)
    R=basis['R_original_from_boundary']
    raw=np.einsum('fgijab,ak,bl->fgijkl',source.raw,R,R,optimize=True)
    by=np.einsum('tfgijab,ak,bl->tfgijkl',source.by_source,R,R,optimize=True)
    bound=replace(stub,raw=raw,by_source=by)
    # Reconstruct a source-derived symbolic channel with one test coefficient
    # that evaluates to the exact independent Ward-complete fixture Q.
    x=sp.Symbol('lambdaT3')
    entries={}
    for A in range(4):
        for B in range(A,4):
            for a in range(len(source.labels)):
                for b in range(a,len(source.labels)):
                    v=float(Q[A,B,a,b])
                    if v:entries[(A,B,a,b)]=sp.Float(v/.3)*x
    source_labels=('NewScalar1',)*(2*DIMS[model[0]][0])+('NewScalar2',)*(2*DIMS[model[0]][1])
    ch=QuarticChannels(model,root/'quartic-seed','f'*64,source_labels,
                       {'lambdaT3':entries},1,1,'test-fixed-real-couplings')
    return ch,bound,snapshot,basis


class SpecialV6Tests(unittest.TestCase):
    def test_eight_models_and_exact_boundary_real_field_layouts(self):
        self.assertEqual(len(SPECIAL),8)
        for model in sorted(SPECIAL):
            labels=expected_special_labels(model)
            self.assertEqual(len(labels),2*sum(DIMS[model[0]]))
            if model in ('B_-1','C_-1'):
                self.assertEqual(labels,('NewScalarX',)*4+('NewScalarY',)*4)
            else:
                self.assertTrue(any('Minus' in name for name in labels))
                self.assertTrue(any('Plus' in name for name in labels))

    def test_special_mass_exact_eigenspectra_and_orthogonality_all_eight(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            for model in sorted(SPECIAL):
                with self.subTest(model=model):
                    _,_,snap,info=fake_special(model,root)
                    R=info['R_original_from_boundary']
                    mass=info['mass2_boundary']
                    self.assertLess(np.max(abs(R.T@R-np.eye(R.shape[0]))),2e-13)
                    self.assertGreater(np.linalg.eigvalsh(mass).min(),0)
                    self.assertAlmostEqual(np.trace(mass),2*DIMS[model[0]][0]*snap.scalar_mass1_sq+
                        2*DIMS[model[0]][1]*snap.scalar_mass2_sq,places=5)

    def test_nonzero_mixed_eps_and_complex_mu_phase_is_preserved(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            _,_,snap,info=fake_special('B_-1',root,mu=2e4+8e3j)
            bas=MixedDoubletMassBasis.build('B_-1',snap.scalar_mass1_sq,snap.scalar_mass2_sq,complex(2e4,8e3))
            R=info['R_original_from_boundary']
            # Comparison to original-fields method from source -- not just eigenvalues.
            for i in (0,2,4,6):
                basis=np.zeros(8);basis[i]=1
                old=R@basis
                X=(basis[:4:2]+1j*basis[1:4:2])/np.sqrt(2)
                Y=(basis[4::2]+1j*basis[5::2])/np.sqrt(2)
                s1,s2=bas.original_fields(X,Y)
                predicted=np.column_stack((np.r_[s1,s2].real,np.r_[s1,s2].imag)).reshape(-1)*np.sqrt(2)
                np.testing.assert_allclose(old,predicted,atol=1e-13)
            self.assertGreater(abs(R[4,2])+abs(R[5,0]),0)

    def test_native_matchete_triplet_split_keeps_su2_indices_unmixed(self):
        # RunT3GeneralSelfMassBasis.wl: S[k]=(Minus[k]+I Plus[k])/Sqrt[2].
        # The spherical-prototype W is NOT a Matchete component map.
        for dim in (1,3):
            R=_self_rotation(dim)
            self.assertEqual(R.shape,(2*dim,2*dim))
            np.testing.assert_allclose(R.T@R,np.eye(2*dim),atol=0,rtol=0)
            for k in range(dim):
                v=np.zeros(2*dim);v[k]=1.0
                expect=np.zeros(2*dim);expect[2*k]=1.0
                np.testing.assert_array_equal(R@v,expect)
                v=np.zeros(2*dim);v[dim+k]=1.0
                expect=np.zeros(2*dim);expect[2*k+1]=1.0
                np.testing.assert_array_equal(R@v,expect)
            if dim==3:
                # Spherical W would mix k=0 and k=2 and complex phases.
                from Numerical.experimental.T3SelfBilinearMassBasis import real_intertwiner
                self.assertFalse(np.allclose(real_intertwiner(3),np.eye(3)))

    def test_native_mass_hessian_matches_potential_for_complex_mu_and_both_signs(self):
        # Matchete L_int=+mu(Sdag)^2+h.c. => V2 has the minus sign.
        rng=np.random.default_rng(190)
        for dim in (1,3):
            for mu in (2700.0, -2700.0, 900.0+2300.0j):
                with self.subTest(dim=dim,mu=mu):
                    m=1.7e6
                    M=_native_self_mass_block(dim,m,mu)
                    ev=np.sort(np.repeat([m-2*abs(mu),m+2*abs(mu)],dim))
                    np.testing.assert_allclose(np.linalg.eigvalsh(M),ev,rtol=2e-14)
                    for _ in range(4):
                        q=rng.normal(size=2*dim)
                        minus=q[:dim];plus=q[dim:]
                        S=(minus+1j*plus)/np.sqrt(2)
                        potential=(m*np.vdot(S,S).real-2*np.real(mu*np.sum(S.conj()**2)))
                        np.testing.assert_allclose(.5*float(q@M@q),potential,rtol=3e-15,atol=1e-8)
                    self.assertEqual(M[0,dim],-2*complex(mu).imag)
                    self.assertEqual(M[0,0],m-2*complex(mu).real)

    def test_negative_mu_fixed_boundary_orientation_and_sorted_eigenvalues(self):
        with tempfile.TemporaryDirectory() as td:
            ch,b,s,_=fake_special('E_-2',Path(td),mu=-3600.0+300j)
            info=build_special_basis('E_-2',s,b)
            R=info['R_original_from_boundary']
            self.assertEqual(R.shape,(12,12))
            np.testing.assert_allclose(R[6:12,6:12],_self_rotation(3))
            # Eigenvalues are correct although muS mass is not diagonal in
            # the muF Matchete named Minus/Plus basis.
            actual=np.linalg.eigvalsh(info['mass2_boundary'])
            expected=np.sort([s.scalar_mass1_sq]*6+
                [s.scalar_mass2_sq-2*abs(-3600+300j)]*3+
                [s.scalar_mass2_sq+2*abs(-3600+300j)]*3)
            np.testing.assert_allclose(actual,expected,atol=1e-8)
            self.assertGreater(abs(info['mass2_boundary'][6,9]),0)

    def test_self_mass_splittings_use_exact_two_times_mu(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            for model in ('A_-2','A_+0','D_-2','D_+0','E_-2','E_+0'):
                _,_,snap,info=fake_special(model,root,mu=2000+300j)
                m=snap.scalar_mass1_sq if kind_for_model(model)=='self_s1' else snap.scalar_mass2_sq
                self.assertAlmostEqual(info['mass_eigenvalues']['plus_gev2']-m,2*abs(2000+300j),places=8)
                self.assertAlmostEqual(m-info['mass_eigenvalues']['minus_gev2'],2*abs(2000+300j),places=8)

    def test_missing_or_tachyonic_mu_never_silently_set_to_zero(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            _,b,s,_=fake_special('E_-2',root)
            s_missing=replace(s,couplings={'lambdaT3':.3})
            with self.assertRaisesRegex(ValueError,'missing muT3GeneralS2Sq'):
                build_special_basis('E_-2',s_missing,b)
            state=dict(s.couplings)
            state['muT3GeneralS2Sq']={'real':2e6,'imag':0}
            with self.assertRaisesRegex(ValueError,'not positive'):
                build_special_basis('E_-2',replace(s,couplings=state),b)

    def test_boundary_labels_must_match_actual_special_release(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            _,b,s,_=fake_special('D_+0',root)
            with self.assertRaisesRegex(ValueError,'ordering differs'):
                build_special_basis('D_+0',s,replace(b,labels=tuple(reversed(b.labels))))

    def test_full_background_and_nonzero_su2_complete_matching_self_and_mixed(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            for model in ('A_-2','B_-1','E_+0'):
                with self.subTest(model=model):
                    ch,b,s,_=fake_special(model,root,mu=0j)
                    s=replace(s,scalar_mass1_sq=1.5e6,scalar_mass2_sq=1.5e6)
                    info=build_special_background(ch,b,s,root/model/'results')
                    background=load_background(Path(info['background_json']),b)
                    self.assertTrue(np.isfinite(background.quartic).all())
                    result=evaluate_ordinary_response(b,Path(info['background_json']))
                    mat=parse_matrix(result['response']['C5_LLSS_total_one_loop_per_N_gev_inverse'])
                    self.assertGreater(float(np.linalg.norm(mat)),1e-19)
                    self.assertLess(result['response']['weinberg_beta_projection_residual'],1e-8)
                    self.assertLess(result['response']['relative_matching_scale_cancellation_residual'],2e-9)
                    self.assertFalse(result['physical_C5_certified'])
                    metrics=numerical_metrics(result)
                    self.assertEqual(len(metrics['conditional_neutrino_mass_singular_values_eV_per_abs_N']),3)
                    self.assertTrue(metrics['conditional_neutrino_mass_singular_values_eV_per_abs_N'][0]>=0)


if __name__=='__main__': unittest.main()
