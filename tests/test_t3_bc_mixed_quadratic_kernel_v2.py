"""Standalone and repository-aware regression tests, no source modifications."""
from __future__ import annotations

import json
import math
from pathlib import Path
from unittest import mock
from types import SimpleNamespace
import tempfile
import unittest

import numpy as np

from Numerical.experimental.T3BCMixedQuadraticKernel import (
    positive_mass_matrix, signed_exact_master, strict_heavy_mass_master,
    hard_soft_master, spectral_kernel, euclidean_matrix_kernel,
    invariant_quartic_contraction,
)
from scripts.RunT3BCMixedQuadraticKernelAudit import (
    read, complex_array, load_result, strict_full_mass_response,
)

ROOT=Path(__file__).resolve().parents[1]


class MasterTests(unittest.TestCase):
    def test_equal_mass_continuous_limit(self):
        self.assertAlmostEqual(signed_exact_master(49.,49.,49.),-0.5,places=13)
        self.assertAlmostEqual(signed_exact_master(1e8,1e6,1e6),
                               signed_exact_master(1e8,1e6+0.01,1e6-0.01),places=7)

    def test_mass_function_symmetry(self):
        a=signed_exact_master(1e10,1e6,1.3e6)
        b=signed_exact_master(1e10,1.3e6,1e6)
        self.assertAlmostEqual(a,b,places=12)
        self.assertLess(a,0)

    def test_hard_soft_three_scales(self):
        mf=1e10; ma=1e6; mb=1.4e6
        strict=strict_heavy_mass_master(mf,ma,mb)
        for q in (1e6,1.4e6,1e10):
            res=hard_soft_master(mf,ma,mb,q)
            self.assertLess(abs(res['hard']+res['soft']-strict),1e-12)
        self.assertLess(abs(signed_exact_master(mf,ma,mb)-strict),0.01)

    def test_hard_soft_equal_scalar_masses(self):
        mf=3e10; m=1e6
        d=hard_soft_master(mf,m,m,1.4e6)
        self.assertAlmostEqual(d['sum'],1+math.log(m/mf),places=12)

    def test_unphysical_mass_rejected(self):
        with self.assertRaises(ValueError):signed_exact_master(-1,1,1)
        with self.assertRaises(ValueError):positive_mass_matrix(np.array([[1,2],[2,1]]))
        with self.assertRaises(ValueError):positive_mass_matrix(np.array([[1,1j],[1j,2]]))


class ResolventTests(unittest.TestCase):
    def setUp(self):
        self.A=np.array([[1.3e6,2.2e4+1.1e4j],[2.2e4-1.1e4j,1e6]],complex)

    def test_integrate_exact_matrix(self):
        out=euclidean_matrix_kernel(self.A,1e10)
        self.assertLess(out['original_basis_vs_spectral_relative'],3e-9)
        self.assertLess(out['dyson_relative_errors']['8'],3e-9)
        self.assertGreater(out['dyson_relative_errors']['0'],1e-4)
        self.assertGreater(out['dyson_relative_errors']['1'],out['dyson_relative_errors']['2'])

    def test_spectral_kernel_and_vertex_basis(self):
        K=spectral_kernel(self.A,1e10)
        self.assertEqual(K.shape,(4,4))
        self.assertLess(np.linalg.norm(K-K.conj().T),1e-12)
        self.assertLess(invariant_quartic_contraction(self.A,1e10),1e-12)

    def test_phase_reversal_conjugates_full_kernel(self):
        x=spectral_kernel(self.A,1e10)
        y=spectral_kernel(self.A.conj(),1e10)
        self.assertLess(np.linalg.norm(y-x.conj()),1e-12)

    def test_diagonal_zero_mix_is_dyson_exact(self):
        A=np.diag([1e6,1.4e6]).astype(complex)
        out=euclidean_matrix_kernel(A,2e10)
        self.assertLess(out['dyson_relative_errors']['0'],3e-8)

    def test_matrix_labels_not_independent_of_components(self):
        K=spectral_kernel(self.A,1e10)
        self.assertGreater(abs(K[0,1]),0)
        self.assertTrue(np.isfinite(K).all())


class SourceContracts(unittest.TestCase):
    def _sample(self):
        A=np.array([[1.355e6,1.3e4+6e3j],[1.3e4-6e3j,1.03e6]])
        eig,U=np.linalg.eigh(A)
        return {'format':'t3-bc-common-uv-scale-restricted-direct-matchete-triangle-v1',
                'status':'RunningCommonScaleRestrictedTriangleEvaluatedNotPhysicalC5',
                'model':'B_-1',
                'common_running_parameters_verified_same_UV_scale':True,
                'physical_C5_certified':False,
                'full_fermion_X_Y_sequential_RG_transport':False,
                'general_extra_mixed_mass_source_matching_certified':False,
                'native_matchete_rgbeta_action_dictionary_certified':False,
                'running_input_direct_vs_sequential_parity':{
                    'status':'RunningInputRestrictedFullMassTriangleParityPassed',
                    'max_relative_residual':1e-15},
                'scalar_mass_basis':{
                    'mass_squared_flavour_matrix_gev2':{
                        'real':A.real.tolist(),'imag':A.imag.tolist()},
                    'unitary_flavour_rotation_U':{
                        'real':U.real.tolist(),'imag':U.imag.tolist()}},
                'muX_gev':float(np.sqrt(eig[0])),
                'muY_gev':float(np.sqrt(eig[1])),
                'MF_takagi_masses_gev':[1e5,1.1e5,1.2e5]}

    def test_valid_source_data(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t)/'B_m1_running_one_scale.json'
            p.write_text(json.dumps(self._sample()),encoding='utf-8')
            file,data,A,U,mf=load_result('B_-1',Path(t))
            self.assertEqual(file,p)
            self.assertEqual(tuple(mf),(1e5,1.1e5,1.2e5))

    def test_cannot_promote_uncertified_source(self):
        s=self._sample();s['physical_C5_certified']=True
        with tempfile.TemporaryDirectory() as t:
            (Path(t)/'B_m1_running_one_scale.json').write_text(json.dumps(s))
            with self.assertRaises(ValueError):load_result('B_-1',Path(t))

    def test_cannot_accept_broken_matching_parity(self):
        s=self._sample();s['running_input_direct_vs_sequential_parity']['max_relative_residual']=0.2
        with tempfile.TemporaryDirectory() as t:
            (Path(t)/'B_m1_running_one_scale.json').write_text(json.dumps(s))
            with self.assertRaises(ValueError):load_result('B_-1',Path(t))

    def test_complex_schema(self):
        with self.assertRaises(ValueError):complex_array({'real':[[0]],'imag':[[0]]},(2,2))

    def test_restricted_matchete_full_vs_strict_fake_source_contract(self):
        # Test the integration layer with separately constructed matched source
        # polynomials; no claim this replaces the live source-dependent check.
        b=self._sample(); A=complex_array(b['scalar_mass_basis']['mass_squared_flavour_matrix_gev2'],(2,2))
        _,U=np.linalg.eigh(A)
        shape=np.array([[.34,.08,.13],[.05,-.27,.12],[.03,.17,-.21]],complex)*1e-3
        y1=shape; y2=shape.T*.7
        b['running_heavy_eigenbasis_y1_3x3']={'real':y1.real.tolist(),'imag':y1.imag.tolist()}
        b['running_heavy_eigenbasis_y2_3x3']={'real':y2.real.tolist(),'imag':y2.imag.tolist()}
        b['conditional_matchete_quartics']={k:{'real':v,'imag':0.}
            for k,v in zip(('lambdaT3','lambdaT3ExtraQ020020I1','lambdaT3ExtraQ020200I1'),(.1,.01,-.02))}
        b['reference_EW_vev_gev']=246.22
        raw=((.6,.02,.5),(.1,.4,.03),(.02,.05,.2))
        names=tuple(b['conditional_matchete_quartics'])
        def coeff(weights,mf,mx,my,strict):
            phi=[(strict_heavy_mass_master if strict else signed_exact_master)(mf**2,mx**2,mx**2),
                 (strict_heavy_mass_master if strict else signed_exact_master)(mf**2,my**2,my**2),
                 (strict_heavy_mass_master if strict else signed_exact_master)(mf**2,mx**2,my**2)]
            return tuple(sum(weights[k][j]*phi[k] for k in range(3))/(16*math.pi**2*mf)
                         for j in range(3))
        fake22=SimpleNamespace(Q_NAMES=names, sequential_weights=lambda *args:raw,
              substitution=lambda *args:{},
              expected_ordered=lambda w,mf,mx,my,v17:coeff(w,mf,mx,my,False))
        fake18=SimpleNamespace(polynomial_values=lambda v,sub: np.asarray(v,complex))
        from scripts import RunT3BCMixedQuadraticKernelAudit as audit
        weights=tuple(np.asarray(v,complex) for v in raw)
        mf=np.asarray(b['MF_takagi_masses_gev'])
        exact=audit._assemble(mf,y1,y2,lambda m:coeff(weights,m,b['muX_gev'],b['muY_gev'],False))
        b['candidate_ordered_source_C5_3x3_gev_inverse_NOT_PHYSICAL']={
          'real':exact.real.tolist(),'imag':exact.imag.tolist()}
        with mock.patch.object(audit,'source_module',side_effect=[fake22,fake18,object(),object()]):
            r=strict_full_mass_response('B_-1',b,U,mf)
        self.assertLess(r['source_matchete_full_mass_parity_relative'],1e-12)
        self.assertGreater(r['leading_1_over_MF_restricted_C5_relative_error'],0)
        self.assertFalse(r['full_renormalized_EFT_counterterms_computed'])

    def test_existing_real_source_when_present(self):
        folder=ROOT/'output/studies/t3_bc_running_one_scale_triangle_v1'
        p=folder/'B_m1_running_one_scale.json'
        if not p.is_file():self.skipTest('needs full Windows project saved B running result')
        _,d,A,U,m=load_result('B_-1',folder)
        self.assertLess(euclidean_matrix_kernel(A,m[0]**2)['original_basis_vs_spectral_relative'],3e-8)


if __name__=='__main__':unittest.main()
