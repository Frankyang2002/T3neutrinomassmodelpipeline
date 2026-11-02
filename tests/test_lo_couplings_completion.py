"""Unit regression for the all-16 HHSS baseline completion and exact extra spans."""
from __future__ import annotations

from types import SimpleNamespace
import unittest
import numpy as np
import sympy as sp

from Numerical.lo_certification.couplings.portal_normalization import rgbeta_portal_columns
from Numerical.lo_certification.couplings.portal_completion import (
    DIMS, exact_portal_completion, completed_portal_numeric)
from Numerical.lo_certification.couplings.extra_hhss_geometry import (
    reference_invariants, real_component_derivatives, model_extra_geometry)
from Numerical.lo_certification.couplings.quartic_completion_batch import run_completion_batch

MODELS=('A_-4','A_-2','A_+0','B_-3','B_-1','B_+1',
        'C_-3','C_-1','C_+1','D_-2','D_+0','D_+2',
        'E_-4','E_-2','E_+0','E_+2')

# Independent static model/monomial inventory copied from the checked
# GeneralScalarVerifiedQuartics field-content counts (HHSS entries ONLY).
SECTORS={
  'A_-4':[], 'A_-2':[(0,2,0,1,0,1),(1,1,0,0,0,2)],
  'A_+0':[(0,2,1,0,1,0),(1,1,0,2,0,0)],
  'B_-3':[(0,2,0,0,0,2)],
  'B_-1':[(0,2,0,0,2,0),(0,2,0,2,0,0),(1,1,0,1,0,1),(1,1,0,1,0,1)],
  'B_+1':[(0,2,2,0,0,0)],
  'C_-3':[(0,2,0,0,0,2)],
  'C_-1':[(0,2,0,0,2,0),(0,2,0,2,0,0),(1,1,0,1,0,1),(1,1,0,1,0,1)],
  'C_+1':[(0,2,2,0,0,0)],
  'D_-2':[(0,2,0,1,0,1),(1,1,0,0,0,2)],
  'D_+0':[(0,2,1,0,1,0),(1,1,0,2,0,0)],
  'D_+2':[],
  'E_-4':[], 'E_-2':[(0,2,0,1,0,1),(1,1,0,0,0,2)],
  'E_+0':[(0,2,1,0,1,0),(1,1,0,2,0,0)],'E_+2':[]}


def pretend_native(model):
    d1,d2=DIMS[model[0]]
    channels={}
    for species,d,offset in ((1,d1,0),(2,d2,2*d1)):
        rgb=rgbeta_portal_columns(d)
        if d==1:
            combinations=[{'singlet':sp.S(1)}]
        elif d==2:
            combinations=[{'singlet':sp.S(1),'adjoint':sp.S(1)},
                          {'singlet':sp.S(0),'adjoint':sp.S(2)}]
        else:
            combinations=[{'singlet':sp.S(1),'adjoint':sp.S(2)}]
        for count,combo in enumerate(combinations,1):
            name=f'lambdaH{species}Inv{count}'
            out={}
            for kind,val in combo.items():
                for key,expr in rgb[kind].items():
                    i,j,a,b=key
                    k=(i,j,offset+a,offset+b)
                    out[k]=out.get(k,0)+val*expr*sp.Symbol(name)
            channels[name]=out
    return SimpleNamespace(model=model,by_name=channels)


class TestPortalCompletedBasis(unittest.TestCase):
    def test_all16_exact_complement_counts_total_14(self):
        assert len(MODELS)==16
        missing=0
        for model in MODELS:
            with self.subTest(model=model):
                ch=pretend_native(model)
                for species,d in enumerate(DIMS[model[0]],1):
                    aud=exact_portal_completion(ch,species)
                    self.assertEqual(aud['completion_count'],1 if d==3 else 0)
                    self.assertEqual(len(aud['augmented_symbols']),len(aud['rgbeta_symbols']))
                    missing+=aud['completion_count']
        self.assertEqual(missing,14)

    def test_numeric_at_two_scales_and_no_double_count(self):
        ch=pretend_native('E_-2')
        stateF={'lambdaH1':0.14,'lambdaH1Adj':0.02,
                'lambdaH2':0.03,'lambdaH2Adj':-0.05}
        stateS={'lambdaH1':0.015,'lambdaH1Adj':-0.09,
                'lambdaH2':0.07,'lambdaH2Adj':0.012}
        doc,arrays=completed_portal_numeric(ch,{'muF':stateF,'muS':stateS})
        self.assertEqual(doc['baseline_missing_native_directions_filled'],2)
        self.assertEqual(len(arrays['muF_channel_labels']),4)
        self.assertEqual(len(arrays['muS_channel_labels']),4)
        self.assertGreater(np.linalg.norm(arrays['muF_portal_sum']-arrays['muS_portal_sum']),0)
        for scale in ('muF','muS'):
            self.assertTrue(np.allclose(arrays[f'{scale}_channel_tensors'].sum(axis=0),
                                        arrays[f'{scale}_portal_sum'],atol=1e-13))
            self.assertLess(doc['scales'][scale]['tensor_relative_reconstruction_residual'],1e-13)

    def test_missing_coupling_writes_no_deceptive_zero_tensor(self):
        ch=pretend_native('E_-4')
        doc,arr=completed_portal_numeric(ch,{'muF':{'lambdaH1':0.01,'lambdaH2':0.03},
                                              'muS':{'lambdaH1':0.01,'lambdaH2':0.03}})
        self.assertIn('lambdaH1Adj',doc['scales']['muF']['missing_rgbeta_couplings'])
        self.assertNotIn('muF_portal_sum',arr)

    def test_bad_component_embedding_fails_closed(self):
        ch=pretend_native('E_-4')
        name='lambdaH1Inv1'
        ch.by_name[name][(0,0,0,0)]+=sp.Symbol(name)*sp.Rational(1,7)
        with self.assertRaisesRegex(ValueError,'outside RGBeta span'):
            exact_portal_completion(ch,1)

    def test_invalid_output_owner_blocked_before_file_reads(self):
        from pathlib import Path
        with self.assertRaisesRegex(ValueError,'ownership'):
            run_completion_batch(Path('/none'),Path('/mnt/elsewhere'))


class TestExtraReferenceGeometry(unittest.TestCase):
    def test_full_all16_independent_registry_rank_24(self):
        self.assertEqual(sum(len(v) for v in SECTORS.values()),24)
        for model,contents in SECTORS.items():
            with self.subTest(model=model):
                reg={'new_rgbeta_symbols':[
                    {'name':f'Extra{i}','counts':list(c)} for i,c in enumerate(contents)],
                     'preexisting_rgbeta_symbols':[]}
                doc=model_extra_geometry(model,reg,include_components=False)
                self.assertEqual(doc['registry_HHSS_extra_symbols'],len(contents))
                self.assertEqual(sum(x['independent_complex_invariant_count']
                                     for x in doc['monomial_sectors']),len(contents))
                self.assertIsNone(doc.get('physical_coupling_map'))

    def test_complex_channel_both_real_and_imag_tensor_phases(self):
        ty,ds,vs=reference_invariants(2,2,(1,1,0,1,0,1))
        self.assertEqual(len(vs),2)
        for v in vs:
            re,im=real_component_derivatives(2,2,ty,ds,v)
            self.assertTrue(re)
            self.assertTrue(im)
            self.assertNotEqual(re,im)

    def test_fourth_derivative_includes_multiplicity(self):
        ty,ds,vs=reference_invariants(2,2,(0,2,0,0,2,0))
        re,im=real_component_derivatives(2,2,ty,ds,vs[0])
        # Typical entries in this source include 2! from repeated Higgs real indices.
        self.assertTrue(any(row[0]==row[1] for row in re))
        self.assertTrue(any(row[2]==row[3] for row in re))

    def test_bad_hypercharge_extra_is_rejected(self):
        reg={'new_rgbeta_symbols':[{'name':'invalid','counts':[2,0,0,0,2,0]}],
             'preexisting_rgbeta_symbols':[]}
        with self.assertRaisesRegex(ValueError,'nonzero hypercharge'):
            model_extra_geometry('B_-1',reg)

    def test_wrong_monomial_rank_fails(self):
        reg={'new_rgbeta_symbols':[{'name':'only_one','counts':[1,1,0,1,0,1]}],
             'preexisting_rgbeta_symbols':[]}
        with self.assertRaisesRegex(ValueError,'exact SU2 boson invariant rank'):
            model_extra_geometry('B_-1',reg)

    def test_extra_exact_component_conversion_not_a_physical_alias(self):
        reg={'new_rgbeta_symbols':[
             {'name':'Candidate1','counts':[0,2,0,0,2,0]}],
             'preexisting_rgbeta_symbols':[]}
        doc=model_extra_geometry('B_-1',reg)
        row=doc['monomial_sectors'][0]
        self.assertIsNone(row['physical_mapping_to_named_RGBeta_coefficients'])
        self.assertTrue(row['exact_reference_real_component_tensor_basis'][0]
                        ['Re_c_unit_tensor_components'])
        self.assertEqual(doc['numerical_extra_couplings_available'],False)

if __name__=='__main__':
    unittest.main()
