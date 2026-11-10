"""Offline tests for conservative full-general completion/finite C5 gates."""
from __future__ import annotations

import math
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest

import numpy as np

from Numerical.parallel_t3.background.general_batch1_api import encode_native
from Numerical.parallel_t3.complete.full_general_completion import (
    CompletionBlocked, C5_FORMAT, OSCILLATION_KEYS,
    from_complex_json, require_oscillation_target, score_nufit,
    verify_completed_c5_document, put_json,
)
from Numerical.parallel_t3.complete.hard_llhh_source import (
    HardSourceError, extract_inventory, numerical_unsubtracted_source,
)


A4_TERMS=6
D2_TERMS=12


def source(model='A_-4'):
    entries=[]
    for y in ('y2', 'y2Conj') if model=='D_-2' else ('y2',):
        lam='lambdaT3ExtraQ020101I1' if y=='y2Conj' else 'lambdaT3'
        for proj in (-1,1):
            for k in ('const','pole','log'):
                factor={'const':'','pole':'ϵ*','log':''}[k]
                term='(2*hbar*Coupling[y1, {Index[p, Flavor], Index[r, NFlavor]}, 0]*'
                term+=f'Coupling[{y}, {{Index[q, Flavor], Index[r, NFlavor]}}, 0]*'
                term+=f'Coupling[{lam}, {{}}, 0]*'
                if k=='log':term+='Log[μbar2/Coupling[MF, {Index[r, NFlavor]}, 0]^2]*'
                term+=f'NCM[DiracProduct[GammaCC, Proj[{proj}]]]/'
                term+=f'(Sqrt[3]*{factor}Coupling[MF, {{Index[r, NFlavor]}}, 0])'
                entries.append({'sector':'CandidateLeptonNumberViolatingLLHH','raw_inputform':term})
    return {'status':'UnrenormalisedOneLoopTermsEnumerated','terms':entries}


def native(model='A_-4',extra=0):
    n=np.eye(3,dtype=complex)*1e-4
    return {'MF':encode_native(np.diag([1e5,1.1e5,1.2e5]).astype(complex)),
            'y1':encode_native(n),'y2':encode_native(n),
            'lambdaT3':encode_native(complex(.1)),
            'y2Conj':encode_native(n*extra)}


class FullCompletionTests(unittest.TestCase):
    def test_ordered_nufit_targets(self):
        target={'ordering':'NO','dataset':'synthetic',
                'central_values':{k:1. for k in OSCILLATION_KEYS},
                'one_sigma_errors':{k:.01 for k in OSCILLATION_KEYS}}
        self.assertEqual(require_oscillation_target(target)['dataset'],'synthetic')
        answer=score_nufit({k:1. for k in OSCILLATION_KEYS}|{'mixing_defined':True},target)
        self.assertEqual(answer['chi2_gaussian_uncorrelated'],0)
        with self.assertRaises(CompletionBlocked):
            score_nufit({'mixing_defined':False},target)

    def test_invalid_target_sigma(self):
        target={'ordering':'NO','central_values':{k:1. for k in OSCILLATION_KEYS},
                'one_sigma_errors':{k:0 for k in OSCILLATION_KEYS}}
        with self.assertRaises(CompletionBlocked):
            require_oscillation_target(target)

    def test_c5_complex_json(self):
        x=np.array([[1,2j,0],[2j,2,0],[0,0,3]],complex)
        actual=from_complex_json(encode_native(x),label='test')
        np.testing.assert_array_equal(actual,x)
        actual=from_complex_json({'real':x.real.tolist(),'imag':x.imag.tolist()},label='test')
        np.testing.assert_array_equal(actual,x)
        with self.assertRaises(CompletionBlocked):
            from_complex_json([[1,2],[3,4]],label='test')

    def test_source_A_inventory_and_numerics(self):
        checked=extract_inventory('A_-4',source())
        self.assertEqual(len(checked),2)
        out=numerical_unsubtracted_source('A_-4',source(),native())
        self.assertEqual(out['source_term_count'],A4_TERMS)
        self.assertFalse(out['physical_C5_computed'])
        result=np.array([[complex(v['real'],v['imag']) for v in row]
                         for row in out['baseline_ordered_holomorphic_regular_GeVinv']])
        self.assertTrue(np.isfinite(result).all())
        self.assertGreater(abs(result[0,0]),0)
        self.assertEqual(result[0,1],0)

    def test_source_D_includes_extra_yukawa_zero_boundary(self):
        out=numerical_unsubtracted_source('D_-2',source('D_-2'),native('D_-2'))
        self.assertEqual(out['source_term_count'],D2_TERMS)
        self.assertEqual(len(out['source_channel_inventory']),4)

    def test_source_D_nonzero_extra_requires_dictionary(self):
        with self.assertRaises(HardSourceError):
            numerical_unsubtracted_source('D_-2',source('D_-2'),native('D_-2',extra=1))

    def test_source_pole_missing_rejected(self):
        obj=source()
        obj['terms'][1]['raw_inputform']=obj['terms'][0]['raw_inputform']
        with self.assertRaises(HardSourceError):
            extract_inventory('A_-4',obj)

    def test_source_wrong_normalization_rejected(self):
        obj=source()
        obj['terms'][0]['raw_inputform']=obj['terms'][0]['raw_inputform'].replace('Sqrt[3]', 'Sqrt[2]')
        with self.assertRaises(HardSourceError):
            extract_inventory('A_-4',obj)

    def test_source_offdiagonal_fermion_mass_rejected(self):
        mat=native()
        mf=np.diag([1e5,1.1e5,1.2e5]).astype(complex)
        mf[0,1]=2e-3
        mat['MF']=encode_native(mf)
        with self.assertRaises(HardSourceError):
            numerical_unsubtracted_source('A_-4',source(),mat)

    def test_no_forged_matching_boundary(self):
        with TemporaryDirectory() as temp:
            root=Path(temp)
            path=root/'fake.json'
            put_json(path,{'format':C5_FORMAT,'model':'A_-4','C5_MS_GeV_inverse':encode_native(np.eye(3))})
            e=SimpleNamespace(model='A_-4',ms_gev=1e3,
                              sources_sha256={'batch1_manifest':'x','raw_manifest':'y','laurent_manifest':'z'})
            with self.assertRaises(CompletionBlocked):
                verify_completed_c5_document(path,e)

    def test_missing_scalar_source_wxf_rejected(self):
        with TemporaryDirectory() as temp:
            root=Path(temp)
            path=root/'fake.json'
            put_json(path,{
                'format':C5_FORMAT,'model':'A_-4',
                'status':'SourceDerivedFiniteFullGeneralC5AtMS',
                'renormalization_scheme':'MSbar',
                'fermion_threshold_subtraction':'UV_hard_minus_EFT1_soft_plus_counterterms',
                'one_loop_full_EFT1_wilson_running_performed':True,
                'one_loop_joint_scalar_threshold_matching_performed':True,
                'integrated_heavy_fields_in_order':[['F'],['S1','S2']],
                'scale_GeV':1e3,
                'input_sha256':{'batch1_manifest':'x','raw_manifest':'y','laurent_manifest':'z'},
                'symbolic_sources':{'scalar_matching_eft_wxf':{'file':'nonexistent.wxf','sha256':'abc'}},
                'C5_MS_GeV_inverse':encode_native(np.eye(3))})
            e=SimpleNamespace(model='A_-4',ms_gev=1e3,
                              sources_sha256={'batch1_manifest':'x','raw_manifest':'y','laurent_manifest':'z'})
            with self.assertRaises(CompletionBlocked):
                verify_completed_c5_document(path,e)


if __name__=='__main__':
    unittest.main()
