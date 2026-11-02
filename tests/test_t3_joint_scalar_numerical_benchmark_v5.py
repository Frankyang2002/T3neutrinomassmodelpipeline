"""Tests for v5 source-backed exploratory numerical matching, no physical claims."""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from dataclasses import replace

import numpy as np
import sympy as sp

from Numerical.joint_scalar_matching.fermion_boundary_bridge import SourceBoundary
from Numerical.joint_scalar_matching.source_quartic_v4 import QuarticChannels,FORMAT,POTENTIAL_CONVENTION
from Numerical.joint_scalar_matching.numerical_benchmark_v5 import (
    ALIASES,REGULAR,SPECIAL,read_v4_cached_channels,load_running_snapshot,
    select_candidate_couplings,inspect_model,build_ordinary_background,
    evaluate_ordinary_response)
from Numerical.joint_scalar_matching.large_batch import load_background


def src(model='A_-4',special=False):
    labels=['NewScalar1']*2+['NewScalar2']*6
    if special:
        labels=['NewScalar1']*2+['T3Scalar2Minus']*3+['T3Scalar2Plus']*3
    a=np.zeros((3,3,2,2,8,8),complex)
    a[0,0,0,0,0,2]=a[0,0,0,0,2,0]=1e-12
    a[0,0,1,1,0,2]=a[0,0,1,1,2,0]=-1e-12
    return SourceBoundary(model,Path('mock-source'),'h'*64,1e5,a,np.stack([a]),
                          (2,),tuple(labels),False,'provisional')


def quartic(model='A_-4'):
    n=8
    labels=('NewScalar1',)*2+('NewScalar2',)*6
    # An actual real-potential neutral+charged HHSS tensor component pattern.
    x=sp.Symbol('lambdaT3')
    entries={(0,0,0,2):x,(1,1,0,2):x,(2,2,0,2):x,(3,3,0,2):x}
    return QuarticChannels(model,Path('seed'), 'a'*64, labels,
                           {'lambdaT3':entries},1,1,'mock source parser')


def running_doc(model='A_-4'):
    return {'model':model,'status':'NumericalRenormalisableRunningOnly',
            'illustrative_benchmark':True,
            'physical_mass_parameters_available':False,
            'eft1_wilson_running_performed':False,
            'scales_gev':{'UV':1e7,'F':1e5,'S':1e3},
            'eft1_at_fermion_threshold':{'mS1Sq':4e6,'mS2Sq':5e6,'lambdaT3':0.9},
            'eft1_at_scalar_threshold':{'mS1Sq':1e6,'mS2Sq':9e6,
                                        'lambdaT3':{'real':.3,'imag':0},
                                        'lambdaH1':.1,'lambdaH2Adj':.2},
            'symbolic_provenance':{'physical_matchete_rgbeta_mapping_certified':False}}


def load_doc(root,model='A_-4',updates=None):
    doc=running_doc(model)
    if updates:doc.update(updates)
    p=root/'running.json';p.write_text(json.dumps(doc));return load_running_snapshot(p,model)


class V5Tests(unittest.TestCase):
    def test_eight_ordinary_and_eight_special_complete(self):
        self.assertEqual(len(REGULAR),8)
        self.assertEqual(len(SPECIAL),8)
        self.assertFalse(REGULAR & SPECIAL)
        self.assertEqual(len(REGULAR|SPECIAL),16)

    def test_snapshot_uses_eft1_at_scalar_threshold_not_earlier(self):
        with tempfile.TemporaryDirectory() as td:
            x=load_doc(Path(td))
            self.assertEqual(x.scalar_mass1_sq,1e6)
            self.assertEqual(x.scalar_mass2_sq,9e6)
            self.assertEqual(x.couplings['lambdaT3']['real'],.3)
            self.assertEqual(x.mu_f,1e5)
            self.assertEqual(x.mu_s,1e3)

    def test_rejects_uv_fallback_and_unverified_nonillustrative_input(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            for edits in ({'eft1_at_scalar_threshold':None},
                          {'illustrative_benchmark':False},
                          {'scales_gev':{'UV':1e7,'F':1e5,'S':1e5}}):
                with self.assertRaises(ValueError):load_doc(root,updates=edits)

    def test_every_named_channel_requires_explicit_numerical_value(self):
        with tempfile.TemporaryDirectory() as td:
            snap=load_doc(Path(td))
            ch=quartic()
            values,mapping=select_candidate_couplings(ch,snap)
            self.assertEqual(mapping,{'lambdaT3':'lambdaT3'})
            self.assertEqual(values['lambdaT3']['real'],.3)
            bad=QuarticChannels(ch.model,ch.seed_path,ch.seed_sha256,ch.source_basis_labels,
                                {'unknown_new_quartic':{(0,0,0,0):sp.Symbol('unknown_new_quartic')}},1,1,'mock')
            with self.assertRaisesRegex(ValueError,'UNMAPPED'):
                select_candidate_couplings(bad,snap)
            self.assertEqual(ALIASES["lambdaH2Inv2"],"lambdaH2Adj")

    def test_ordinary_background_roundtrip_and_conditional_matrix(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            # The toy tensor has no physical SU2 Ward completion. Use zero
            # Wilson source to test the numerical/metadata transport only.
            # Nonzero SU2-complete responses are already covered by v2/v3.
            raw=src()
            boundary=replace(raw,raw=np.zeros_like(raw.raw),
                             by_source=np.zeros_like(raw.by_source))
            ch=quartic()
            snapshot=load_doc(root)
            info=inspect_model(ch,boundary,snapshot)
            self.assertFalse(info['special_mixing_or_split_required'])
            out=build_ordinary_background(ch,boundary,snapshot,root/'output')
            bg=load_background(Path(out['background_json']),boundary)
            self.assertTrue(np.allclose(np.diag(bg.mass2),[1e6]*2+[9e6]*6))
            self.assertEqual(bg.mu_s,1000)
            self.assertEqual(bg.quartic.shape,(4,4,8,8))
            answer=evaluate_ordinary_response(boundary,Path(out['background_json']))
            self.assertFalse(answer['physical_C5_certified'])
            self.assertIsNone(answer['source_to_master_spinor_factor'])
            self.assertTrue(np.isfinite(bg.mass2).all())
            with self.assertRaisesRegex(ValueError,'refusing to overwrite'):
                build_ordinary_background(ch,boundary,snapshot,root/'output')

    def test_nonzero_su2_complete_conditional_response(self):
        # Reuse the independently Ward-complete LLSS/Q construction from v2,
        # unlike the deliberately incomplete toy tensor in the roundtrip test.
        from test_t3_joint_scalar_integrated_batch import clean_synthetic_source
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            source,Q,_=clean_synthetic_source('A_-4',root/'synthetic')
            source=replace(source,labels=tuple(['NewScalar1']*2+['NewScalar2']*6))
            x=sp.Symbol('lambdaT3')
            entries={}
            for a in range(4):
                for b in range(a,4):
                    for c in range(8):
                        for d in range(c,8):
                            if Q[a,b,c,d]!=0:
                                entries[(a,b,c,d)]=float(Q[a,b,c,d])/.3*x
            self.assertGreater(len(entries),2)
            ch=QuarticChannels('A_-4',Path('source_seed'), 'f'*64,
                 tuple(['NewScalar1']*2+['NewScalar2']*6),
                 {'lambdaT3':entries},1,1,'constructed Ward-complete unit fixture')
            snap=load_doc(root)
            snap.couplings['mS2Sq']=1e6
            # Re-read the modified endpoint, not a new synthetic scalar running.
            snap=replace(snap,scalar_mass2_sq=1e6)
            out=build_ordinary_background(ch,source,snap,root/'backgrounds')
            actual=evaluate_ordinary_response(source,Path(out['background_json']))
            from Numerical.joint_scalar_matching.large_batch import parse_matrix
            C=parse_matrix(actual['response']['C5_LLSS_total_one_loop_per_N_gev_inverse'])
            self.assertGreater(float(np.linalg.norm(C)),1e-20)
            self.assertLess(actual['response']['weinberg_beta_projection_residual'],1e-8)
            self.assertFalse(actual['physical_C5_certified'])

    def test_mixed_split_blocked_without_real_basis_conversion(self):
        with tempfile.TemporaryDirectory() as td:
            snap=load_doc(Path(td),model='A_-2')
            boundary=src('A_-2',special=True)
            ch=quartic('A_-2')
            status=inspect_model(ch,boundary,snap)
            self.assertTrue(status['special_mixing_or_split_required'])
            with self.assertRaisesRegex(ValueError,'ordinary ORIGINAL=BOUNDARY'):
                build_ordinary_background(ch,boundary,snap,Path(td)/'out')

    def test_cached_quartics_require_current_hash_and_correct_basis(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            seed=root/'seed.json';seed.write_text('original seed bytes')
            import hashlib
            dig=hashlib.sha256(seed.read_bytes()).hexdigest()
            b=src()
            doc={'format':FORMAT,'model':'A_-4','source_boundary_sha256':b.hash,
                 'seed_quartic_sha256':dig,
                 'source_real_scalar_labels':['NewScalar1']*2+['NewScalar2']*6,
                 'fermion_boundary_real_scalar_labels':list(b.labels),
                 'potential_tensor_normalization':POTENTIAL_CONVENTION,
                 'couplings_required_at_muS':['lambdaT3'],
                 'hhss_native_terms':1,'all_quartic_native_terms':2,
                 'channels':{'lambdaT3':[
                    {'component_zero_based_HA_HB_Sa_Sb':[2,2,0,2],
                     'exact_expr':'3**(3/4)*(-lambdaT3 - conjugate(lambdaT3))/6'}]}}
            cache=root/'v4.json';cache.write_text(json.dumps(doc))
            ch=read_v4_cached_channels(cache,seed,b)
            self.assertIn('lambdaT3',ch.by_name)
            self.assertTrue(ch.by_name['lambdaT3'][(2,2,0,2)].has(sp.Symbol('lambdaT3')))
            seed.write_text('tampered')
            with self.assertRaisesRegex(ValueError,'hashes'):
                read_v4_cached_channels(cache,seed,b)


if __name__=='__main__':unittest.main()
