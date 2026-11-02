"""Unit tests for v4 source-backed HHSS channel extraction and muS handoff."""
from __future__ import annotations
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import sympy as sp

from Numerical.joint_scalar_matching.fermion_boundary_bridge import SourceBoundary
from Numerical.joint_scalar_matching.large_batch import load_background
from Numerical.joint_scalar_matching.source_quartic_v4 import (
    QuarticChannels,build_source_channels,evaluate_channels,
    build_background,source_seed_path,symbolic_export,ENDPOINT_FORMAT)


def fake_source(model='A_-4',special=False):
    label1,label2=('NewScalarX','NewScalarY') if special else ('NewScalar1','NewScalar2')
    return SourceBoundary(model,Path('synthetic-source.json'),'h'*64,1e5,
                          np.zeros((3,3,2,2,8,8),complex),
                          np.zeros((1,3,3,2,2,8,8),complex),
                          (2,),tuple([label1]*2+[label2]*6),False,'test-source')


def channels(path):
    x=sp.Symbol('x')
    y=sp.Symbol('y')
    return QuarticChannels('A_-4',path,'f'*64,
                           tuple(['NewScalar1']*2+['NewScalar2']*6),
                           {'x':{(2,2,0,0):x+sp.conjugate(x)},
                            'y':{(0,0,0,1):2*y}},2,3,'fake verified tensor parser')


def endpoint(source,extra=None):
    base={'format':ENDPOINT_FORMAT,'model':'A_-4',
          'seed_quartic_sha256':'f'*64,'source_boundary_sha256':source.hash,
          'renormalization_scheme':'MSbar','scalar_mass_origin':'unbroken_electroweak',
          'both_scalar_multiplets_active_until_single_muS':True,
          'masses_and_quartics_at_same_muS':True,
          'mu_f_gev':1e5,'mu_s_gev':1e3,
          'mass_squared_gev2_original_real_basis_muS':(np.eye(8)*1e6).tolist(),
          'quartic_couplings_muS':{'x':{'real':0.5,'imag':0.3},'y':0.2},
          'hhss_source_operator_set_reviewed':True,
          'provenance':{'scalar_coupling_trajectory':'independent muS couplings',
                        'scalar_mass_trajectory':'independent muS masses',
                        'hhss_source_review':'source review'}}
    if extra:base.update(extra)
    return base


class V4Tests(unittest.TestCase):
    def test_source_path_and_independent_channels(self):
        self.assertTrue(str(source_seed_path(Path('output/full16_symbolic'),'A_-4'))
                        .endswith('T3_A_alpha_m4/data/eft1_after_F_scalar_quartic_seed.json'))
        x=sp.Symbol('x')
        class FakeSparse:
            def __init__(self,entries):self.entries=entries
            def nonzero_items(self):return iter(self.entries.items())
        def backend(seed,**kwargs):
            self.assertEqual(kwargs,{'d_s1':1,'d_s2':3,'shared_scalar':False})
            coupling='x' if 'x' in seed['ScalarQuarticTerms'][0]['TermInputForm'] else 'y'
            return FakeSparse({(3,3,5,5):sp.Symbol(coupling),
                               (1,1,1,1):sp.Symbol(coupling)})
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/'source.json'
            seed={'ScalarQuarticTermCount':3,'IntegratedFields':['F'],
                  'ActiveHeavyFields':['S1','S2'],'CGRegistry':[],
                  'ScalarQuarticTerms':[
                      {'HiggsFieldCount':2,'Scalar1FieldCount':2,'Scalar2FieldCount':0,
                       'TermInputForm':'Coupling[x, {}, 0]*rest'},
                      {'HiggsFieldCount':2,'Scalar1FieldCount':1,'Scalar2FieldCount':1,
                       'TermInputForm':'Coupling[y, {}, 0]*other'},
                      {'HiggsFieldCount':4,'Scalar1FieldCount':0,'Scalar2FieldCount':0,
                       'TermInputForm':'Coupling[z, {}, 0]*higgs'}]}
            path.write_text(json.dumps(seed))
            c=build_source_channels(path,'A_-4',builder=backend)
            self.assertEqual(set(c.by_name),{'x','y'})
            self.assertEqual(c.seed_hhss_term_count,2)
            self.assertEqual(c.by_name['x'][(2,2,0,0)],x)
            self.assertEqual(len(c.by_name['x']),1) # H^4 not copied as HHSS
            self.assertEqual(symbolic_export(c,fake_source())['source_basis_matches_fermion_boundary'],True)

    def test_numeric_couplings_real_potential_and_explicit_zero(self):
        c=channels(Path('example'))
        t,n=evaluate_channels(c,{'x':{'real':0.5,'imag':0.3},'y':0})
        self.assertEqual(n,('x','y'))
        self.assertAlmostEqual(t[0,2,2,0,0],1.0)
        self.assertEqual(t[1,0,0,0,1],0.0)
        self.assertEqual(t[1,0,0,1,0],0.0)
        with self.assertRaisesRegex(ValueError,'missing muS quartic'):
            evaluate_channels(c,{'x':0.5})
        with self.assertRaisesRegex(ValueError,'component is complex'):
            evaluate_channels(QuarticChannels('A_-4',Path('x'),'f'*64,c.source_basis_labels,
                                {'x':{(0,0,0,0):sp.Symbol('x')}},1,1,'mock'),
                              {'x':{'real':0,'imag':0.5}})

    def test_ordinary_background_writes_v3_compatible_npz(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            ch=channels(root/'source.json');src=fake_source()
            ep=root/'muS.json';ep.write_text(json.dumps(endpoint(src)))
            r=build_background(ch,src,ep,root/'backgrounds')
            bg=load_background(root/'backgrounds'/'A_m4.json',src)
            self.assertEqual(bg.quartic.shape,(4,4,8,8))
            self.assertAlmostEqual(bg.quartic[2,2,0,0],1.)
            self.assertAlmostEqual(bg.quartic[0,0,0,1],.4)
            self.assertTrue(np.allclose(bg.mass2,np.eye(8)*1e6))
            self.assertFalse(r['physical_C5_certified'])
            with self.assertRaisesRegex(ValueError,'refusing to overwrite'):
                build_background(ch,src,ep,root/'backgrounds')

    def test_mixed_split_needs_provenance_and_explicit_orthogonal_map(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);ch=channels(root/'seed.json');src=fake_source(special=True)
            ep=root/'endpoint.json'
            ep.write_text(json.dumps(endpoint(src)))
            with self.assertRaisesRegex(ValueError,'independently sourced'):
                build_background(ch,src,ep,root/'out1')
            x=endpoint(src);x['provenance']['original_to_boundary_basis_map']='separate external SU2/phase CG calculation'
            x['special_scalar_basis_transform_externally_validated']=True
            x['original_real_equals_boundary_real_R']=(np.eye(8)*2).tolist()
            ep.write_text(json.dumps(x))
            with self.assertRaisesRegex(ValueError,'orthogonal'):
                build_background(ch,src,ep,root/'out2')
            R=np.eye(8);R[[0,1]]=R[[1,0]]
            x['original_real_equals_boundary_real_R']=R.tolist()
            ep.write_text(json.dumps(x))
            build_background(ch,src,ep,root/'out3')
            with np.load(root/'out3'/'A_m4.npz') as z:
                self.assertAlmostEqual(z['lambda_HHSS_muS'][2,2,1,1],1.)
                self.assertAlmostEqual(z['lambda_HHSS_muS'][0,0,0,1],.4)
                self.assertTrue(np.allclose(z['source_to_matching_real_U'],np.eye(8)))

    def test_uv_not_substituted_for_muS_and_missing_provenance(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);ch=channels(root/'s');src=fake_source();ep=root/'ep.json'
            e=endpoint(src,{'masses_and_quartics_at_same_muS':False})
            ep.write_text(json.dumps(e))
            with self.assertRaisesRegex(ValueError,'common muS'):
                build_background(ch,src,ep,root/'out')
            e=endpoint(src);e['provenance'].pop('scalar_coupling_trajectory')
            ep.write_text(json.dumps(e))
            with self.assertRaisesRegex(ValueError,'mass AND coupling'):
                build_background(ch,src,ep,root/'out')
            e=endpoint(src);e['mass_squared_gev2_original_real_basis_muS'][0][0]=-1e6
            ep.write_text(json.dumps(e))
            with self.assertRaises((RuntimeError,ValueError)):
                build_background(ch,src,ep,root/'out')


if __name__=='__main__':unittest.main()
