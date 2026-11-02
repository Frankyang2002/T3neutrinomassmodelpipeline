"""Run with: python -X utf8 -m unittest discover -s tests -p test_t3_general_wilson_boundary.py"""
from __future__ import annotations

import importlib.util
from pathlib import Path
import tempfile
import json
from unittest.mock import patch
import unittest

import numpy as np

SOURCE = Path(__file__).resolve().parents[1] / 'scripts' / 'RunT3GeneralWilsonBoundary.py'
spec = importlib.util.spec_from_file_location('RunT3GeneralWilsonBoundary', SOURCE)
mod = importlib.util.module_from_spec(spec)
import sys
sys.modules[spec.name] = mod
spec.loader.exec_module(mod)

P = 'd$$1'
Q = 'd$$2'
R = 'd$$1'

def y(name, ext, bar=True):
    expr = f'Coupling[{name}, {{Index[{ext}, Flavor], Index[{R}, NFlavor]}}, 0]'
    return f'Bar[{expr}]' if bar else expr

def matchete(y_a='y1',y_b='y2',bar=True):
    return (f'({y(y_a,P,bar)}*{y(y_b,Q,bar)}*'
            'CG[T3MixCG, {Index[d$$1, SU2L[fund]]}]*'
            'Field[NewScalar1, Scalar, {}, {}]*Field[NewScalar2, Scalar, {}, {}]*'
            'NCM[DiracProduct[GammaCC, Proj[-1]]])/Coupling[MF, {}, 0]')

class BoundaryTest(unittest.TestCase):
    def test_16_model_registry(self):
        self.assertEqual(len(mod.MODELS),16)
        self.assertEqual(len(set(mod.MODELS)),16)

    def test_factorisation_preserves_source_scalar_cg_and_phase(self):
        expr=matchete()
        x=mod._normalise_term(expr)
        template=x['operator_template_unit_flavor_inputform']
        self.assertIn('CG[T3MixCG',template)
        self.assertIn('Field[NewScalar2',template)
        self.assertNotIn('Coupling[y1',template)
        self.assertNotIn('Coupling[MF',template)
        self.assertEqual([a['name'] for a in x['factors']],['y1','y2'])
        self.assertTrue(all(a['conjugated'] for a in x['factors']))

    def test_complex_flavour_matrix_preserves_bar(self):
        x=mod._normalise_term(matchete())
        a=np.eye(3,dtype=complex)*(1+2j)
        b=np.eye(3,dtype=complex)*(3+4j)
        m=np.array([2.,4.,8.])
        z=mod.flavour_matrix(x,m,{'y1':a,'y2':b})
        expected=(1-2j)*(3-4j)/m
        self.assertTrue(np.allclose(z,np.diag(expected)))

    def test_nonbar_yukawa_and_extra_channel(self):
        x=mod._normalise_term(matchete('y1','y2Conj',bar=False))
        a=np.diag([1+1j,2+1j,3+1j])
        b=np.diag([1-2j,2-2j,3-2j])
        z=mod.flavour_matrix(x,np.array([10.,20.,30.]),{'y1':a,'y2Conj':b})
        self.assertAlmostEqual(z[0,0],(1+1j)*(1-2j)/10)
        self.assertTrue(all(not f['conjugated'] for f in x['factors']))

    def test_bad_heavy_dummy_blocked(self):
        expr=matchete().replace(f'Index[{Q}, Flavor], Index[{R}, NFlavor]',
                               f'Index[{Q}, Flavor], Index[d$$5, NFlavor]')
        with self.assertRaisesRegex(ValueError,'common heavy'):
            mod._normalise_term(expr)

    def test_missing_source_yukawa_blocked(self):
        expr=matchete().replace('Coupling[y2,','Coupling[y3,')
        with self.assertRaisesRegex(ValueError,'two source Yukawa'):
            mod._normalise_term(expr)

    def test_mass_denominator_rejected_when_ambiguous(self):
        expr=matchete().replace(')/Coupling[MF, {}, 0]',')*Coupling[MF, {}, 0]')
        with self.assertRaisesRegex(ValueError,'mass-denominator'):
            mod._normalise_term(expr)

    def test_real_e_family_sqrt3_mass_denominator_preserved(self):
        source=matchete().replace(')/Coupling[MF, {}, 0]',')/(Sqrt[3]*Coupling[MF, {}, 0])')
        template=mod._normalise_term(source)['operator_template_unit_flavor_inputform']
        self.assertIn('/(Sqrt[3]*1)',template)

    def test_seed_inventory_and_chirality(self):
        x=matchete()
        payload={'IntegratedFields':['F'],'ActiveHeavyFields':['S1','S2'],
                 'TreeWilsonTermCount':2,
                 'TreeWilsonTerms':[{'Index':1,'Chirality':'PL','CGNames':['T3MixCG'],'TermInputForm':x},
                                    {'Index':2,'Chirality':'PR','CGNames':['T3MixCG'],'TermInputForm':x.replace('Proj[-1]','Proj[1]')}]}
        items=mod._validate_seed(payload,'A_-4')
        self.assertEqual(len(items),1)
        self.assertEqual(items[0]['source_index'],1)

    def test_broken_seed_blocked(self):
        with self.assertRaises(ValueError):
            mod._validate_seed({'IntegratedFields':['F'], 'ActiveHeavyFields':['S1','S2'],
                                'TreeWilsonTermCount':3,'TreeWilsonTerms':[]},'B_+1')

    def test_nonmixed_file_route_with_mocked_heavy_basis(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            folder=root/'source'/'A_m4'
            (folder/'data').mkdir(parents=True)
            (folder/'status.json').write_text(json.dumps({
                'model':'A_-4','status':'FermionThresholdTreeWilsonExported'
            }),encoding='utf8')
            example={'IntegratedFields':['F'],'ActiveHeavyFields':['S1','S2'],
                     'TreeWilsonTermCount':1,
                     'TreeWilsonTerms':[{'Index':1,'Chirality':'PL','CGNames':['T3MixCG'],
                                         'TermInputForm':matchete()}],
                     'CGRegistry':[{'Name':'T3MixCG','TensorInputForm':'not-evaluated-here'}]}
            (folder/'data'/'eft1_after_F_wilson_seed.json').write_text(json.dumps(example))
            run=root/'running'/'A_m4'
            run.mkdir(parents=True)
            (run/'running.json').write_text(json.dumps({'model':'A_-4'}))
            yy={'y1':np.eye(3,dtype=complex)*0.02,
                'y2':np.eye(3,dtype=complex)*0.03}
            with patch.object(mod,'ROOT',root),patch.object(mod,'_running_basis',
                    return_value=({'scales_gev':{'F':1e5}},np.array([1.5e5,1.8e5,2e5]),yy,{})):
                row=mod.process('A_-4',(Path('source'),),(Path('running'),),root/'out')
            self.assertTrue(row['status'].startswith('SourceFactored'))
            result=mod.read(root/'out'/'A_m4'/'wilson_boundary.json')
            self.assertEqual(result['ordered_flavour_factor_count'],1)
            self.assertTrue(result['numeric_flavour_contractions_computed'])
            self.assertFalse(result['eft1_wilson_rge_computed'])

    def test_mixed_is_source_only_not_misrepresented_as_numeric(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            folder=root/'source'/'B_m1'
            folder.mkdir(parents=True)
            raw=root/'mixed_llss.txt'
            raw.write_text('Source Matchete operator with NewScalarX and NewScalarY')
            (folder/'status.json').write_text(json.dumps({
                'model':'B_-1','status':'FermionThresholdMixedMassTreeLLSSExported',
                'tree_wilson_term_count':20,'llss_operator_file':str(raw),
                'mass_basis':'validated experimental mixed-doublet X,Y'
            }))
            with patch.object(mod,'ROOT',root):
                row=mod.process('B_-1',(Path('source'),),(),root/'out')
            self.assertTrue(row['status'].startswith('VerifiedMixed'))
            info=mod.read(root/'out'/'B_m1'/'wilson_boundary.json')
            self.assertFalse(info['numerical_wilson_boundary_computed'])
            self.assertEqual(info['operator_term_count'],20)

    def test_encoded_matrix_full_3x3(self):
        j=mod.encode_matrix(np.eye(3,dtype=complex)*1e-5)
        self.assertEqual(len(j['real']),3)
        self.assertAlmostEqual(j['frobenius_norm'],np.sqrt(3)*1e-5)

if __name__=='__main__': unittest.main()
