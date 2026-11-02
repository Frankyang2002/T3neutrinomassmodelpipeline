"""Regressions for numerical B/C mixed matching; no live Wolfram required."""
from __future__ import annotations

import json
from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np
import sympy as sp

from Numerical.running.GeneralT3MixedWilsonTensor import (
    split_matchete_sum, factorise_matched_term, build_mixed_source_terms,
    project_mixed_terms, _cg_used,
)
from scripts.RunT3GeneralMixedWilsonTensor import instrument_source, CG_STAGE_ANCHOR


class MixedWilsonTests(unittest.TestCase):
    def test_splits_signed_nested_matchete_plus(self):
        source = '(CG[x, {1, 2}]*GammaCC)/(2*MF) - (CG[y, {1,2}]*GammaCC)/MF + (GammaCC*3)'
        parts=split_matchete_sum(source)
        self.assertEqual(len(parts),3)
        self.assertTrue(parts[1].startswith('-'))
        with self.assertRaises(ValueError):
            split_matchete_sum('(GammaCC + (Broken)')

    def test_instrument_exact_current_source_only(self):
        data='T3MEPerform[]\n llssTerms=Select[something];\n'+CG_STAGE_ANCHOR
        result=instrument_source(data)
        self.assertIn('fermion-first-LLSS-native-CG-registry-exported',result)
        self.assertIn('$CGtensors[name]',result)
        self.assertEqual(result.count(CG_STAGE_ANCHOR),1)
        with self.assertRaises(ValueError):
            instrument_source('changed upstream source')

    def test_reject_unmatched_registry(self):
        terms=[{'chirality':'PL','operator_template_unit_flavor_inputform':'CG[Demo, {}]*GammaCC'}]
        with self.assertRaises(ValueError):
            project_mixed_terms('B_-1',terms,{'format':'bad','model':'B_-1','CGRegistry':[]})

    def test_barred_spurion_recognition(self):
        from Numerical.running.GeneralT3MixedWilsonTensor import _SPURION
        z=complex(0.18,0.24)
        self.assertEqual(_SPURION.findall('Bar[Coupling[T3MEU21, {}, 0]]'), [('Bar[','T3MEU21',']')])
        self.assertAlmostEqual(abs(z),0.3)

    def test_source_chirality_matches_native_runmatching_convention(self):
        from Numerical.running.GeneralT3MixedWilsonTensor import factorise_matched_term
        bare = ("Bar[Coupling[y1, {Index[d$$1, Flavor], Index[d$$9, NFlavor]}, 0]]"
                "*Bar[Coupling[y2, {Index[d$$2, Flavor], Index[d$$9, NFlavor]}, 0]]"
                "*NCM[Transp[Field[l, Fermion, {}, {}]],"
                "DiracProduct[GammaCC, Proj[-1]], Field[l, Fermion, {}, {}]]"
                "/Coupling[MF, {}, 0]")
        conj = ("Coupling[y1, {Index[d$$1, Flavor], Index[d$$9, NFlavor]}, 0]"
                "*Coupling[y2, {Index[d$$2, Flavor], Index[d$$9, NFlavor]}, 0]"
                "*NCM[Bar[Field[l, Fermion, {}, {}]],"
                "DiracProduct[GammaCC, Proj[1]], Transp[Bar[Field[l, Fermion, {}, {}]]]]"
                "/Coupling[MF, {}, 0]")
        self.assertEqual(factorise_matched_term(bare, 1, {})['chirality'], 'PL')
        self.assertEqual(factorise_matched_term(conj, 2, {})['chirality'], 'PR')
        with self.assertRaisesRegex(ValueError, 'projector does not have'):
            factorise_matched_term(bare.replace('Proj[-1]', 'Proj[1]'), 3, {})

    def test_synchronised_benchmark_sources_if_present(self):
        root=Path(__file__).resolve().parents[1]/'output'/'general_observable'
        folder=root/'mixed_f_synchronised_v1'
        if not folder.is_dir():
            self.skipTest('Synchronized F-threshold benchmark outputs not installed')
        for tag,model,count in [('B_m1','B_-1',20),('C_m1','C_-1',24)]:
            with self.subTest(model=model):
                text=(folder/tag/'fermion_tree'/f'{tag}_final_eft1_fermion_first_tree_llss_tree.inputform.txt').read_text()
                basis=json.loads((folder/tag/'mass_basis'/f'{tag}_final_mass_basis.json').read_text())
                run=json.loads((root/'mixed_adapter_v3'/tag/'running.json').read_text())
                def matrix(name):
                    return np.array([[complex(v['real'],v['imag']) for v in r]
                                     for r in run['uv_at_fermion_threshold'][name]])
                mass=matrix('MF')
                self.assertLess(np.max(abs(mass-np.diag(np.diag(mass)))),1e-9)
                self.assertLess(np.max(abs(np.diag(mass).imag)),1e-9)
                terms=build_mixed_source_terms(model,text,basis,np.diag(mass).real,
                                               {'y1':matrix('y1'),'y2':matrix('y2')},
                                               expected_count=count)
                self.assertEqual(len(terms),count)
                self.assertEqual(sum(x['chirality']=='PL' for x in terms),count//2)
                self.assertEqual(sum(x['chirality']=='PR' for x in terms),count//2)
                self.assertTrue(all('Proj[-1]' in t['source_term_inputform']
                                    for t in terms if t['chirality']=='PL'))
                self.assertTrue(all('Proj[1]' in t['source_term_inputform']
                                    for t in terms if t['chirality']=='PR'))
                self.assertEqual(_cg_used(terms),
                    {'T3MENativeY1CG','T3MENativeY2CG'} if model=='B_-1' else
                    {'T3ExtraCGQ020020I1','T3MEBridgeBB','eps[SU2L]'})
                # Test only parser and basis mapping with CG STUBS. These
                # cannot establish numerical Matchete CG phases/norms.
                cg_names=_cg_used(terms)
                fake_registry=[{'Name':name} for name in sorted(cg_names)]
                with patch('RGE.running.intermediate.ScalarOnlyTensorAdapters.load_wilson_cg_registry',
                           return_value={x:object() for x in cg_names}), \
                     patch('RGE.running.intermediate.ScalarOnlyTensorAdapters._wilson_cg_value',
                           return_value=sp.S.One):
                    projected=project_mixed_terms(model,terms,{
                        'format':'t3-mixed-matchete-native-cg-registry-v1',
                        'model':model,'CGRegistry':fake_registry})
                self.assertEqual(len(projected),count//2)
                self.assertGreater(sum(x['component_count'] for x in projected),0)


if __name__=='__main__':
    unittest.main()
