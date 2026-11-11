"""Synthetic structural tests for provisional Matchete overlay."""
import json
import tempfile
import unittest
from pathlib import Path
import numpy as np
from Numerical.running.T3ExploratoryFullQuartic import build_overlay, ExploratoryQuarticError

class TestOverlay(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        p=Path(self.tmp.name)
        self.registry=p/'registry.json';self.candidates=p/'candidates.json'
        self.registry.write_text(json.dumps({'models':{'A_+0':{
            'new_rgbeta_symbols':[
                {'name':'lambdaExtraX','counts':[0,0,2,2,0,0]},
                {'name':'lambdaExtraH','counts':[1,1,1,1,0,0]}],
            'preexisting_rgbeta_symbols':[], 'total_additional_quartics':2}}}))
        def rec(name, entries):
            return {'symbol':name,'cartesian_dimensions':12,
                'responses':[{'component':'real','entries':entries},
                             {'component':'imag','entries':[]}]}
        records=[rec('lambdaT3ExtraQ002200I1',
                      [{'indices':[4,4,5,5], 'lambda_abcd':'8'}]),
                 rec('lambdaT3ExtraQ111100I1',
                      [{'indices':[0,0,4,4], 'lambda_abcd':'7'}])]
        self.candidates.write_text(json.dumps({
            'format':'t3-matchete-candidate-cartesian-quartics-v1',
            'rgbeta_mapping_certified':False,
            'models':{'A_+0':{'records':records}}}))
        self.hhss=lambda state: np.zeros((4,4,8,8))
        self.baseline=lambda model,state,h: np.zeros((12,)*4)
    def tearDown(self): self.tmp.cleanup()
    def build(self):
        return build_overlay('A_+0',self.candidates,self.registry,
                             self.baseline,self.hhss,
                             required_quartics=('lambdaExtraX','lambdaExtraH'))
    def test_nonhhss_changes_with_running(self):
        basis=self.build()
        a=basis.evaluate({'lambdaExtraX':0.3,'lambdaExtraH':2.})
        b=basis.evaluate({'lambdaExtraX':0.6,'lambdaExtraH':2.})
        self.assertAlmostEqual(a[4,4,5,5],2.4)
        self.assertAlmostEqual(b[4,4,5,5],4.8)
    def test_hhss_not_double_counted(self):
        basis=self.build()
        a=basis.evaluate({'lambdaExtraX':1.,'lambdaExtraH':5.})
        self.assertEqual(a[0,0,4,4],0.)
        self.assertEqual(basis.validate_hhss({'lambdaExtraX':1.,'lambdaExtraH':5.},
                                                self.hhss({})),0.)
    def test_mixed_snapshot_nested_complex_running(self):
        basis=self.build()
        z={'lambdaExtraX':{'real':0.3,'imag':0.25},
           'lambdaExtraH':{'real':2.0,'imag':0.0}}
        one=basis.evaluate(z)
        self.assertAlmostEqual(one[4,4,5,5],2.4)
        z['lambdaExtraX']={'real':0.7,'imag':0.25}
        two=basis.evaluate(z)
        self.assertAlmostEqual(two[4,4,5,5],5.6)

    def test_malformed_nested_complex_fails(self):
        basis=self.build()
        with self.assertRaises(ExploratoryQuarticError):
            basis.evaluate({'lambdaExtraX':{'real':1},
                            'lambdaExtraH':2.})

    def test_fails_missing_running_coupling(self):
        with self.assertRaises(ExploratoryQuarticError):
            self.build().evaluate({'lambdaExtraH':1.})
    def test_fails_ambiguous_basis(self):
        x=json.loads(self.candidates.read_text())
        x['models']['A_+0']['records'][0]['symbol']='lambdaT3ExtraQ002200I2'
        self.candidates.write_text(json.dumps(x))
        with self.assertRaises(ExploratoryQuarticError):self.build()
    def test_refuses_certified_mislabel(self):
        x=json.loads(self.candidates.read_text());x['rgbeta_mapping_certified']=True
        self.candidates.write_text(json.dumps(x))
        with self.assertRaises(ExploratoryQuarticError):self.build()

if __name__=='__main__':unittest.main()
