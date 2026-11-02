"""No physical T3 matrices are certified by these synthetic adapter checks."""
import json
import tempfile
import unittest
from pathlib import Path
import numpy as np
import sympy as sp
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from RunT3DirectOneLoop import models, CLASSES
from BuildT3DirectOneLoopInputs import construct
from Numerical.direct_one_loop.core import evaluate


class AdapterChecks(unittest.TestCase):
    def test_all_16_catalog(self):
        items=list(models())
        self.assertEqual(len(items),16)
        self.assertEqual(len({i for i,_ in items}),16)
        self.assertEqual([len(z[3]) for z in CLASSES.values()],[3,3,3,3,4])

    def test_explicitly_uncertified_cg_input(self):
        # Only nonzero entries relevant to the neutral B_-1 subdiagram; not
        # a source CG registry or the actual physical benchmark.
        case={'model':{'dS1':2,'dS2':2,'dF':1,'alpha':-1},'vertices':{
            'T3Y1CG':{'dimensions':[2,1,2],
                'components':[{'indices':[1,1,1],'cg':'1'}]},
            'T3Y2CG':{'dimensions':[2,1,2],
                'components':[{'indices':[1,1,2],'cg':'1'}]},
            'T3MixCG':{'dimensions':[2,2,2,2],
                'components':[{'indices':[2,2,1,2],'cg':'1'}]}}}
        y={'real':(np.ones((3,3))*0.0001).tolist(),'imag':np.zeros((3,3)).tolist()}
        bench={'representation':{'d_s1':2,'d_s2':2,'d_f':1,'alpha':-1},
               'numerical':{'vev_gev':246.22},
               'base_state':{'mu_gev':1e7,'ordinary':{
                    'y1':y,'y2':y,
                    'MF':{'real':np.diag([1e5,1.1e5,1.2e5]).tolist(),
                          'imag':np.zeros((3,3)).tolist()},
                    'mS1Sq':1306122.44898,'mS2Sq':1e6,
                    'lambdaH1':0.02,'lambdaH2':0.03,
                    'lambdaT3':{'real':.1,'imag':0}}}}
        with tempfile.TemporaryDirectory() as td:
            source=Path(td)/'synthetic.json';source.write_text(json.dumps(bench))
            obj=construct(case,bench,{1:sp.eye(1),2:sp.eye(2),3:sp.eye(3)},
                          'B_-1',source,source,source,True)
            self.assertFalse(obj['provenance']['validated_uv_component_conventions'])
            result=evaluate(obj,retain_individual_terms=False)
            self.assertTrue(result['diagnostics']['uv_finite'])
            self.assertGreater(max(result['total']['neutrino_masses_ev']),0.)
            obj['scalar_A_gev2']['real'][0][3]=0.
            obj['scalar_A_gev2']['real'][3][0]=0.
            zero=evaluate(obj,retain_individual_terms=False)
            self.assertLess(zero['diagnostics']['total_frobenius_gev'],1e-17)

if __name__=='__main__': unittest.main()

class RunnerAuditTests(unittest.TestCase):
    def test_missing_physical_components_never_fabricate_matrices(self):
        from RunT3DirectOneLoop import run
        from tempfile import TemporaryDirectory
        with TemporaryDirectory() as td:
            root=Path(td)
            bench=root/'Numerical'/'benchmarks'/'fitted_reference'
            bench.mkdir(parents=True)
            for short,key in models():
                rep=CLASSES[short[0]]
                data={'representation':{'d_s1':rep[0],'d_s2':rep[1],
                                        'd_f':rep[2],'alpha':int(short[2:])}}
                (bench/(key+'.json')).write_text(json.dumps(data))
            dest=root/'output'/'general_observable'/'lo_direct_one_loop'
            data=run(root,root/'input'/'lo_direct_one_loop',dest)
            self.assertEqual(data['count'],16)
            self.assertEqual(data['evaluated'],0)
            self.assertEqual(len(list(dest.glob('*.json'))),1)
            self.assertTrue(all(r['status']=='input_missing' for r in data['models']))
