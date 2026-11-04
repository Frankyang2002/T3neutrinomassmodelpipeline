"""Unit tests of uniform-0.01 UV boundary and fit output bookkeeping.

These are isolated tests; no actual physics integrations are run here.
"""
from __future__ import annotations

import copy
import math
from unittest import TestCase

import numpy as np

from scripts import RunT3UniformScalarBenchmark as study


class UniformBoundaryTests(TestCase):
    def setUp(self):
        self.cfg = {'base_state':{'sm':{'lambdaH':0.25,'g2':0.59493},
                                 'ordinary':{
            'lambdaT3':{'real':0.1,'imag':0.0},
            'lambdaS1':0.1,'lambdaS2':0.12,'lambdaH1':0.02,
            'lambdaH2':0.03,'lambda12':0.01,
            'lambdaH1Adj':0.0,
            'lambda12Cross':{'real':0.0,'imag':0.0},
            'mS1Sq':1.306122448979592e6,'mS2Sq':1e6,
            'y1':{'real':np.eye(3).tolist(),'imag':np.zeros((3,3)).tolist()},
            'y2':{'real':np.eye(3).tolist(),'imag':np.zeros((3,3)).tolist()}}}}

    def test_all_ordinary_bsm_quartics_001_sm_unchanged(self):
        before=copy.deepcopy(self.cfg)
        result,keys=study.uniform_config(self.cfg,0.0005)
        self.assertEqual(result['base_state']['ordinary']['lambdaT3'],
                         {'real':0.0005,'imag':0.0})
        self.assertEqual(result['base_state']['sm'],before['base_state']['sm'])
        self.assertEqual(result['base_state']['ordinary']['mS1Sq'],
                         before['base_state']['ordinary']['mS1Sq'])
        self.assertEqual(result['base_state']['ordinary']['y1'],
                         before['base_state']['ordinary']['y1'])
        self.assertEqual(set(keys),{'lambdaS1','lambdaS2','lambdaH1',
                                    'lambdaH2','lambda12','lambdaH1Adj','lambda12Cross'})
        for key in keys:
            x=result['base_state']['ordinary'][key]
            self.assertEqual(x if isinstance(x,float) else x['real'],0.01)
            if isinstance(x,dict):
                self.assertEqual(x['imag'],0.0)
        self.assertEqual(self.cfg,before)

    def test_zero_conjugate_yukawa_layout_without_touching_ordinary(self):
        mat=[[{'real':1.0,'imag':0.2} for _ in range(3)] for _ in range(3)]
        fixed={'lambdaH':0.25,'lambdaT3':{'real':0.1,'imag':0},
               'lambdaS2':0.12,'lambdaS1Adj':0.0,
               'lambdaExtra':{'real':0.001,'imag':0.0},
               'y1Conj':mat,'y1':mat}
        before=copy.deepcopy(fixed)
        out,quartics,conjs=study.uniform_overrides(fixed)
        self.assertEqual(quartics,['lambdaExtra','lambdaS1Adj','lambdaS2'])
        self.assertEqual(conjs,['y1Conj'])
        self.assertTrue(all(z=={'real':0.0,'imag':0.0}
                            for row in out['y1Conj'] for z in row))
        self.assertEqual(out['y1'],fixed['y1'])
        self.assertEqual(out['lambdaExtra'],{'real':0.01,'imag':0.0})
        self.assertEqual(out['lambdaH'],0.25)
        self.assertEqual(fixed,before)

    def test_absolute_bc_quartics_independent_of_lambda(self):
        rec={'new':('lambdaNew',),'reused':tuple('lambdaR'+str(i) for i in range(6))}
        extra=study.fixed_bc_extra(rec)
        self.assertEqual(len(extra),7)
        self.assertTrue(all(z==0.01+0j for z in extra.values()))
        for lam in (0.0005,0.005,0.5):
            # No backend quartic fraction has been used or capped here.
            self.assertTrue(math.isclose(extra['lambdaNew'].real,0.01))

    def test_bc_quartic_inventory_failure_loud(self):
        with self.assertRaisesRegex(ValueError,'seven'):
            study.fixed_bc_extra({'new':('lambdaA',),'reused':('lambdaB',)})
        with self.assertRaisesRegex(ValueError,'seven'):
            study.fixed_bc_extra({'new':('lambdaA',),'reused':('lambdaA',)*6})

    def test_scalar_boundary_rejects_invalid_encodings(self):
        with self.assertRaisesRegex(ValueError,'encoding'):
            study.real_scalar_001({'something':0.01})
        with self.assertRaisesRegex(ValueError,'positive'):
            study.uniform_config(self.cfg,0.0)


class SummaryTests(TestCase):
    def test_mixed_accepted_seed_and_refit_are_not_conflated(self):
        rows=[{'model':'B_-1','lambdaT3':lam,'status':status,
               'yukawa_norms':{'combined_frobenius':norm}}
              for lam,status,norm in zip(study.base.BENCHMARKS,
                  ('AcceptedRescaledSeed','AcceptedApproxFit','AcceptedApproxFit'),
                  (0.01,0.003,0.001))]
        out=study.model_summary(rows)
        self.assertEqual(out['status'],'Complete')
        self.assertEqual(out['accepted_count'],3)
        self.assertEqual(out['accepted_rescaled_seeds'],1)
        self.assertIsNotNone(out['log_log_yukawa_norm_vs_lambda_slope_from_accepted_endpoints'])
        rows[1]={'model':'B_-1','lambdaT3':0.005,'status':'Failed'}
        self.assertEqual(study.model_summary(rows)['status'],'Incomplete')


if __name__=='__main__':
    import unittest
    unittest.main()
