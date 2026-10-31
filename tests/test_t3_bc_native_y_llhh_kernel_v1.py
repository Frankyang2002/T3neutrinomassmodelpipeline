"""Regression and negative tests for the native finite Y LLHH component reduction."""
from __future__ import annotations
from fractions import Fraction
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np
from Numerical.experimental import T3BCNativeYLLHHFlavorKernel as k
from scripts.ProjectT3BCNativeYLLHH import (
    FINITE, POLE, ROOT, run, FORMAT, STATUS,
)


class TestComponentLogic(unittest.TestCase):
    def test_epsilon_orientation(self):
        self.assertEqual(k._epsilon_value(1,2),1)
        self.assertEqual(k._epsilon_value(2,1),-1)
        self.assertEqual(k._epsilon_value(1,1),0)

    def test_epsilon_out_of_range_rejected(self):
        with self.assertRaises(k.KernelError):k._epsilon_value(3,1)

    def test_unrecognized_source_fails(self):
        with self.assertRaises(k.KernelError):k.reduce_finite_y_llhh('hbar','B_-1')

    def test_wrong_model_fails(self):
        with self.assertRaises(k.KernelError):k._record_from_term('0','B_0')

    def test_valid_conjugate_pair(self):
        a={'chirality':'PR','quartic':'lambdaT3','quartic_conjugated':True,
           'yukawa_i':'y1','yukawa_j':'y2','yukawas_conjugated':False,
           'mixing':[{'name':'T3MEU11','conjugated':False,'exponent':2}],
           'sqrt3_exponent':0,'coefficient_numerator':2,'coefficient_denominator':3}
        b={**a,'chirality':'PL','quartic_conjugated':False,
           'yukawas_conjugated':True,
           'mixing':[{'name':'T3MEU11','conjugated':True,'exponent':2}]}
        self.assertTrue(k.check_exact_chirality_conjugacy({'monomials':[a,b]}))

    def test_wrong_relative_sign_rejected(self):
        a={'chirality':'PR','quartic':'lambdaT3','quartic_conjugated':True,
           'yukawa_i':'y1','yukawa_j':'y1','yukawas_conjugated':False,
           'mixing':[{'name':'T3MEU11','conjugated':False,'exponent':4}],
           'sqrt3_exponent':0,'coefficient_numerator':2,'coefficient_denominator':1}
        b={**a,'chirality':'PL','quartic_conjugated':False,
           'yukawas_conjugated':True,'coefficient_numerator':-2,
           'mixing':[{'name':'T3MEU11','conjugated':True,'exponent':4}]}
        with self.assertRaises(k.KernelError):k.check_exact_chirality_conjugacy({'monomials':[a,b]})

    def test_numeric_evaluation_preserves_chiralities(self):
        a={'chirality':'PR','quartic':'lambdaT3','quartic_conjugated':False,
           'yukawa_i':'y1','yukawa_j':'y2','yukawas_conjugated':False,
           'mixing':[{'name':'T3MEU11','conjugated':False,'exponent':4}],
           'sqrt3_exponent':0,'coefficient_numerator':2,'coefficient_denominator':1}
        b={**a,'chirality':'PL','quartic_conjugated':True,
           'yukawas_conjugated':True,
           'mixing':[{'name':'T3MEU11','conjugated':True,'exponent':4}]}
        d={'model':'B_-1','monomials':[a,b]}
        y1=np.array([1+1j,2,0],complex)
        y2=np.array([0,3-2j,4],complex)
        out=k.evaluate_spurion_kernel(d,U=np.eye(2),y1=y1,y2=y2,
                                      quartics={'lambdaT3':0.5+0.7j},MF=2,hbar=0.25)
        expected=(2*(.5+.7j)*(.25/2))*np.outer(y1,y2)
        np.testing.assert_allclose(out['PR'],expected,atol=1e-14)
        np.testing.assert_allclose(out['PL'],expected.conj(),atol=1e-14)

    def test_invalid_U_rejected(self):
        with self.assertRaises(k.KernelError):
            k.evaluate_spurion_kernel({'model':'B_-1','monomials':[{'quartic':'lambdaT3'}]},
               U=np.ones((2,2)),y1=np.ones(3),y2=np.ones(3),
               quartics={'lambdaT3':1},MF=1,hbar=1)


class TestSourceAndDriver(unittest.TestCase):
    def test_certification_false(self):
        p=ROOT/'Numerical/experimental/T3BCNativeYLLHHFlavorKernel.py'
        source=p.read_text(encoding='utf8')
        self.assertIn("'native_Y_finite_C5_physical_certified':False",source)
        self.assertIn('def check_exact_chirality_conjugacy',source)

    def test_no_untrusted_expression_execution(self):
        source=(ROOT/'Numerical/experimental/T3BCNativeYLLHHFlavorKernel.py').read_text(encoding='utf8')
        self.assertNotIn('eval(',source)
        self.assertNotIn('sympify(',source)

    def test_paths_in_project(self):
        self.assertIn('output/audits',str(FINITE))
        self.assertIn('output/audits',str(POLE))
        self.assertTrue(FORMAT.startswith('t3-bc-'))
        self.assertIn('NoPhysical',STATUS)

    def test_real_archived_evidence_if_available(self):
        f=ROOT/FINITE/'summary.json';p=ROOT/POLE/'summary.json'
        if not f.is_file() or not p.is_file():
            self.skipTest('Prior v1/v2 audits are not mounted')
        with TemporaryDirectory() as td:
            result=run(('B_-1','C_-1'),finite_dir=f.parent,pole_dir=p.parent,
                       out_dir=Path(td),root=ROOT)
            self.assertFalse(result['physical_sequential_matching_certified'])
            for m in ('B_m1','C_m1'):
                r=json.loads((Path(td)/(m+'_native_y_llhh_neutral_flavor_kernel.json')).read_text())
                self.assertEqual(r['source_additive_terms'],32)
                self.assertEqual(r['chirality_monomial_counts'],{'PR':16,'PL':16})
                self.assertTrue(k.check_exact_chirality_conjugacy(r))


if __name__=='__main__':unittest.main()
