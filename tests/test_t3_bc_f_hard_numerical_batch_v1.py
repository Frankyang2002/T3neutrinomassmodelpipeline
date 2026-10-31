import sys
import types
import unittest
from unittest.mock import patch
import numpy as np

from Numerical.experimental import T3BCFHardNumericalBoundary as core

Q={'lambdaT3':1.2+0.4j,'lambdaT3ExtraQ020020I1':0.3,'lambdaT3ExtraQ020200I1':-0.2}
U=np.eye(2,dtype=complex)
M=np.array([2e4,3e4,5e4])
A=np.array([[0.1,0.2,0.3],[0.04,0.02,0.07],[0.03,0.01,0.02]],complex)
B=np.array([[0.3,0.2,0.5],[0.1,0.2,0.4],[0.1,0.2,0.4]],complex)

def mon(ch,scale):
    return {'chirality':ch,'quartic':'lambdaT3',
       'quartic_conjugated':ch=='PR',
       'yukawa_i':'y1','yukawa_j':'y1',
       'yukawas_conjugated':ch=='PL',
       'mixing':[{'name':'T3MEU11','conjugated':ch=='PR','exponent':4}],
       'coefficient_numerator':scale,'coefficient_denominator':1,'sqrt3_exponent':0}

def ledger(scale):
    return {'model':'B_-1','monomials':[mon('PL',scale),mon('PR',scale)]}

def evaluator(record,*,U,quartics,y1,y2,MF,hbar):
    result={'PL':np.zeros((3,3),complex),'PR':np.zeros((3,3),complex)}
    for x in record['monomials']:
        yy=np.conj(y1) if x['yukawas_conjugated'] else y1
        coupling=quartics[x['quartic']]
        if x['quartic_conjugated']:
            coupling=np.conj(coupling)
        result[x['chirality']]+=x['coefficient_numerator']*coupling*hbar/MF*np.outer(yy,yy)
    return result

class TestNativeFormalNumerics(unittest.TestCase):
    def setUp(self):
        fake=types.ModuleType('Numerical.experimental.T3BCNativeYLLHHFlavorKernel')
        fake.evaluate_spurion_kernel=evaluator
        self.patcher=patch.dict(sys.modules,{'Numerical.experimental.T3BCNativeYLLHHFlavorKernel':fake})
        self.patcher.start()
    def tearDown(self):self.patcher.stop()
    def args(self):return dict(masses=M,y1=A,y2=B,U=U,quartics=Q,mu_f=1e5)
    def test_3_heavy_masses_and_warsaw_symmetric(self):
        c,parts,hc=core.evaluate_sector(ledger(2),U=U,quartics=Q,masses=M,y1=A,y2=B)
        self.assertEqual(parts.shape,(3,3,3))
        self.assertLess(np.linalg.norm(c-c.T),1e-24)
        self.assertLess(hc,1e-14)
        self.assertGreater(np.linalg.norm(c),0)
    def test_pole_log_and_log_derivative(self):
        raw=core.evaluate_three({'pole':ledger(2),'log':ledger(2),'constant':ledger(4)},**self.args())
        self.assertLess(raw['pole_log_matrix_relative_residual'],1e-14)
        self.assertLess(raw['finite_difference_log_derivative_relative_residual'],2e-9)
        self.assertFalse(raw['renormalized_boundary'])
        self.assertIsNone(raw['physical_C5'])
    def test_individual_MF_log(self):
        raw=core.evaluate_three({'pole':ledger(2),'log':ledger(2),'constant':ledger(4)},**self.args())
        self.assertTrue(np.allclose(raw['matching_scale_log_per_heavy'],2*np.log(1e5/M)))
    def test_log_pole_tamper_rejected(self):
        with self.assertRaisesRegex(ValueError,'pole/log'):
            core.evaluate_three({'pole':ledger(2),'log':ledger(3),'constant':ledger(4)},**self.args())
    def test_scale_scan(self):
        x,scan=core.scale_scan({'pole':ledger(2),'log':ledger(2),'constant':ledger(4)},**self.args())
        self.assertEqual(len(scan),3)
        self.assertTrue(np.allclose(scan[1]['matrix'],x['formal_at_mu_f']))
        self.assertLess(scan[0]['frobenius_norm'],scan[2]['frobenius_norm'])
    def test_bad_masses(self):
        with self.assertRaisesRegex(ValueError,'masses'):
            core.normalise_inputs([0,3,4],A,B,U,Q)
    def test_bad_yukawa(self):
        with self.assertRaisesRegex(ValueError,'Yukawas'):
            core.normalise_inputs(M,A[:2],B,U,Q)
    def test_bad_unitary(self):
        with self.assertRaisesRegex(ValueError,'rotation'):
            core.normalise_inputs(M,A,B,2*U,Q)
    def test_unexpected_quartics(self):
        with self.assertRaisesRegex(ValueError,'quartics'):
            core.normalise_inputs(M,A,B,U,{'lambdaT3':1})
    def test_incomplete_sector(self):
        with self.assertRaisesRegex(ValueError,'incomplete'):
            core.evaluate_three({'pole':ledger(2),'log':ledger(2)},**self.args())
    def test_missing_numeric_physical_matching(self):
        import inspect
        from pathlib import Path
        source=Path(__file__).resolve().parents[1]/'scripts/RunT3BCFHardNumericalBatch.py'
        text=source.read_text()
        self.assertIn("'independent_UV_EFT1_counterterms_computed':False",text)
        self.assertIn("'physical_C5':None",text)
        self.assertNotIn('EFT2_evolution',text)

if __name__=='__main__':unittest.main()
