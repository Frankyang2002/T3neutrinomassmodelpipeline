from __future__ import annotations
import unittest
import numpy as np
from scripts.RunT3WorkerBEndpointAudit import audit_state, analytic_spectrum

class WorkerBEndpointAuditTests(unittest.TestCase):
    def test_a4_degeneracies(self):
        state={'mS1Sq':10.,'mS2Sq':20.}
        r=audit_state('A_-4',1e3,state)
        self.assertTrue(r['internal_pass'])
        self.assertEqual(r['source_mass_eigenvalues_gev2'],[10.,10.,20.,20.,20.,20.,20.,20.])
    def test_d2_complex_eigenvalues(self):
        state={'mS1Sq':20.,'mS2Sq':10.,'muT3GeneralS2Sq':{'real':3.,'imag':4.}}
        r=audit_state('D_-2',1e3,state)
        self.assertTrue(r['internal_pass'])
        self.assertEqual(r['source_mass_eigenvalues_gev2'],[0.,20.,20.,20.,20.,20.,20.,20.])
    def test_d2_missing_mu_fails(self):
        with self.assertRaises(ValueError):
            analytic_spectrum('D_-2',{'mS1Sq':1.,'mS2Sq':2.})
    def test_negative_masses_reported(self):
        r=audit_state('A_-4',1e3,{'mS1Sq':-1.,'mS2Sq':4.})
        self.assertTrue(r['tachyonic'])
        self.assertTrue(r['internal_pass'])
    def test_historical_reference(self):
        s={'mS1Sq':10.,'mS2Sq':20.}
        r=audit_state('A_-4',1e3,s,reference=[20]*6+[10]*2)
        self.assertEqual(r['historical_max_abs_residual_gev2'],0.)
    def test_bad_historical_shape(self):
        with self.assertRaises(ValueError):
            audit_state('A_-4',1e3,{'mS1Sq':1.,'mS2Sq':2.},reference=np.ones(2))

if __name__=='__main__':unittest.main()
