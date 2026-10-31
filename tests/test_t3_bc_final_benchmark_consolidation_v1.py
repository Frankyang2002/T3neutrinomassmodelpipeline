import unittest
import numpy as np
from Numerical.experimental import T3BCFinalBenchmarkConsolidation as c

class TestChecks(unittest.TestCase):
    def example(self):
        row={'model':'B_-1','physical_C5_certified':False,'full_matching_certified':False,
          'finite_scalar_threshold_correction':'UNKNOWN_NOT_ZERO','mu_F_gev':1e5,
          'mu_S_grouped_gev':1e3,'candidate_masses_ev_descending':[.04,.02,.003]}
        return {'format':'t3-bc-single-eft1-grouped-scalar-conditional-v1',
          'status':'HistoricalEFT1CoupledWilsonGroupedScalarCandidateNoFiniteMatching',
          'same_trajectory_as_batch1_five_mass_running_verified':False,
          'batch1_boundary_numerically_inserted':False,
          'one_intermediate_EFT_only':True,'finite_scalar_threshold_matching_computed':False,
          'physical_C5_certified':False,'models':{'B_-1':row}}
    def test_historical_acceptance(self):
        r=c.validate_grouped('B_-1',self.example())
        self.assertFalse(r['same_five_mass_trajectory'])
        self.assertEqual(r['mu_S_gev'],1000.)
    def test_forged_same_trajectory_rejected(self):
        d=self.example();d['same_trajectory_as_batch1_five_mass_running_verified']=True
        with self.assertRaises(ValueError):c.validate_grouped('B_-1',d)
    def test_forged_physical_matching_rejected(self):
        d=self.example();d['models']['B_-1']['full_matching_certified']=True
        with self.assertRaises(ValueError):c.validate_grouped('B_-1',d)
    def test_spectrum_checks(self):
        np.testing.assert_array_equal(c.sorted_masses([.03,.02,.001],'a'),[.03,.02,.001])
        with self.assertRaises(ValueError):c.sorted_masses([.01,.03,.02],'bad')
    def test_cannot_infer_finite_scalar_matching(self):
        d=self.example();d['finite_scalar_threshold_matching_computed']=True
        with self.assertRaises(ValueError):c.validate_grouped('B_-1',d)

if __name__=='__main__':unittest.main()
