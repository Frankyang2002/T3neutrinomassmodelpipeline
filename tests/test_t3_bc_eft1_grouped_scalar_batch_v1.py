import unittest
import numpy as np
from scripts import RunT3BCEFT1GroupedScalarBatch as m

class TestGrouped(unittest.TestCase):
    def fixture(self):
        from scripts import RunT3BCCoupledWilsonSequential as coupled
        c=np.eye(3,dtype=complex)*1e-16
        return {'format':coupled.FORMAT,'status':coupled.STATUS,'model':'B_-1',
           'one_generation_gauge_index_wilson_evolution_completed':True,
           'eft1_samples':65,'wilson_rtol':5e-7,
           'finite_scalar_threshold_matching_included':False,
           'full_general_scalar_matching_certified':False,
           'full_sequential_matching_performed':False,
           'additional_quartic_mapping_physically_certified':False,
           'full_three_generation_wilson_anomalous_dimension_included':False,
           'coupled_C5_at_scalar_threshold_gev_inverse':{'real':c.real.tolist(),'imag':c.imag.tolist()},
           'hard_threshold_grouped_scale_gev':1e5,'scalar_threshold_gev':1e3,
           'takagi_masses_ev_descending':[.03,.02,.003],
           'scalar_threshold_weinberg_residual':1e-9,
           'source_normalization_checks':{'direct_source_parity':1e-8,'stage17_extra_source_parity':1e-8},
           'scalar_ode_nfev':42,'wilson_ode_nfev':53}
    def test_accepts_conditional(self):
        x=m.validate_output('B_-1',self.fixture(),expected_samples=65,expected_rtol=5e-7)
        self.assertFalse(x['physical_C5_certified'])
        self.assertEqual(x['finite_scalar_threshold_correction'],'UNKNOWN_NOT_ZERO')
    def test_rejects_physical_relabel(self):
        f=self.fixture();f['finite_scalar_threshold_matching_included']=True
        with self.assertRaises(ValueError):m.validate_output('B_-1',f,expected_samples=65,expected_rtol=5e-7)
    def test_rejects_wrong_scale(self):
        f=self.fixture();f['scalar_threshold_gev']=2e5
        with self.assertRaises(ValueError):m.validate_output('B_-1',f,expected_samples=65,expected_rtol=5e-7)
    def test_rejects_wrong_sample_count(self):
        with self.assertRaises(ValueError):m.validate_output('B_-1',self.fixture(),expected_samples=33,expected_rtol=5e-7)
    def test_rejects_bad_audit(self):
        with self.assertRaises(ValueError):m.check_boundary({'status':'SomethingElse'})
if __name__=='__main__':unittest.main()
