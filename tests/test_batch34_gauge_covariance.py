"""Fail-closed expected-failure Ward regression for synthetic SM toy."""
import unittest
from scripts.CheckBatch34GaugeCovariance import diagnose, TOL
class WardGateTest(unittest.TestCase):
    def test_nonzero_yukawa_quartic_and_weinberg_inputs_are_invariant(self):
        rep=diagnose()
        for val in rep['input_max_absolute_ward_residuals'].values():
            self.assertLess(val,TOL)
    def test_both_crossed_conventions_fail_gauge_ward(self):
        rep=diagnose()
        self.assertEqual(rep['status'],'BLOCKED_GAUGE_WARD_VIOLATION')
        self.assertEqual(len(rep['failed_variants']),2)
        for val in rep['variants'].values():
            self.assertGreater(val['full_beta_relative_ward_residual'],0.1)
            self.assertGreater(val['per_term_absolute_ward_residual']['crossed_yukawa'],1e-3)
    def test_diagnostic_is_not_certification(self):
        self.assertFalse(diagnose()['physical_Aminus4_certified'])
if __name__=='__main__':unittest.main()
