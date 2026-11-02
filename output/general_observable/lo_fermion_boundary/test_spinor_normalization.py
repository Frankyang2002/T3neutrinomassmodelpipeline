"""Dedicated offline matrix/normalization regression (does not require Matchete)."""
from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest

THIS = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('t3_spinor_probe', THIS/'probe_spinor_normalization.py')
module = importlib.util.module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(module)

class SpinorProbeTests(unittest.TestCase):
    def test_exact_clifford_and_charge_conjugation(self):
        self.assertTrue(module.check_exact_algebra()['independent_matrix_identities_pass'])

    def test_conditional_candidate_signs_are_distinct(self):
        r=module.check_exact_algebra()
        self.assertEqual(r['candidate_if_GammaCC_is_C_eq_i_gamma2_gamma0'], '+1')
        self.assertEqual(r['candidate_if_GammaCC_is_C_inverse'], '-1')

    def test_master_orbit_multiplicities(self):
        self.assertEqual(module.check_master_orbit_normalization()['distinct_orbit_sizes_verified'], [1,2,4])

    def test_does_not_certify_physical_factor(self):
        src=(THIS/'probe_spinor_normalization.py').read_text(encoding='utf-8')
        self.assertIn("'certified': False, 'value': None", src)
        self.assertIn("'ready_for_physical_eft1_running': False", src)

if __name__=='__main__':
    unittest.main(verbosity=2)
