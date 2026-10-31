"""Read-only checks of existing or freshly reproduced B/C EFT2 conditional reports.

Run the ORIGINAL legacy runner with its ORIGINAL inputs and a NEW output-id
before this test if you need fresh calculation parity. Never regenerate or
replace a missing report from a guessed C5 or matching coefficient.
"""
from pathlib import Path
import json
import os
import unittest

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
REPORT_DIR=Path(os.environ.get('T3_BC_EFT2_REPORT_DIR', str(
    ROOT/'output/studies/general_scalar_extended/continuous_full_scalar'/
    'continuous_v1/bc_eft2_sequential_v1')))
# Values retrieved read-only from the public main-branch saved original reports.
EXPECTED={
    'B_-1':(0.03853362292521312,0.01561382145082619,0.00045848474457153946),
    'C_-1':(0.046405741900074446,0.045923847867146636,0.021684644219742994),
}


class TestExistingConditionalBCEFT2Reports(unittest.TestCase):
    def _check(self,model):
        file=REPORT_DIR/(model.replace('+','p').replace('-','m')+'_eft2_sequential.json')
        if not file.is_file():
            self.skipTest('Actual B/C saved report unavailable; not reproduced: '+str(file))
        report=json.loads(file.read_text(encoding='utf-8-sig'))
        self.assertEqual(report['model'],model)
        self.assertEqual(report['format'],'t3-bc-conditional-sequential-hx-to-sm-v1')
        self.assertEqual(report['status'],
                         'EFT2HXRunningXCrossingSMWeinbergConditionalFiniteScalarMatchingMissing')
        self.assertIs(report['physical_neutrino_mass_prediction_certified'],False)
        self.assertIs(report['finite_Y_threshold_delta_C_computed'],False)
        self.assertIs(report['finite_X_threshold_delta_C5_computed'],False)
        self.assertEqual(report['finite_threshold_corrections_status'],
                         'BOTH_MISSING_NOT_PHYSICAL_ZERO')
        self.assertTrue(report['mu_Y_gev']>report['mu_X_corrected_gev']>100.)
        masses=np.asarray(report['takagi_masses_ev_descending'],dtype=float)
        self.assertEqual(masses.shape,(3,))
        self.assertTrue(np.isfinite(masses).all())
        np.testing.assert_allclose(masses,np.array(EXPECTED[model]),rtol=0,atol=1e-10)

    def test_B_minus_1(self):
        self._check('B_-1')

    def test_C_minus_1(self):
        self._check('C_-1')


if __name__=='__main__':
    unittest.main()
