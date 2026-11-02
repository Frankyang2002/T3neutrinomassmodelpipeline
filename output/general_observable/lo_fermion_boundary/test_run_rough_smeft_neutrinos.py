"""Small offline smoke tests for the rough benchmark handoff (no source validation)."""
from __future__ import annotations

import unittest
import numpy as np

import run_rough_smeft_neutrinos as pilot


def enc(a):
    return {'real': np.asarray(a).real.tolist(), 'imag': np.asarray(a).imag.tolist()}


def scalar_fixture():
    lead = np.diag([2.e-15, 3.e-15, 4.e-15])
    fin = np.diag([-1.e-16, -2.e-16, -3.e-16])
    return {'format': 't3-conditional-c5-per-unknown-gammacc-mixed-split-v6',
            'model': 'B_-1',
            'response': {'model': 'B_-1',
                         'status': 'ConditionalOneLoopLLSSResponse_PerUnitUncertifiedSpinorNormalization',
                         'renormalization_scheme': 'MSbar',
                         'mu_f_gev': 100000., 'mu_s_gev': 1000.,
                         'C5_LLSS_total_one_loop_per_N_gev_inverse': enc(lead+fin),
                         'C5_EFT1_induced_leading_log_per_N_gev_inverse': enc(lead),
                         'C5_joint_scalar_finite_bubble_per_N_gev_inverse': enc(fin)}}


class TestRoughPilot(unittest.TestCase):
    def test_total_c5_and_signed_factor(self):
        x, extra = pilot.take_scalar_c5(scalar_fixture(), 'B_-1', -2.)
        self.assertTrue(np.allclose(x, np.diag([-3.8e-15, -5.6e-15, -7.4e-15]), rtol=1e-14, atol=0))
        self.assertEqual(extra['mu_S_gev'], 1000)

    def test_missing_bubble_is_not_double_counted(self):
        s = scalar_fixture()
        s['response']['C5_joint_scalar_finite_bubble_per_N_gev_inverse'] = enc(np.zeros((3,3)))
        with self.assertRaisesRegex(ValueError, 'differs from saved total'):
            pilot.take_scalar_c5(s, 'B_-1', 1.)

    def test_reject_complex_nonsymmetric_c5(self):
        s = scalar_fixture()
        bad = np.diag([1.e-15, 2.e-15, 3.e-15]).astype(complex)
        bad[0,1] = 2.e-15
        s['response']['C5_LLSS_total_one_loop_per_N_gev_inverse'] = enc(bad)
        s['response']['C5_EFT1_induced_leading_log_per_N_gev_inverse'] = enc(bad)
        s['response']['C5_joint_scalar_finite_bubble_per_N_gev_inverse'] = enc(np.zeros((3,3)))
        with self.assertRaisesRegex(ValueError, 'not complex symmetric'):
            pilot.take_scalar_c5(s, 'B_-1', 1.)

    def test_reject_zero_spinor_factor(self):
        with self.assertRaisesRegex(ValueError, 'finite nonzero'):
            pilot.take_scalar_c5(scalar_fixture(), 'B_-1', 0.)

    def test_historical_complex_running_matrix(self):
        source = [[{'real': 1. if a==b else 0., 'imag': 0.} for b in range(3)] for a in range(3)]
        self.assertTrue(np.array_equal(pilot.decoded_state_matrix(source, 'ye'), np.eye(3)))

    def test_c5_json_roundtrip(self):
        z = np.array([[1+2j, 3j, 0], [3j, -4, 0], [0, 0, 5]], complex)
        self.assertTrue(np.array_equal(pilot.decode_matrix(pilot.encode_matrix(z)), z))


if __name__ == '__main__':
    unittest.main(verbosity=2)
