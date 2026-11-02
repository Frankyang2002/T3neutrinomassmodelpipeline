"""Offline non-certification and schema checks for native Feynman rule evidence."""
from __future__ import annotations
import unittest
from native_feynman_rules_check import validate_report


def dummy():
    obj = {
        'format': 't3-native-matchete-fr-evidence-v1',
        'model': 'A_-4',
        'status': 'NativeUV3pointAndEFT1ContactVerticesExtracted_AmplitudeOpen',
        'matchete_gamma_cc_to_weyl_phase_certified': False,
        'independent_connected_uv_four_point_amplitude_computed': False,
        'physical_master_wilson_normalization_certified': False,
        'original_16_model_release_changed': False,
        'uv_y1_with_hc': 2, 'uv_y2_with_hc': 2,
        'uv_bsm_yukawa_terms': 4,
        'native_source_index': 2,
        'native_source_inputform': 'abc',
        'uv_three_point_rules_extracted': True,
        'eft1_four_point_rules_extracted': True,
        'heavy_propagator_and_both_exchange_orderings_contracted': False,
    }
    for f in ('uv_yukawa_lagrangian', 'uv_three_point_rules',
              'eft1_four_point_rules', 'uv_mass_lagrangian_source'):
        obj[f] = {'inputform_characters': 3,
                  'inputform_complete': True, 'inputform': 'abc',
                  'inputform_sha256': 'a' * 64}
    return obj


class NativeFeynmanTests(unittest.TestCase):
    def test_schema_checks(self):
        self.assertIn('NoAmplitudeCertification', validate_report(dummy(), 'abc'))

    def test_source_mismatch_blocked(self):
        with self.assertRaises(AssertionError):
            validate_report(dummy(), 'different native source')

    def test_amplitude_claim_blocked(self):
        obj = dummy()
        obj['independent_connected_uv_four_point_amplitude_computed'] = True
        with self.assertRaises(AssertionError):
            validate_report(obj, 'abc')

    def test_phase_claim_blocked(self):
        obj = dummy()
        obj['matchete_gamma_cc_to_weyl_phase_certified'] = True
        with self.assertRaises(AssertionError):
            validate_report(obj, 'abc')

    def test_overflow_not_mislabeled(self):
        obj = dummy()
        obj['uv_three_point_rules']['inputform_complete'] = False
        with self.assertRaises(AssertionError):
            validate_report(obj, 'abc')

    def test_blocked_diagnostic_retained(self):
        obj = dummy()
        obj['status'] = 'Blocked'
        obj['blocked_phase'] = 'UV Feynman rules'
        obj['reason'] = 'an error'
        self.assertEqual(validate_report(obj, 'abc'), 'Blocked')


if __name__ == '__main__':
    unittest.main(verbosity=2)
