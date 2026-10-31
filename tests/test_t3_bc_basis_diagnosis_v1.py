import unittest
from Numerical.experimental.T3BCBasisDiagnosis import analyse_report


def fixture(mf, scalar):
    d = {'MF': mf, 'y1': 2e-8, 'y2': 4e-8, 'mS1Sq': scalar}
    z = {key: 2e-16 for key in d}
    return {'model': 'B_-1', 'status': 'SourceStudy',
            'native_RGBeta_heavy_basis_stress': {
                'status': 'NativeRGBetaComplexHeavyBasisCovarianceComputed',
                'beta_tests': {
                    'source_CConj_orientation': {'beta_relative_residual_by_parameter': d},
                    'same_U_negative_control': {'beta_relative_residual_by_parameter': z},
                },
            }}


class DiagnosisTests(unittest.TestCase):
    def test_expected_failed_blocks(self):
        r = analyse_report(fixture(0.78, 0.001))
        self.assertEqual(set(r['CConj_orientation_failed_blocks_at_1e-6']), {'MF', 'mS1Sq'})
        self.assertTrue(r['same_U_passes_at_1e-6'])
        self.assertFalse(r['physical_C5_certified'])

    def test_no_false_failure_for_tiny_blocks(self):
        r = analyse_report(fixture(1e-8, 1e-9))
        self.assertTrue(r['CConj_orientation_passes_at_1e-6'])

    def test_missing_data_rejected(self):
        with self.assertRaises(ValueError):
            analyse_report({'model': 'B_-1'})


if __name__ == '__main__':
    unittest.main()
