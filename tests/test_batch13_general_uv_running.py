"""Regression tests for candidate full-inventory UV integration.

Intended path: tests/test_batch13_general_uv_running.py
Run: python -m unittest discover -s tests -p test_batch13_general_uv_running.py -v
"""
import unittest
import numpy as np

from Numerical.experimental.Batch13GeneralUVRunning import (
    EXPECTED_CONVENTION, EXPECTED_FORMAT, EXPECTED_STATUS,
    illustrative_initial_state, run_uv, serialise_state, validate_uv_export,
)


def synthetic_payload(model='A_+0'):
    # A(+0) has dS1=1, dS2=3, dF=2, optional H2/S2 adjoint terms,
    # candidate self-S1 quadratic, and five native extra quartics.
    names = ('gY','g2','g3','yu','yd','ye','lambdaH','y1','y2','MF',
             'mS1Sq','mS2Sq','lambdaS1','lambdaS2','lambdaH1','lambdaH2',
             'lambda12','lambdaT3','lambdaH2Adj','lambdaS2Adj',
             'muT3GeneralS1Sq','lambdaT3GeneralQ000400I1')
    betas = {n: ('0*'+n if n in {'yu','yd','ye','y1','y2','MF'} else '0')
             for n in names}
    return {
        'format': EXPECTED_FORMAT, 'status': EXPECTED_STATUS,
        'physical_mapping_certified': False,
        'metadata': {
            'GeneralScalarModel': model, 'GeneralScalarStage': 'UV',
            'GeneralScalarBasisStatus': 'CandidateNativeRGBetaUnmappedToMatchete',
            'ReportBetaConvention': EXPECTED_CONVENTION,
            'dS1': 1, 'dS2': 3, 'dF': 2, 'alpha': 0,
            'QuadraticKind': 'self_s1',
            'QuadraticMassSymbol': 'muT3GeneralS1Sq',
            'GeneralScalarNewBetaKeys': ['lambdaT3GeneralQ000400I1'],
        },
        'report_betas': betas,
    }


class TestBatch13GeneralUVRunning(unittest.TestCase):
    def test_accepts_exact_batch13_wolfram_export_format(self):
        payload = synthetic_payload()
        self.assertEqual(
            payload['format'],
            'general-t3-quadratic-rgbeta-batch13-restricted-yukawa-v2',
        )
        self.assertEqual(payload['status'], 'SuccessCandidateMassBeta')
        self.assertIn('y1', validate_uv_export(payload, 'A_+0'))
        self.assertIn('y2', validate_uv_export(payload, 'A_+0'))

    def test_rejects_old_export_format_with_diagnostic(self):
        payload = synthetic_payload()
        payload['format'] = 'general-t3-quadratic-rgbeta-experimental-v1'
        with self.assertRaisesRegex(ValueError, 'format=.*experimental-v1'):
            validate_uv_export(payload, 'A_+0')

    def test_inventory_accepts_restricted_yukawa(self):
        payload=synthetic_payload()
        self.assertIn('muT3GeneralS1Sq',validate_uv_export(payload,'A_+0'))

    def test_reject_conjugate_yukawa(self):
        payload=synthetic_payload()
        payload['report_betas']['y1Conj']='0*y1'
        with self.assertRaisesRegex(ValueError, 'conjugate Yukawa'):
            validate_uv_export(payload,'A_+0')

    def test_reject_inconsistent_representation(self):
        payload=synthetic_payload()
        with self.assertRaisesRegex(ValueError, 'metadata'):
            validate_uv_export(payload, 'D_+0')

    def test_initial_state_covers_all_exported_couplings(self):
        payload=synthetic_payload()
        initial=illustrative_initial_state(payload,'A_+0')
        self.assertEqual(set(serialise_state(initial)),set(payload['report_betas']))
        self.assertEqual(initial.general_couplings['muT3GeneralS1Sq'],0j)

    def test_zero_beta_uv_evolution_preserves_all_couplings(self):
        payload=synthetic_payload()
        initial=illustrative_initial_state(payload,'A_+0')
        output=run_uv(payload,'A_+0',initial,mu_final_gev=1e5)
        final=output['state']
        self.assertEqual(output['status'],'CandidateUVIntegrated')
        self.assertFalse(output['matchete_rgbeta_normalisation_certified'])
        self.assertEqual(set(serialise_state(final)),set(payload['report_betas']))
        self.assertEqual(final.mu_gev,1e5)
        self.assertTrue(np.array_equal(final.MF,initial.MF))
        self.assertEqual(final.general_couplings['muT3GeneralS1Sq'],0j)

    def test_missing_quadratic_beta_is_error(self):
        payload=synthetic_payload()
        del payload['report_betas']['muT3GeneralS1Sq']
        with self.assertRaisesRegex(ValueError,'quadratic beta missing'):
            validate_uv_export(payload,'A_+0')


if __name__=='__main__': unittest.main()
