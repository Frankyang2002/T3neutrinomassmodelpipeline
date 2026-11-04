"""Dependency-light unit tests for the standalone qualitative scalar refitter.

PowerShell: python -X utf8 -m unittest discover -s tests/python -p 'test_t3_qualitative_scalar_refits.py' -v
These tests DO NOT run the physical UV/EFT1/Matchete backends.
"""
from __future__ import annotations

import unittest

import numpy as np

from scripts import RunT3QualitativeScalarRefits as study


class QuarticBoundaryTests(unittest.TestCase):
    def test_bc_lambda_only_preserves_physical_extra_quartics(self):
        q = study.quartic_fraction(0.1, 0.1, 0.5, 'lambda-only')
        self.assertAlmostEqual(q, 0.02)
        self.assertAlmostEqual(q * 0.5, 0.1 * 0.1)

    def test_bc_lambda_only_rejects_both_small_lambdas(self):
        for lam, required in [(0.0005, 20), (0.005, 2)]:
            with self.subTest(lam=lam):
                with self.assertRaisesRegex(study.UnsupportedScalarBoundary, 'required fraction') as cm:
                    study.quartic_fraction(0.1, 0.1, lam, 'lambda-only')
                self.assertIn(f'{required:g}', str(cm.exception))

    def test_bc_joint_changes_extra_quartics_explicitly(self):
        for lam in study.BENCHMARKS:
            q = study.quartic_fraction(0.1, 0.1, lam, 'joint-quartics')
            self.assertAlmostEqual(q, 0.1)
            self.assertAlmostEqual(q * lam, 0.1 * lam)

    def test_nonmixed_extra_only_and_non_mutating(self):
        cfg = {'base_state': {
            'sm': {'lambdaH': 0.25},
            'ordinary': {'lambdaT3': {'real': 0.1, 'imag': 0},
                         'lambdaS1': 0.1, 'lambdaS2Adj': 0.01,
                         'y1': {}, 'y2': {}}}}
        overrides = {
            'lambdaT3': {'real':0.1,'imag':0}, 'lambdaH':0.25,
            'lambdaS1':0.1, 'lambdaS2Adj':0.01,
            'lambdaIndependent': {'real':0.005, 'imag':0.001},
            'lambdaRealExtra':0.004,
            'y1Conj': [[{'real':0.2,'imag':0}]],
        }
        scaled, names = study.scale_extra_quartics(overrides, cfg, 5.0)
        self.assertEqual(names, ['lambdaIndependent','lambdaRealExtra'])
        self.assertAlmostEqual(scaled['lambdaRealExtra'],0.02)
        self.assertEqual(scaled['lambdaIndependent'], {'real':0.025,'imag':0.005})
        self.assertEqual(scaled['lambdaS1'],0.1)
        self.assertEqual(scaled['lambdaT3'],overrides['lambdaT3'])
        self.assertEqual(scaled['y1Conj'],overrides['y1Conj'])
        self.assertEqual(overrides['lambdaIndependent']['real'],0.005)


class FitLogicTests(unittest.TestCase):
    def setUp(self):
        self.config = {'base_state': {'ordinary': {
            'y1': {'real': (np.arange(1,10)*1e-4).reshape(3,3).tolist(),
                   'imag': np.zeros((3,3)).tolist()},
            'y2': {'real': (np.arange(10,19)*1e-4).reshape(3,3).tolist(),
                   'imag': np.zeros((3,3)).tolist()},
            'lambdaT3': {'real':0.1, 'imag':0}}}}

    def test_seed_is_eighteen_real_values_in_y1_y2_order(self):
        x = study.fitted_vector(self.config)
        self.assertEqual(x.shape,(18,))
        self.assertAlmostEqual(x[0], 1e-4)
        self.assertAlmostEqual(x[17], 18e-4)
        self.assertEqual(len(study.parameter_dict(x)), 18)
        self.assertEqual(study.yukawa_matrices(x)['y2']['real'][2][2], x[-1])

    def test_reject_nonzero_imaginary_seed(self):
        self.config['base_state']['ordinary']['y1']['imag'][0][0] = 1e-3
        with self.assertRaisesRegex(ValueError,'18 real'):
            study.fitted_vector(self.config)

    def test_independent_least_squares_invoked_for_all_real_entries(self):
        source = study.fitted_vector(self.config)
        goal = source * 1.5
        seen = []

        class FakeResult:
            ordering = 'NO'
            observable_names = study.OBS
            def __init__(self, x):
                self.whitened_residual = (x[:5] - goal[:5]) / 1e-4
                self.chi2 = float(np.dot(self.whitened_residual,self.whitened_residual))
                self.prediction = np.array([7.537e-5, 0.002511, 0.3088, 0.02248, 0.47])

        def evaluate(params):
            seen.append(tuple(params))
            x = np.array([params[n] for n in study.NAMES], dtype=float)
            return FakeResult(x)

        fit, optimum, record = study.optimize(
            evaluate, source, 0.1, 0.1, bound=0.1,
            max_nfev=40, regularization=1e-4)
        self.assertEqual(len(seen),record['n_evaluator_attempts'])
        self.assertTrue(all(len(s)==18 for s in seen))
        self.assertEqual(optimum.shape,(18,))
        self.assertLess(fit.chi2,1e-3)
        self.assertNotEqual(optimum[0],source[0])
        self.assertTrue(record['all_18_real_yukawas_active'])
        self.assertGreater(record['n_evaluator_attempts'],18)

    def test_summary_handles_only_one_supported_point(self):
        rows = [
            {'model':'B_-1','lambdaT3':0.0005,'status':'UnsupportedScalarBoundary'},
            {'model':'B_-1','lambdaT3':0.005,'status':'UnsupportedScalarBoundary'},
            {'model':'B_-1','lambdaT3':0.5,'status':'AcceptedApproxFit',
             'yukawa_norms':{'combined_frobenius':0.001}},
        ]
        out = study.summarize(rows)
        self.assertEqual(out['accepted_count'],1)
        self.assertEqual(out['status'], 'CompleteWithUnsupported')
        self.assertEqual(out['benchmark_status']['0.5'], 'AcceptedApproxFit')
        self.assertIsNone(out['log_log_yukawa_norm_vs_lambda_slope_from_accepted_endpoints'])

    def test_summary_status_for_complete_and_failed_runs(self):
        rows = [
            {'model':'A_-4', 'lambdaT3':lam, 'status':'AcceptedApproxFit',
             'yukawa_norms':{'combined_frobenius':norm}}
            for lam, norm in zip(study.BENCHMARKS, (0.003,0.001,0.0001))
        ]
        complete = study.summarize(rows)
        self.assertEqual(complete['status'], 'Complete')
        self.assertEqual(complete['accepted_count'], 3)
        self.assertIsNotNone(complete['log_log_yukawa_norm_vs_lambda_slope_from_accepted_endpoints'])
        rows[1] = {'model':'A_-4', 'lambdaT3':0.005, 'status':'Failed'}
        partial = study.summarize(rows)
        self.assertEqual(partial['status'], 'Incomplete')
        self.assertEqual(partial['accepted_count'], 2)


if __name__ == '__main__':
    unittest.main()
