"""Unit tests for running scalar quartic response reconstruction."""
import unittest
import numpy as np
from Numerical.running.T3RunningQuarticTensor import (
    RunningQuarticBasis, QuarticEvidenceError
)


def payload(model='A_-4'):
    n = 12
    t = np.zeros((n,)*4)
    t[4,4,4,4] = 24.
    h = np.zeros((n,)*4)
    for a in range(4):
        for b in range(4, n):
            for perm in set(__import__('itertools').permutations((a,a,b,b))):
                h[perm] = 2.
    return dict(format='t3-running-real-quartic-response-v1',
                model=model,
                convention='V4=lambda_abcd phi_a phi_b phi_c phi_d/4!',
                basis='fixed-original-Cartesian',
                source='synthetic-test-only',
                constant=np.zeros((n,)*4).tolist(),
                responses=[
                    dict(symbol='lambdaS1',component='real',tensor=t.tolist()),
                    dict(symbol='lambdaH1',component='real',tensor=h.tolist()),
                ])


class TestRunningQuartic(unittest.TestCase):
    def test_linear_running_and_hhss(self):
        model = RunningQuarticBasis.from_payload(
            payload(), {'lambdaS1','lambdaH1'})
        one = model.evaluate({'lambdaS1':0.3,'lambdaH1':0.2})
        two = model.evaluate({'lambdaS1':0.9,'lambdaH1':0.4})
        self.assertAlmostEqual(one[4,4,4,4], 7.2)
        self.assertAlmostEqual(two[4,4,4,4], 21.6)
        self.assertAlmostEqual(one[0,0,4,4], 0.4)
        self.assertAlmostEqual(two[0,0,4,4], 0.8)
        self.assertAlmostEqual(model.validate_hhss(
            {'lambdaS1':0.3,'lambdaH1':0.2},one[:4,:4,4:,4:]),0.)

    def test_missing_direction_fails(self):
        x = payload()
        x['responses'].pop()
        with self.assertRaises(QuarticEvidenceError):
            RunningQuarticBasis.from_payload(x, {'lambdaS1','lambdaH1'})

    def test_wrong_convention_fails(self):
        x = payload()
        x['convention'] = 'unknown'
        with self.assertRaises(QuarticEvidenceError):
            RunningQuarticBasis.from_payload(x, {'lambdaS1','lambdaH1'})

    def test_complex_imaginary_direction(self):
        x = payload()
        x['responses'][0]['component'] = 'imag'
        x['responses'].append(dict(
            symbol='lambdaS1', component='real',
            tensor=np.zeros((12,)*4).tolist()))
        b = RunningQuarticBasis.from_payload(x, {'lambdaS1','lambdaH1'})
        out = b.evaluate({'lambdaS1':{'real':0.2,'imag':0.6},
                          'lambdaH1':0.})
        self.assertAlmostEqual(out[4,4,4,4],14.4)

    def test_bad_basis_fails(self):
        x = payload()
        x['basis'] = 'moving-mass-basis'
        with self.assertRaises(QuarticEvidenceError):
            RunningQuarticBasis.from_payload(x, {'lambdaS1','lambdaH1'})


if __name__ == '__main__':
    unittest.main()
