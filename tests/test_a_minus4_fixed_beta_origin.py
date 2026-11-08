"""Synthetic algebra regression for full-beta inclusion/exclusion diagnostic.

These tests check decomposition mathematics, not an actual native beta export.
"""
from __future__ import annotations

from unittest import TestCase, main
from unittest.mock import patch
import numpy as np

from scripts import InspectAminus4FixedFullBetaDifferences as mod


class BetaOriginAlgebra(TestCase):
    def setUp(self):
        self.point = {
            'lambdaH2Adj': 0.2, 'gY': 0.37, 'g2': 0.63, 'g3': 1.1,
            'yu': np.eye(3, dtype=complex) * 0.9,
            'yd': np.eye(3, dtype=complex) * 0.1,
            'ye': np.eye(3, dtype=complex) * 0.03,
        }
        self.fixed = {'label': 'fixed'}
        self.diagnostic = {'label': 'diag'}

    @staticmethod
    def synthetic_beta(beta, p, _):
        x = p['lambdaH2Adj']
        # Same scalar-only beta for fixed and earlier diagnostic.
        base = 12 * x*x
        if beta.get('label') == 'diag':
            return base
        return (base + 7*x*p['g2']**2 + 4*p['gY']**4
                + 3*x*np.trace(p['yu'].conj().T @ p['yu'])
                + 5*p['g2']*np.real(p['yu'][0,0]))

    def _decompose(self):
        with patch.object(mod.eft, '_point', side_effect=lambda p: dict(p)):
            with patch.object(mod,'evaluated',side_effect=self.synthetic_beta):
                return mod.decomposed_difference(self.fixed,self.diagnostic,
                                                   self.point,'lambdaH2Adj')

    def test_four_marginal_probes_reconstruct_full_difference(self):
        d=self._decompose()
        self.assertAlmostEqual(d['scalar_only'].real, 0.0, places=12)
        self.assertAlmostEqual(d['gauge_dependent'].real,
                               7*0.2*0.63**2+4*0.37**4,places=12)
        self.assertAlmostEqual(d['yukawa_dependent'].real,
                               3*0.2*3*0.9**2,places=12)
        self.assertAlmostEqual(d['joint_gauge_yukawa_interference'].real,
                               5*0.63*0.9,places=12)
        self.assertAlmostEqual(d['full'].real,
                               sum(d[k].real for k in ('scalar_only',
                                    'gauge_dependent', 'yukawa_dependent',
                                    'joint_gauge_yukawa_interference')),places=12)

    def test_zero_gauge_and_yukawa_does_not_change_scalar_reference(self):
        with patch.object(mod.eft,'_point',side_effect=lambda p: dict(p)):
            env = mod.input_at(self.point,gauge=False,yukawa=False)
        for name in ('gY','g2','g3'):
            self.assertEqual(env[name],0.0)
        for name in ('yu','yd','ye'):
            self.assertTrue(np.array_equal(env[name],np.zeros((3,3),complex)))
        self.assertEqual(env['lambdaH2Adj'],0.2)
        self.assertEqual(self.point['g2'],0.63)


if __name__ == '__main__':
    main()
