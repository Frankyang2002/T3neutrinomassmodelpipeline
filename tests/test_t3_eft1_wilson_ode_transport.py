"""Synthetic ODE implementation tests; NOT physical T3 certification."""
from __future__ import annotations

import math
import sys
import types
import unittest
from unittest.mock import patch

import numpy as np

from Numerical.running.T3EFT1WilsonODETransport import (
    LOOP, SourceBridgeError, integrate_master_tensor,
)


class MasterODETests(unittest.TestCase):
    def test_coupled_running_and_zero_boundary_generated_c5(self):
        stub = types.ModuleType('RGE.running.intermediate.GeneralEFT1C5Mixing')

        def fake_master_beta(ctx, C):
            # Seven toy terms: LLSS self-running and LLSS -> HH mixing.
            one = ctx['gamma'] * C
            one = one.copy()
            one[0, 0, 2, 2] += ctx['lambda'] * C[0, 0, 4, 4]
            return {'term'+str(i): one if i == 0 else np.zeros_like(C)
                    for i in range(7)}

        stub.calculate_master_beta = fake_master_beta
        C0 = np.zeros((2, 2, 6, 6), dtype=complex)
        C0[0, 0, 4, 4] = 1e-7 + 2e-8j
        mf, ms = 1e5, 1e3
        grid = []
        def running(mu):
            grid.append(mu)
            return {'gamma': 0.2 + 0.03*math.log(mu/mf),
                    'lambda': 0.04 + 0.002*math.log(mu/mf)}

        with patch.dict(sys.modules, {'RGE.running.intermediate.GeneralEFT1C5Mixing':stub}):
            C, info = integrate_master_tensor(C0, running, mf, ms,
                                               rtol=2e-10, atol=1e-22)
        self.assertGreater(len(grid), 5)
        self.assertEqual(info['one_loop_beta_terms'], 7)
        self.assertNotEqual(C[0, 0, 2, 2], 0)
        self.assertNotEqual(C[0, 0, 4, 4], C0[0, 0, 4, 4])
        # Leading log in scalar-only source: O(1 loop) sign from t_s < t_f.
        leading = LOOP * 0.04 * math.log(ms/mf) * C0[0, 0, 4, 4]
        self.assertLess(abs((C[0, 0, 2, 2] - leading)/leading), 0.3)
        self.assertTrue(np.isfinite(C).all())

    def test_bad_scales_rejected(self):
        with self.assertRaises(SourceBridgeError):
            integrate_master_tensor(np.zeros((2,2,6,6)), lambda _:None,
                                    100, 100)


if __name__=='__main__':
    unittest.main()
