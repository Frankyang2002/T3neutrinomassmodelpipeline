"""Structural tests for the real mixed production entry-point override."""
import unittest
from unittest.mock import patch
import json
import tempfile
from pathlib import Path

import numpy as np
from scripts import RunT3ExploratoryMixedBCWithQuartics as runner
from scripts import RunT3SimultaneousEFT1All16 as general


class FakeOverlay:
    source = 'synthetic-test'
    directions = (('x','real'),('x','imag'))*3
    def evaluate(self, state):
        return np.ones((12,)*4)*state['x']
    def validate_hhss(self, state, hhss):
        return 0.


class TestMixedCLI(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.registry = Path(self.tmp.name)/'registry.json'
        self.registry.write_text(json.dumps({'models': {'B_-1': {
            'new_rgbeta_symbols': [{'name':'x'}], 'preexisting_rgbeta_symbols': []}}}))
    def tearDown(self):
        self.tmp.cleanup()
    def test_override_receives_instantaneous_state(self):
        with patch.object(runner, 'build_overlay', return_value=FakeOverlay()):
            substitute, obj = runner.make_quartic_override(
                'B_-1', lambda *_: None, __file__, self.registry)
            self.assertAlmostEqual(substitute('B_-1', {'x':2}, np.zeros((4,4,8,8)))[0,0,0,0],2.)

    def test_global_override_restored(self):
        previous = general.general_quartic
        with patch.object(runner, 'make_quartic_override',
                          return_value=(lambda *_: 1., FakeOverlay())):
            with runner.full_mixed_quartic('B_-1'):
                self.assertIsNot(general.general_quartic, previous)
            self.assertIs(general.general_quartic, previous)

    def test_wrong_model_fails(self):
        with patch.object(runner, 'build_overlay', return_value=FakeOverlay()):
            substitute, _ = runner.make_quartic_override(
                'B_-1',lambda *_:None,__file__,self.registry)
            with self.assertRaises(RuntimeError):
                substitute('C_-1', {'x':1}, np.zeros((4,4,8,8)))


if __name__ == '__main__':
    unittest.main()
