"""Regression checks for fixed-order B/C complex matrix input encoding."""
from __future__ import annotations
import importlib.util
from pathlib import Path
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('bc_fixed_order', ROOT/'scripts/RunT3BCFixedOrderMasses.py')
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)

class TestComplexInput(unittest.TestCase):
    def test_archived_real_imag_rotation(self):
        z = {'real': [[1., 0.], [0., 0.]], 'imag': [[0., 0.], [0., 1.]]}
        np.testing.assert_allclose(MODULE.matrix(z, (2,2)), np.diag([1, 1j]))
    def test_complex_scalar(self):
        self.assertEqual(MODULE.as_complex(1+2j), 1+2j)
    def test_complex_pairs(self):
        z = [[[1, 2], [3, 4]], [[5, 6], [7, 8]]]
        np.testing.assert_allclose(MODULE.matrix(z, (2,2)), [[1+2j,3+4j],[5+6j,7+8j]])
    def test_bad_shape(self):
        with self.assertRaises(ValueError): MODULE.matrix([[1]], (2,2))

if __name__ == '__main__': unittest.main()
