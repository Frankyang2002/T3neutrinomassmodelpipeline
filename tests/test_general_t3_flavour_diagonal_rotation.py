"""Regression tests for the positive-diagonal vectorlike MF basis fast path.

Project path: tests/test_general_t3_flavour_diagonal_rotation.py
Run with python -X utf8 -m unittest discover -s tests -p test_general_t3_flavour_diagonal_rotation.py
"""
from types import SimpleNamespace
import unittest

import numpy as np

from Numerical.matching.GeneralT3FlavourC5 import (
    _diagonal_positive_identity_right_rotation,
    prepare_heavy_flavour,
)
from Numerical.running.FinalC5TrajectoryAdapter import _heavy_mass_basis_data


class VectorlikeDiagonalRotationTests(unittest.TestCase):
    @staticmethod
    def state(mf=None):
        z = np.diag([1.5e5, 1.8e5, 2.0e5]).astype(complex) if mf is None else mf
        y1 = np.array([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6], [0.7, 0.8, 0.9]], dtype=complex)
        y2 = np.array([[0.2j, 0.1, 0.3], [0.3, 0.4, 0.5], [0.1, 0.2, 0.3]], dtype=complex)
        conj = np.array([[0.11, 0.22, 0.33], [0.44, 0.55, 0.66], [0.77, 0.88, 0.99]], dtype=complex)
        return SimpleNamespace(
            MF=z, y1=y1, y2=y2,
            general_couplings={'y2Conj': conj},
            representation=SimpleNamespace(majorana_fermion=False),
        )

    def test_real_adapter_positive_diagonal_has_no_right_rotation(self):
        state = self.state()
        basis = _heavy_mass_basis_data(
            state, offdiagonal_atol=1e-7, imaginary_atol=1e-9,
            takagi_residual_rtol=1e-8, takagi_residual_atol=1e-7,
        )
        self.assertIsNone(basis.right_rotation)
        np.testing.assert_allclose(basis.rotation, np.eye(3), atol=1e-12)
        np.testing.assert_allclose(basis.masses_gev, [150000, 180000, 200000])

    def test_extra_vectorlike_y2conj_uses_identity_in_diagonal_basis(self):
        state = self.state()
        masses, yukawas, diagnostics = prepare_heavy_flavour(
            state, {'y1', 'y2', 'y2Conj'})
        np.testing.assert_allclose(masses, [150000, 180000, 200000])
        np.testing.assert_allclose(yukawas['y1'], state.y1)
        np.testing.assert_allclose(yukawas['y2'], state.y2)
        np.testing.assert_allclose(yukawas['y2Conj'], state.general_couplings['y2Conj'])
        self.assertTrue(diagnostics['vectorlike_right_rotation_identity_fast_path'])
        self.assertEqual(diagnostics['extra_yukawa_rotation'],
                         'same-heavy-leg-as-parent-assumption')

    def test_non_diagonal_mass_uses_actual_right_svd_matrix(self):
        mf = np.array([[1.5e5, 1e4, 0], [0, 1.8e5, 4e3], [0, 0, 2e5]], dtype=complex)
        state = self.state(mf)
        masses, yukawas, diag = prepare_heavy_flavour(state, {'y1', 'y2', 'y2Conj'})
        basis = _heavy_mass_basis_data(
            state, offdiagonal_atol=1e-7, imaginary_atol=1e-9,
            takagi_residual_rtol=1e-8, takagi_residual_atol=1e-7,
        )
        self.assertIsNotNone(basis.right_rotation)
        np.testing.assert_allclose(yukawas['y2Conj'],
                                   state.general_couplings['y2Conj'] @ basis.right_rotation,
                                   rtol=1e-12, atol=1e-12)
        self.assertFalse(diag['vectorlike_right_rotation_identity_fast_path'])
        self.assertTrue(np.all(masses > 0))

    def test_missing_nontrivial_right_rotation_fails_closed(self):
        state = self.state(np.array([[150000, 1000, 0], [0, 180000, 0], [0, 0, 200000]], dtype=complex))
        basis = SimpleNamespace(rotation=np.eye(3), masses_gev=np.array([150000, 180000, 200000]))
        with self.assertRaisesRegex(ValueError, 'non-diagonal MF'):
            _diagonal_positive_identity_right_rotation(state, basis)

    def test_complex_or_nonpositive_diagonal_fails_closed(self):
        basis = SimpleNamespace(rotation=np.eye(3), masses_gev=np.array([150000, 180000, 200000]))
        complex_mf = np.diag([150000+1j, 180000, 200000])
        with self.assertRaisesRegex(ValueError, 'non-positive-real MF'):
            _diagonal_positive_identity_right_rotation(self.state(complex_mf), basis)
        negative_mf = np.diag([-150000, 180000, 200000])
        with self.assertRaisesRegex(ValueError, 'non-positive-real MF'):
            _diagonal_positive_identity_right_rotation(self.state(negative_mf), basis)


if __name__ == '__main__':
    unittest.main()
