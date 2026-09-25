from __future__ import annotations

import unittest

import numpy as np
import sympy as sp

from Numerical.WeinbergRunning import neutrino_mass_matrix
from physics.NeutrinoMass import build_neutrino_mass_matrix


class NeutrinoMassConventionTests(unittest.TestCase):
    """Keep symbolic and numerical C5 -> m_nu conventions identical."""

    def setUp(self) -> None:
        self.vev_gev = 246.22
        self.c5 = np.array(
            [
                [1.0e-14 + 2.0e-15j, 3.0e-15 - 1.0e-15j, -2.0e-15 + 4.0e-16j],
                [3.0e-15 - 1.0e-15j, 8.0e-15 + 0.0j, 5.0e-16 + 7.0e-16j],
                [-2.0e-15 + 4.0e-16j, 5.0e-16 + 7.0e-16j, 1.2e-14 - 3.0e-15j],
            ],
            dtype=complex,
        )

    def test_numeric_conversion_is_minus_v_squared_over_two_c5(self) -> None:
        expected = -0.5 * self.vev_gev**2 * self.c5
        result = neutrino_mass_matrix(self.c5, vev_gev=self.vev_gev)
        np.testing.assert_allclose(result, expected, rtol=1.0e-14, atol=0.0)

    def test_symbolic_and_numeric_conversions_agree(self) -> None:
        c5_symbolic = sp.Matrix(self.c5.tolist())
        symbolic_mass = build_neutrino_mass_matrix(
            c5_symbolic,
            vev=sp.Float(self.vev_gev),
        )
        symbolic_numeric = np.asarray(
            symbolic_mass.evalf().tolist(),
            dtype=complex,
        )
        numerical_mass = neutrino_mass_matrix(
            self.c5,
            vev_gev=self.vev_gev,
        )
        np.testing.assert_allclose(
            numerical_mass,
            symbolic_numeric,
            rtol=1.0e-13,
            atol=1.0e-30,
        )

    def test_numeric_conversion_preserves_symmetry(self) -> None:
        result = neutrino_mass_matrix(self.c5, vev_gev=self.vev_gev)
        np.testing.assert_allclose(
            result,
            result.T,
            rtol=1.0e-14,
            atol=0.0,
        )


if __name__ == "__main__":
    unittest.main()
