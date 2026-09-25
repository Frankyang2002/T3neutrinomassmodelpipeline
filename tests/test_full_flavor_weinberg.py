from __future__ import annotations

import unittest
import numpy as np

from Numerical.WeinbergRunning import (
    LOOP,
    SMInitialConditions,
    _beta,
    _pack,
    _unpack,
    neutrino_mass_matrix,
)
from Numerical.FinalSMBoundary import (
    FinalSMBoundaryState,
    build_weinberg_initial_conditions,
)
from Numerical.State import SMNumericalState
from Numerical.WeinbergTrajectory import charged_lepton_mass_basis_matrix


class FullFlavorWeinbergTests(unittest.TestCase):
    def _initial(self) -> SMInitialConditions:
        ye = np.array(
            [
                [3.0e-6, 1.0e-5, 0.0],
                [2.0e-5, 6.0e-4, 3.0e-5],
                [0.0, 4.0e-5, 1.0e-2],
            ],
            dtype=complex,
        )
        yu = np.diag([1.0e-5, 7.0e-3, 0.75]).astype(complex)
        yd = np.diag([2.0e-5, 4.0e-4, 1.8e-2]).astype(complex)
        K = np.array(
            [
                [1.0, 0.2 + 0.1j, 0.1],
                [0.2 + 0.1j, 1.5, -0.3j],
                [0.1, -0.3j, 2.0],
            ],
            dtype=complex,
        ) * 1.0e-14
        return SMInitialConditions(
            gY=0.36,
            g2=0.64,
            g3=0.60,
            lambdaH=0.20,
            ye=ye,
            yu=yu,
            yd=yd,
            K=K,
        ).validated()

    def test_pack_round_trip_preserves_full_complex_yukawas(self) -> None:
        initial = self._initial()
        packed = _pack(initial)
        self.assertEqual(packed.shape, (76,))
        values = _unpack(packed)
        np.testing.assert_allclose(values[4], initial.ye)
        np.testing.assert_allclose(values[5], initial.yu)
        np.testing.assert_allclose(values[6], initial.yd)
        np.testing.assert_allclose(values[7], initial.K)

    def test_diagonal_limit_matches_previous_component_formula(self) -> None:
        ye = np.diag([3.0e-6, 6.0e-4, 1.0e-2]).astype(complex)
        yu = np.diag([1.0e-5, 7.0e-3, 0.75]).astype(complex)
        yd = np.diag([2.0e-5, 4.0e-4, 1.8e-2]).astype(complex)
        K = np.eye(3, dtype=complex) * 1.0e-14
        initial = SMInitialConditions(
            0.36, 0.64, 0.60, 0.20, ye, yu, yd, K
        ).validated()
        derivative = _beta(0.0, _pack(initial))
        _, _, _, _, beta_ye, beta_yu, beta_yd, _ = _unpack(
            derivative * LOOP
        )

        ye_diag = np.diag(ye).real
        yu_diag = np.diag(yu).real
        yd_diag = np.diag(yd).real
        T = np.sum(ye_diag**2) + 3*np.sum(yu_diag**2) + 3*np.sum(yd_diag**2)

        expected_yu = yu_diag * (
            1.5*(yu_diag**2 - yd_diag**2)
            + T
            - (17.0/12.0)*0.36**2
            - (9.0/4.0)*0.64**2
            - 8.0*0.60**2
        )
        expected_yd = yd_diag * (
            1.5*(yd_diag**2 - yu_diag**2)
            + T
            - (5.0/12.0)*0.36**2
            - (9.0/4.0)*0.64**2
            - 8.0*0.60**2
        )
        expected_ye = ye_diag * (
            1.5*ye_diag**2
            + T
            - (15.0/4.0)*0.36**2
            - (9.0/4.0)*0.64**2
        )

        np.testing.assert_allclose(np.diag(beta_yu).real, expected_yu)
        np.testing.assert_allclose(np.diag(beta_yd).real, expected_yd)
        np.testing.assert_allclose(np.diag(beta_ye).real, expected_ye)
        np.testing.assert_allclose(beta_yu - np.diag(np.diag(beta_yu)), 0.0)
        np.testing.assert_allclose(beta_yd - np.diag(np.diag(beta_yd)), 0.0)
        np.testing.assert_allclose(beta_ye - np.diag(np.diag(beta_ye)), 0.0)

    def test_final_boundary_accepts_full_yukawa_matrices(self) -> None:
        initial = self._initial()
        sm = SMNumericalState(
            gY=initial.gY,
            g2=initial.g2,
            g3=initial.g3,
            lambdaH=initial.lambdaH,
            ye=initial.ye,
            yu=initial.yu,
            yd=initial.yd,
        ).validated()
        boundary = FinalSMBoundaryState(
            mu_gev=1.0e9,
            sm=sm,
        ).validated()
        built = build_weinberg_initial_conditions(
            boundary,
            initial.K,
        )
        np.testing.assert_allclose(built.ye, initial.ye)
        np.testing.assert_allclose(built.yu, initial.yu)
        np.testing.assert_allclose(built.yd, initial.yd)

    def test_charged_lepton_basis_rotation_diagonalizes_left_hermitian(self) -> None:
        rng = np.random.default_rng(10)
        z = rng.normal(size=(3, 3)) + 1j*rng.normal(size=(3, 3))
        u, _ = np.linalg.qr(z)
        singular = np.array([3.0e-6, 6.0e-4, 1.0e-2])
        ye = u @ np.diag(singular)

        mass = np.array(
            [
                [1.0, 0.2, 0.1j],
                [0.2, 1.5, 0.3],
                [0.1j, 0.3, 2.0],
            ],
            dtype=complex,
        ) * 1.0e-11

        rotated = charged_lepton_mass_basis_matrix(mass, ye)
        self.assertTrue(
            np.allclose(
                rotated,
                rotated.T,
                rtol=1e-12,
                atol=1e-20,
            )
        )

        left, singular_values, _ = np.linalg.svd(ye)
        order = np.argsort(singular_values)
        left = left[:, order]
        hermitian = left.conj().T @ (ye @ ye.conj().T) @ left
        np.testing.assert_allclose(
            hermitian,
            np.diag(np.diag(hermitian)),
            rtol=1e-10,
            atol=1e-18,
        )

    def test_project_neutrino_mass_normalisation(self) -> None:
        K = np.eye(3, dtype=complex) * 1.0e-14
        m = neutrino_mass_matrix(K, vev_gev=246.22)
        np.testing.assert_allclose(
            m,
            -(246.22**2 / 2.0) * K,
        )


if __name__ == "__main__":
    unittest.main()
