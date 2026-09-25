from __future__ import annotations

import unittest

import numpy as np

from Numerical.WeinbergTrajectory import (
    neutrino_mass_from_c5,
    run_weinberg_trajectory,
    scale_dependent_neutrino_observables,
)
from Numerical.WeinbergRunning import SMInitialConditions


def _initial() -> SMInitialConditions:
    c5 = np.array(
        [
            [1.0e-15, 2.0e-16, 0.0],
            [2.0e-16, 3.0e-15, 1.0e-16],
            [0.0, 1.0e-16, 5.0e-15],
        ],
        dtype=complex,
    )

    return SMInitialConditions(
        gY=0.36,
        g2=0.63,
        g3=0.90,
        lambdaH=0.24,
        ye=np.diag([3.0e-6, 6.0e-4, 1.0e-2]).astype(complex),
        yu=np.diag([1.0e-5, 7.0e-3, 0.75]).astype(complex),
        yd=np.diag([2.0e-5, 4.0e-4, 1.8e-2]).astype(complex),
        K=c5,
    ).validated()


class WeinbergTrajectoryTests(unittest.TestCase):
    def test_requested_scales_are_retained(self) -> None:
        initial = _initial()

        result = run_weinberg_trajectory(
            initial,
            1.0e5,
            1.0e3,
            save_scales_gev=[1.0e5, 1.0e4, 1.0e3],
        )

        np.testing.assert_allclose(
            result.mu_gev,
            np.array([1.0e5, 1.0e4, 1.0e3]),
            rtol=0.0,
            atol=0.0,
        )
        self.assertEqual(result.n_points, 3)

    def test_c5_stays_symmetric_at_saved_points(self) -> None:
        result = run_weinberg_trajectory(
            _initial(),
            1.0e5,
            1.0e3,
            save_scales_gev=[1.0e5, 1.0e4, 1.0e3],
        )

        for index in range(result.n_points):
            c5 = result.point_at_index(index).c5
            np.testing.assert_allclose(c5, c5.T)

    def test_neutrino_mass_normalisation(self) -> None:
        c5 = _initial().K
        vev = 246.22

        mass = neutrino_mass_from_c5(
            c5,
            vev_gev=vev,
        )

        np.testing.assert_allclose(
            mass,
            -(vev**2 / 2.0) * c5,
        )

    def test_observables_are_available_at_every_scale(self) -> None:
        result = run_weinberg_trajectory(
            _initial(),
            1.0e5,
            1.0e3,
            save_scales_gev=[1.0e5, 1.0e4, 1.0e3],
        )

        points = scale_dependent_neutrino_observables(
            result,
            ordering="NO",
        )

        self.assertEqual(len(points), 3)

        for point in points:
            self.assertEqual(point.observables.masses_ev.shape, (3,))
            self.assertEqual(point.observables.pmns_abs.shape, (3, 3))
            self.assertLess(point.observables.takagi_residual, 1.0e-7)
            self.assertTrue(
                np.all(point.observables.masses_ev >= 0.0)
            )

    def test_running_changes_c5(self) -> None:
        result = run_weinberg_trajectory(
            _initial(),
            1.0e5,
            1.0e3,
            save_scales_gev=[1.0e5, 1.0e3],
        )

        self.assertFalse(
            np.allclose(
                result.initial_point.c5,
                result.final_point.c5,
                rtol=1.0e-12,
                atol=0.0,
            )
        )


if __name__ == "__main__":
    unittest.main()
