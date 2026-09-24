from __future__ import annotations

import unittest

import numpy as np

from Numerical.EFT1State import (
    SharedT3EFT1State,
    T3EFT1State,
)
from Numerical.FinalSMBoundary import (
    build_weinberg_initial_conditions,
    project_after_scalar_threshold,
    scalar_threshold_masses,
)
from Numerical.State import (
    SMNumericalState,
    T3Representation,
)


def _sm() -> SMNumericalState:
    return SMNumericalState(
        gY=0.36,
        g2=0.65,
        g3=1.05,
        lambdaH=0.25,
        yu=np.diag([1.0e-5, 7.0e-3, 0.85]).astype(complex),
        yd=np.diag([2.0e-5, 4.0e-4, 1.8e-2]).astype(complex),
        ye=np.diag([3.0e-6, 6.0e-4, 1.0e-2]).astype(complex),
    )


def _ordinary() -> T3EFT1State:
    return T3EFT1State(
        mu_gev=7.0e9,
        representation=T3Representation(1, 3, 2, 0),
        sm=_sm(),
        mS1Sq=(8.0e9) ** 2,
        mS2Sq=(7.0e9) ** 2,
        lambdaS1=0.10,
        lambdaS2=0.12,
        lambdaH1=0.02,
        lambdaH2=0.03,
        lambda12=0.01,
        lambdaT3=0.005 + 0.002j,
        lambdaH2Adj=0.004,
        lambdaS2Adj=0.006,
    ).validated()


class FinalSMBoundaryTests(unittest.TestCase):
    def test_split_scalar_masses_are_extracted(self) -> None:
        state = _ordinary()
        diagnostic = scalar_threshold_masses(state)

        np.testing.assert_allclose(
            diagnostic.masses_gev,
            np.array([8.0e9, 7.0e9]),
            rtol=0.0,
            atol=0.0,
        )
        self.assertFalse(diagnostic.shared_scalar)
        self.assertFalse(diagnostic.degenerate)

    def test_shared_scalar_mass_is_extracted(self) -> None:
        state = SharedT3EFT1State(
            mu_gev=6.0e9,
            representation=T3Representation(
                2, 2, 1, -1, shared_scalar=True
            ),
            sm=_sm(),
            mSSq=(6.0e9) ** 2,
            lambdaS=0.1,
            lambda3=0.01,
            lambda4=0.02,
            lambda5=0.03,
        ).validated()

        diagnostic = scalar_threshold_masses(state)

        self.assertTrue(diagnostic.shared_scalar)
        self.assertTrue(diagnostic.degenerate)
        self.assertEqual(diagnostic.masses_gev, (6.0e9,))

    def test_negative_mass_squared_is_rejected_for_threshold(self) -> None:
        state = _ordinary()

        broken = T3EFT1State(
            mu_gev=state.mu_gev,
            representation=state.representation,
            sm=state.sm,
            mS1Sq=-1.0,
            mS2Sq=state.mS2Sq,
            lambdaS1=state.lambdaS1,
            lambdaS2=state.lambdaS2,
            lambdaH1=state.lambdaH1,
            lambdaH2=state.lambdaH2,
            lambda12=state.lambda12,
            lambdaT3=state.lambdaT3,
            lambdaH2Adj=state.lambdaH2Adj,
            lambdaS2Adj=state.lambdaS2Adj,
        ).validated()

        with self.assertRaisesRegex(
            ValueError,
            "positive to define a physical scalar threshold",
        ):
            scalar_threshold_masses(broken)

    def test_projection_keeps_only_sm_boundary_data(self) -> None:
        state = _ordinary()
        boundary = project_after_scalar_threshold(state)

        self.assertEqual(boundary.mu_gev, state.mu_gev)
        self.assertEqual(boundary.sm.gY, state.sm.gY)
        self.assertFalse(hasattr(boundary, "mS1Sq"))
        self.assertFalse(hasattr(boundary, "lambdaT3"))

    def test_projection_rejects_scale_jump(self) -> None:
        state = _ordinary()

        with self.assertRaisesRegex(
            ValueError,
            "must first be evolved",
        ):
            project_after_scalar_threshold(
                state,
                matching_scale_gev=5.0e9,
            )

    def test_build_weinberg_initial_conditions(self) -> None:
        boundary = project_after_scalar_threshold(_ordinary())

        c5 = np.array(
            [
                [1.0e-14, 2.0e-15 + 1.0e-15j, 0.0],
                [2.0e-15 + 1.0e-15j, 3.0e-14, 1.0e-15],
                [0.0, 1.0e-15, 5.0e-14],
            ],
            dtype=complex,
        )

        initial = build_weinberg_initial_conditions(
            boundary,
            c5,
        )

        np.testing.assert_allclose(initial.K, c5)
        np.testing.assert_allclose(
            initial.yu,
            boundary.sm.yu,
        )
        np.testing.assert_allclose(
            initial.yd,
            boundary.sm.yd,
        )
        np.testing.assert_allclose(
            initial.ye,
            boundary.sm.ye,
        )

    def test_nondiagonal_yukawa_is_preserved(self) -> None:
        state = _ordinary()
        yu = state.sm.yu.copy()
        yu[0, 1] = 1.0e-4 + 2.0e-5j

        modified_sm = SMNumericalState(
            gY=state.sm.gY,
            g2=state.sm.g2,
            g3=state.sm.g3,
            lambdaH=state.sm.lambdaH,
            yu=yu,
            yd=state.sm.yd,
            ye=state.sm.ye,
        ).validated()

        modified = T3EFT1State(
            mu_gev=state.mu_gev,
            representation=state.representation,
            sm=modified_sm,
            mS1Sq=state.mS1Sq,
            mS2Sq=state.mS2Sq,
            lambdaS1=state.lambdaS1,
            lambdaS2=state.lambdaS2,
            lambdaH1=state.lambdaH1,
            lambdaH2=state.lambdaH2,
            lambda12=state.lambda12,
            lambdaT3=state.lambdaT3,
            lambdaH2Adj=state.lambdaH2Adj,
            lambdaS2Adj=state.lambdaS2Adj,
        ).validated()

        boundary = project_after_scalar_threshold(modified)
        initial = build_weinberg_initial_conditions(
            boundary,
            np.eye(3, dtype=complex) * 1.0e-14,
        )

        np.testing.assert_allclose(initial.yu, yu)


if __name__ == "__main__":
    unittest.main()
