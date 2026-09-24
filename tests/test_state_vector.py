from __future__ import annotations

import unittest

import numpy as np

from Numerical.State import (
    SMNumericalState,
    SharedT3UVState,
    T3Representation,
    T3UVState,
)
from Numerical.StateVector import (
    pack_uv_state,
    unpack_uv_state,
)


def _sm_state() -> SMNumericalState:
    return SMNumericalState(
        gY=0.36,
        g2=0.65,
        g3=1.1,
        lambdaH=0.26,
        yu=np.diag([1.0e-5, 7.0e-3, 0.9]).astype(complex),
        yd=np.diag([2.0e-5, 4.0e-4, 2.0e-2]).astype(complex),
        ye=np.diag([3.0e-6, 6.0e-4, 1.0e-2]).astype(complex),
    )


def _assert_matrix_equal(
    testcase: unittest.TestCase,
    left: np.ndarray,
    right: np.ndarray,
) -> None:
    np.testing.assert_allclose(left, right, rtol=0.0, atol=0.0)


class StateVectorTests(unittest.TestCase):
    def test_class_a_round_trip(self) -> None:
        state = T3UVState(
            mu_gev=1.0e12,
            representation=T3Representation(1, 3, 2, 0),
            sm=_sm_state(),
            y1=np.eye(3, dtype=complex) * (0.1 + 0.02j),
            y2=np.eye(3, dtype=complex) * (0.2 - 0.01j),
            MF=np.diag([1.0e10, 1.1e10, 1.2e10]).astype(complex),
            mS1Sq=(8.0e9) ** 2,
            mS2Sq=(7.0e9) ** 2,
            lambdaS1=0.1,
            lambdaS2=0.2,
            lambdaH1=0.01,
            lambdaH2=0.02,
            lambda12=0.03,
            lambdaT3=0.04 + 0.01j,
            lambdaH2Adj=0.005,
            lambdaS2Adj=0.006,
        ).validated()

        vector, layout = pack_uv_state(state)
        rebuilt = unpack_uv_state(vector, layout)

        self.assertEqual(layout.names[:7], (
            "gY", "g2", "g3", "yu", "yd", "ye", "lambdaH"
        ))
        self.assertEqual(rebuilt.mu_gev, state.mu_gev)
        self.assertEqual(rebuilt.representation, state.representation)
        _assert_matrix_equal(self, rebuilt.y1, state.y1)
        _assert_matrix_equal(self, rebuilt.y2, state.y2)
        _assert_matrix_equal(self, rebuilt.MF, state.MF)
        _assert_matrix_equal(self, rebuilt.sm.yu, state.sm.yu)
        self.assertEqual(rebuilt.lambdaT3, state.lambdaT3)

    def test_doublet_majorana_round_trip_keeps_complex_quartics(self) -> None:
        state = T3UVState(
            mu_gev=1.0e11,
            representation=T3Representation(2, 2, 1, -1),
            sm=_sm_state(),
            y1=np.eye(3, dtype=complex) * 0.1,
            y2=np.eye(3, dtype=complex) * 0.2,
            MF=np.diag([1.0e9, 1.1e9, 1.2e9]).astype(complex),
            mS1Sq=1.0e18,
            mS2Sq=0.9e18,
            lambdaS1=0.1,
            lambdaS2=0.2,
            lambdaH1=0.01,
            lambdaH2=0.02,
            lambda12=0.03,
            lambdaT3=0.04 + 0.03j,
            lambdaH1Adj=0.001,
            lambdaH2Adj=0.002,
            lambda12Adj=0.003,
            lambdaHHdagS2S2=0.01 + 0.02j,
            lambdaHHdagS1barS1bar=0.02 - 0.01j,
            lambdaS1bar2S2bar2=0.03 + 0.04j,
            lambdaS1barS2S2bar2=0.04 - 0.02j,
            lambdaS1S1bar2S2bar=0.05 + 0.01j,
            lambdaHHdagS1barS2barCross=0.06 - 0.03j,
        ).validated()

        vector, layout = pack_uv_state(state)
        rebuilt = unpack_uv_state(vector, layout)

        self.assertIn("lambdaHHdagS2S2", layout.names)
        self.assertEqual(
            rebuilt.lambdaHHdagS2S2,
            state.lambdaHHdagS2S2,
        )
        self.assertEqual(
            rebuilt.lambdaHHdagS1barS2barCross,
            state.lambdaHHdagS1barS2barCross,
        )

    def test_shared_scalar_round_trip(self) -> None:
        state = SharedT3UVState(
            mu_gev=1.0e12,
            representation=T3Representation(
                2, 2, 1, -1, shared_scalar=True
            ),
            sm=_sm_state(),
            h=np.eye(3, dtype=complex) * (0.1 + 0.02j),
            MF=np.diag([1.0e10, 1.1e10, 1.2e10]).astype(complex),
            mSSq=(8.0e9) ** 2,
            lambdaS=0.1,
            lambda3=0.01,
            lambda4=0.02,
            lambda5=0.03,
        ).validated()

        vector, layout = pack_uv_state(state)
        rebuilt = unpack_uv_state(vector, layout)

        self.assertTrue(layout.shared_scalar)
        self.assertEqual(
            layout.names[-7:],
            ("h", "MF", "mSSq", "lambdaS", "lambda3", "lambda4", "lambda5"),
        )
        _assert_matrix_equal(self, rebuilt.h, state.h)
        _assert_matrix_equal(self, rebuilt.MF, state.MF)

    def test_unpack_can_update_scale_without_changing_couplings(self) -> None:
        state = SharedT3UVState(
            mu_gev=1.0e12,
            representation=T3Representation(
                2, 2, 3, -1, shared_scalar=True
            ),
            sm=_sm_state(),
            h=np.eye(3, dtype=complex) * 0.1,
            MF=np.diag([1.0e10, 1.1e10, 1.2e10]).astype(complex),
            mSSq=1.0e18,
            lambdaS=0.1,
            lambda3=0.01,
            lambda4=0.02,
            lambda5=0.03,
        ).validated()

        vector, layout = pack_uv_state(state)
        rebuilt = unpack_uv_state(
            vector,
            layout,
            mu_gev=5.0e10,
        )

        self.assertEqual(rebuilt.mu_gev, 5.0e10)
        _assert_matrix_equal(self, rebuilt.h, state.h)

    def test_wrong_vector_length_is_rejected(self) -> None:
        state = SharedT3UVState(
            mu_gev=1.0e12,
            representation=T3Representation(
                2, 2, 1, -1, shared_scalar=True
            ),
            sm=_sm_state(),
            h=np.eye(3, dtype=complex) * 0.1,
            MF=np.eye(3, dtype=complex),
            mSSq=1.0,
            lambdaS=0.1,
            lambda3=0.01,
            lambda4=0.02,
            lambda5=0.03,
        ).validated()

        vector, layout = pack_uv_state(state)

        with self.assertRaisesRegex(ValueError, "layout expects"):
            unpack_uv_state(vector[:-1], layout)


if __name__ == "__main__":
    unittest.main()
