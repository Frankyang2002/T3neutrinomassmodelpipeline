from __future__ import annotations

import unittest

import numpy as np

from Numerical.State import (
    SMNumericalState,
    SharedT3UVState,
    T3Representation,
    T3UVState,
)


def _sm_state() -> SMNumericalState:
    return SMNumericalState(
        gY=0.36,
        g2=0.65,
        g3=1.1,
        lambdaH=0.26,
        yu=np.diag([1.0e-5, 7.0e-3, 0.9]),
        yd=np.diag([2.0e-5, 4.0e-4, 2.0e-2]),
        ye=np.diag([3.0e-6, 6.0e-4, 1.0e-2]),
    )


class T3NumericalStateTests(unittest.TestCase):
    def test_ordinary_class_a_state_validates(self) -> None:
        rep = T3Representation(d_s1=1, d_s2=3, d_f=2, alpha=0)

        state = T3UVState(
            mu_gev=1.0e12,
            representation=rep,
            sm=_sm_state(),
            y1=np.eye(3, dtype=complex) * 0.1,
            y2=np.eye(3, dtype=complex) * 0.2,
            MF=np.diag([1.0e10, 1.1e10, 1.2e10]),
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

        self.assertEqual(state.y1.shape, (3, 3))
        self.assertEqual(state.MF.shape, (3, 3))
        self.assertFalse(state.representation.majorana_fermion)

    def test_class_e_requires_rg_closed_quartics(self) -> None:
        rep = T3Representation(d_s1=3, d_s2=3, d_f=2, alpha=0)

        with self.assertRaisesRegex(ValueError, "lambdaH1Adj is required"):
            T3UVState(
                mu_gev=1.0e12,
                representation=rep,
                sm=_sm_state(),
                y1=np.zeros((3, 3), dtype=complex),
                y2=np.zeros((3, 3), dtype=complex),
                MF=np.eye(3),
                mS1Sq=1.0,
                mS2Sq=1.0,
                lambdaS1=0.0,
                lambdaS2=0.0,
                lambdaH1=0.0,
                lambdaH2=0.0,
                lambda12=0.0,
                lambdaT3=0.0,
            ).validated()

    def test_shared_scalar_state_validates(self) -> None:
        rep = T3Representation(
            d_s1=2,
            d_s2=2,
            d_f=1,
            alpha=-1,
            shared_scalar=True,
        )

        state = SharedT3UVState(
            mu_gev=1.0e12,
            representation=rep,
            sm=_sm_state(),
            h=np.eye(3, dtype=complex) * 0.1,
            MF=np.diag([1.0e10, 1.1e10, 1.2e10]),
            mSSq=(8.0e9) ** 2,
            lambdaS=0.1,
            lambda3=0.01,
            lambda4=0.02,
            lambda5=0.03,
        ).validated()

        self.assertTrue(state.representation.shared_scalar)
        self.assertTrue(state.representation.majorana_fermion)

    def test_irrelevant_conditional_quartic_is_rejected(self) -> None:
        rep = T3Representation(d_s1=1, d_s2=3, d_f=2, alpha=0)

        with self.assertRaisesRegex(ValueError, "lambdaH1Adj is not present"):
            T3UVState(
                mu_gev=1.0e12,
                representation=rep,
                sm=_sm_state(),
                y1=np.zeros((3, 3), dtype=complex),
                y2=np.zeros((3, 3), dtype=complex),
                MF=np.eye(3),
                mS1Sq=1.0,
                mS2Sq=1.0,
                lambdaS1=0.0,
                lambdaS2=0.0,
                lambdaH1=0.0,
                lambdaH2=0.0,
                lambda12=0.0,
                lambdaT3=0.0,
                lambdaH1Adj=0.0,
                lambdaH2Adj=0.0,
                lambdaS2Adj=0.0,
            ).validated()


if __name__ == "__main__":
    unittest.main()
