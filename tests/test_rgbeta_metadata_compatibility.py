from __future__ import annotations

import unittest

import numpy as np

from Numerical.IntermediateScalarRGBetaEvaluator import (
    _validate_intermediate_scalar_metadata,
)
from Numerical.RGBetaEvaluator import _validate_payload_metadata
from Numerical.IntermediateScalarState import T3IntermediateScalarState
from Numerical.State import (
    SMNumericalState,
    T3Representation,
    T3UVState,
)


def _sm() -> SMNumericalState:
    return SMNumericalState(
        gY=0.36,
        g2=0.65,
        g3=0.90,
        lambdaH=0.25,
        yu=np.diag([1.0e-5, 7.0e-3, 0.75]).astype(complex),
        yd=np.diag([2.0e-5, 4.0e-4, 1.8e-2]).astype(complex),
        ye=np.diag([3.0e-6, 6.0e-4, 1.0e-2]).astype(complex),
    ).validated()


def _rep() -> T3Representation:
    return T3Representation(
        d_s1=2,
        d_s2=2,
        d_f=1,
        alpha=-1,
        shared_scalar=False,
    ).validated()


def _uv_state() -> T3UVState:
    zero = 0.0 + 0.0j
    return T3UVState(
        mu_gev=1.0e12,
        representation=_rep(),
        sm=_sm(),
        y1=np.diag([0.05, 0.07, 0.09]).astype(complex),
        y2=np.diag([0.04, 0.06, 0.08]).astype(complex),
        MF=np.diag([1.0e10, 1.1e10, 1.2e10]).astype(complex),
        mS1Sq=6.4e19,
        mS2Sq=4.9e19,
        lambdaS1=0.10,
        lambdaS2=0.12,
        lambdaH1=0.02,
        lambdaH2=0.03,
        lambda12=0.01,
        lambdaT3=0.005 + 0.002j,
        lambdaH1Adj=0.0,
        lambdaH2Adj=0.0,
        lambda12Adj=0.0,
        lambdaHHdagS2S2=zero,
        lambdaHHdagS1barS1bar=zero,
        lambdaS1bar2S2bar2=zero,
        lambdaS1barS2S2bar2=zero,
        lambdaS1S1bar2S2bar=zero,
        lambdaHHdagS1barS2barCross=zero,
    ).validated()


def _intermediate_state() -> T3IntermediateScalarState:
    zero = 0.0 + 0.0j
    return T3IntermediateScalarState(
        mu_gev=7.0e9,
        representation=_rep(),
        sm=_sm(),
        mS1Sq=6.4e19,
        mS2Sq=4.9e19,
        lambdaS1=0.10,
        lambdaS2=0.12,
        lambdaH1=0.02,
        lambdaH2=0.03,
        lambda12=0.01,
        lambdaT3=0.005 + 0.002j,
        lambdaH1Adj=0.0,
        lambdaH2Adj=0.0,
        lambda12Adj=0.0,
        lambdaHHdagS2S2=zero,
        lambdaHHdagS1barS1bar=zero,
        lambdaS1bar2S2bar2=zero,
        lambdaS1barS2S2bar2=zero,
        lambdaS1S1bar2S2bar=zero,
        lambdaHHdagS1barS2barCross=zero,
    ).validated()


class RGBetaMetadataCompatibilityTests(unittest.TestCase):
    def test_ordinary_uv_payload_may_omit_sharedscalar(self) -> None:
        payload = {
            "metadata": {
                "dS1": 2,
                "dS2": 2,
                "dF": 1,
                "alpha": -1,
                "ReportBetaConvention": (
                    "16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g"
                ),
            }
        }
        _validate_payload_metadata(payload, _uv_state())

    def test_ordinary_intermediate_payload_may_omit_sharedscalar(self) -> None:
        payload = {
            "metadata": {
                "dS1": 2,
                "dS2": 2,
                "dF": 1,
                "alpha": -1,
                "IntegratedField": "F",
                "ActiveBSMFields": ["S1", "S2"],
                "ReportBetaConvention": (
                    "16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g"
                ),
            }
        }
        _validate_intermediate_scalar_metadata(payload, _intermediate_state())

    def test_explicit_sharedscalar_true_is_rejected_for_ordinary_state(self) -> None:
        payload = {
            "metadata": {
                "SharedScalar": True,
                "dS1": 2,
                "dS2": 2,
                "dF": 1,
                "alpha": -1,
                "ReportBetaConvention": (
                    "16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g"
                ),
            }
        }

        with self.assertRaisesRegex(ValueError, "SharedScalar"):
            _validate_payload_metadata(payload, _uv_state())


if __name__ == "__main__":
    unittest.main()
