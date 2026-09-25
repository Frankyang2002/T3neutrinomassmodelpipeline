"""Evaluate renormalisable RGBeta beta functions in the scalar-only intermediate EFT.

The parser itself is shared with ``Numerical.RGBetaEvaluator`` because the
current UV and intermediate Wolfram runners serialise expressions using the
same ``InputForm`` structures. This module supplies the intermediate scalar
state environment and stage-specific metadata checks.
"""

from __future__ import annotations

from typing import Any, Mapping

import numpy as np

from Numerical.BetaVector import beta_values_to_derivative, validate_beta_keys
from Numerical.IntermediateScalarState import (
    SharedT3IntermediateScalarState,
    T3IntermediateScalarState,
)
from Numerical.IntermediateScalarStateVector import (
    IntermediateScalarState,
    pack_intermediate_scalar_state,
)
from Numerical.RGBetaEvaluator import evaluate_inputform_expression


EXPECTED_CONVENTION = (
    "16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g"
)


def intermediate_scalar_state_environment(
    state: IntermediateScalarState,
) -> dict[str, Any]:
    """Return the numerical symbol environment for one scalar-only state."""

    state = state.validated()

    environment: dict[str, Any] = {
        "gY": state.sm.gY,
        "g2": state.sm.g2,
        "g3": state.sm.g3,
        "yu": state.sm.yu,
        "yd": state.sm.yd,
        "ye": state.sm.ye,
        "lambdaH": state.sm.lambdaH,
    }

    if isinstance(state, SharedT3IntermediateScalarState):
        environment.update(
            {
                "mSSq": state.mSSq,
                "lambdaS": state.lambdaS,
                "lambda3": state.lambda3,
                "lambda4": state.lambda4,
                "lambda5": state.lambda5,
            }
        )
        return environment

    environment.update(
        {
            "mS1Sq": state.mS1Sq,
            "mS2Sq": state.mS2Sq,
            "lambdaS1": state.lambdaS1,
            "lambdaS2": state.lambdaS2,
            "lambdaH1": state.lambdaH1,
            "lambdaH2": state.lambdaH2,
            "lambda12": state.lambda12,
            "lambdaT3": state.lambdaT3,
        }
    )

    for name in (
        "lambdaH1Adj",
        "lambdaH2Adj",
        "lambdaS1Adj",
        "lambdaS2Adj",
        "lambda12Adj",
        "lambda12Cross",
        "lambdaHHdagS2S2",
        "lambdaHHdagS1barS1bar",
        "lambdaS1bar2S2bar2",
        "lambdaS1barS2S2bar2",
        "lambdaS1S1bar2S2bar",
        "lambdaHHdagS1barS2barCross",
    ):
        value = getattr(state, name)
        if value is not None:
            environment[name] = value

    return environment


def _validate_intermediate_scalar_metadata(
    payload: Mapping[str, Any],
    state: IntermediateScalarState,
) -> None:
    metadata = payload.get("metadata", {})
    if not isinstance(metadata, Mapping):
        raise ValueError("RGBeta intermediate payload metadata must be an object.")

    rep = state.representation.validated()

    # Historical ordinary-T3 exports omitted SharedScalar. Absence is
    # interpreted as False for the ordinary branch. A shared-scalar payload
    # must still declare SharedScalar=True explicitly.
    actual_shared = bool(metadata.get("SharedScalar", False))

    expected = {
        "SharedScalar": rep.shared_scalar,
        "dS1": rep.d_s1,
        "dS2": rep.d_s2,
        "dF": rep.d_f,
        "alpha": rep.alpha,
        "IntegratedField": "F",
    }

    actual_values = {
        "SharedScalar": actual_shared,
        "dS1": metadata.get("dS1"),
        "dS2": metadata.get("dS2"),
        "dF": metadata.get("dF"),
        "alpha": metadata.get("alpha"),
        "IntegratedField": metadata.get("IntegratedField"),
    }

    for key, expected_value in expected.items():
        if actual_values[key] != expected_value:
            raise ValueError(
                "RGBeta intermediate payload does not match the numerical state: "
                f"{key}: payload={actual_values[key]!r}, "
                f"state={expected_value!r}"
            )

    expected_active = ["S"] if rep.shared_scalar else ["S1", "S2"]
    if list(metadata.get("ActiveBSMFields", [])) != expected_active:
        raise ValueError(
            "RGBeta intermediate payload has unexpected ActiveBSMFields: "
            f"{metadata.get('ActiveBSMFields')!r}."
        )

    if metadata.get("ReportBetaConvention") != EXPECTED_CONVENTION:
        raise ValueError(
            "Unsupported RGBeta intermediate report-beta convention."
        )


def evaluate_intermediate_scalar_rgbeta_payload(
    payload: Mapping[str, Any],
    state: IntermediateScalarState,
) -> dict[str, Any]:
    """Evaluate all renormalisable scalar-only ``report_betas``."""

    if payload.get("status") != "Success":
        raise ValueError("RGBeta intermediate payload status is not Success.")

    state = state.validated()
    _validate_intermediate_scalar_metadata(payload, state)

    report_betas = payload.get("report_betas")
    if not isinstance(report_betas, Mapping):
        raise ValueError(
            "RGBeta intermediate payload does not contain report_betas."
        )

    environment = intermediate_scalar_state_environment(state)

    evaluated = {
        str(name): evaluate_inputform_expression(
            str(expression),
            environment,
        )
        for name, expression in report_betas.items()
    }

    _, layout = pack_intermediate_scalar_state(state)
    validate_beta_keys(evaluated, layout)

    return evaluated


def intermediate_scalar_derivative_from_payload(
    payload: Mapping[str, Any],
    state: IntermediateScalarState,
) -> tuple[np.ndarray, Any]:
    """Return ``dy/dln(mu)`` for the renormalisable scalar-only state."""

    state = state.validated()
    _, layout = pack_intermediate_scalar_state(state)
    evaluated = evaluate_intermediate_scalar_rgbeta_payload(payload, state)

    derivative = beta_values_to_derivative(
        evaluated,
        layout,
        divide_by_loop_factor=True,
    )

    return derivative, layout


# Transitional aliases for existing callers/configuration tests.  Serialized
# EFT1 metadata keys are intentionally unchanged for report compatibility.
EFT1State = IntermediateScalarState
eft1_state_environment = intermediate_scalar_state_environment
evaluate_eft1_rgbeta_payload = evaluate_intermediate_scalar_rgbeta_payload
eft1_derivative_from_payload = intermediate_scalar_derivative_from_payload
