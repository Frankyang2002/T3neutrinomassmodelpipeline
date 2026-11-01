"""Evaluate standard and general scalar-only EFT1 RGBeta payloads."""
from __future__ import annotations

from typing import Any, Mapping

from Numerical.core.BetaVector import beta_values_to_derivative, validate_beta_keys
from Numerical.core.RGBetaEvaluator import evaluate_inputform_expression
from Numerical.running.IntermediateScalarState import (
    SharedT3IntermediateScalarState, T3IntermediateScalarState,
)
from Numerical.running.IntermediateScalarStateVector import (
    IntermediateScalarState, pack_intermediate_scalar_state,
)

EXPECTED_CONVENTION = '16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g'
GENERAL_FORMAT = 'general-t3-quadratic-rgbeta-experimental-v1'
FERMION_PARAMETERS = frozenset({'y1', 'y2', 'MF', 'y1Conj', 'y2Conj'})


def intermediate_scalar_state_environment(state: IntermediateScalarState) -> dict[str, Any]:
    state = state.validated()
    env = {
        'gY': state.sm.gY, 'g2': state.sm.g2, 'g3': state.sm.g3,
        'yu': state.sm.yu, 'yd': state.sm.yd, 'ye': state.sm.ye,
        'lambdaH': state.sm.lambdaH,
    }
    if isinstance(state, SharedT3IntermediateScalarState):
        env.update({n: getattr(state, n) for n in
                    ('mSSq', 'lambdaS', 'lambda3', 'lambda4', 'lambda5')})
        return env
    for n in (
        'mS1Sq', 'mS2Sq', 'lambdaS1', 'lambdaS2', 'lambdaH1',
        'lambdaH2', 'lambda12', 'lambdaT3', 'lambdaH1Adj',
        'lambdaH2Adj', 'lambdaS1Adj', 'lambdaS2Adj', 'lambda12Adj',
        'lambda12Cross', 'lambdaHHdagS2S2', 'lambdaHHdagS1barS1bar',
        'lambdaS1bar2S2bar2', 'lambdaS1barS2S2bar2',
        'lambdaS1S1bar2S2bar', 'lambdaHHdagS1barS2barCross',
    ):
        v = getattr(state, n)
        if v is not None:
            env[n] = v
    env.update(state.general_couplings)
    return env


def _validate_intermediate_scalar_metadata(
    payload: Mapping[str, Any], state: IntermediateScalarState,
) -> None:
    meta = payload.get('metadata', {})
    if not isinstance(meta, Mapping):
        raise ValueError('EFT1 RGBeta metadata must be a mapping')
    rep = state.representation.validated()
    general = payload.get('format') == GENERAL_FORMAT

    expected = {
        'SharedScalar': rep.shared_scalar,
        'dS1': rep.d_s1, 'dS2': rep.d_s2,
        'dF': rep.d_f, 'alpha': rep.alpha,
    }
    actual = {**meta, 'SharedScalar': bool(meta.get('SharedScalar', False))}
    for key, value in expected.items():
        if actual.get(key) != value:
            raise ValueError(
                f'EFT1 RGBeta {key}: payload={actual.get(key)!r}, state={value!r}'
            )

    if general:
        # The experimental general exporter identifies the stage by
        # GeneralScalarStage, but does not emit IntegratedField or
        # ActiveBSMFields. Its EFT1 theory is built by T3RGBetaBuildEFT1.
        if not isinstance(state, T3IntermediateScalarState):
            raise ValueError('General EFT1 requires ordinary split scalar state')
        if meta.get('GeneralScalarStage') != 'EFT1':
            raise ValueError('General EFT1 stage metadata must be EFT1')
        expected_model = None
        dims = (rep.d_s1, rep.d_s2, rep.d_f)
        families = {
            (1, 3, 2): 'A', (2, 2, 1): 'B', (2, 2, 3): 'C',
            (3, 1, 2): 'D', (3, 3, 2): 'E',
        }
        if dims in families:
            expected_model = f'{families[dims]}_{rep.alpha:+d}'
        if meta.get('GeneralScalarModel') != expected_model:
            raise ValueError('General EFT1 model metadata does not match state')
        if meta.get('IntegratedField') not in (None, 'F'):
            raise ValueError('General EFT1 IntegratedField contradicts F threshold')
        if 'ActiveBSMFields' in meta and list(meta['ActiveBSMFields']) != ['S1', 'S2']:
            raise ValueError('General EFT1 ActiveBSMFields contradicts scalar-only EFT')
    else:
        if meta.get('IntegratedField') != 'F':
            raise ValueError(
                f"EFT1 RGBeta IntegratedField: payload={meta.get('IntegratedField')!r}, state='F'"
            )
        fields = ['S'] if rep.shared_scalar else ['S1', 'S2']
        if list(meta.get('ActiveBSMFields', [])) != fields:
            raise ValueError('Unexpected EFT1 ActiveBSMFields')

    if meta.get('ReportBetaConvention') != EXPECTED_CONVENTION:
        raise ValueError('Unexpected EFT1 RGBeta convention')


def evaluate_intermediate_scalar_rgbeta_payload(
    payload: Mapping[str, Any], state: IntermediateScalarState,
) -> dict[str, Any]:
    fmt = payload.get('format')
    status = payload.get('status')
    if fmt == GENERAL_FORMAT:
        if status != 'SuccessCandidateMassBeta':
            raise ValueError('General EFT1 RGBeta payload status mismatch')
        if not isinstance(state, T3IntermediateScalarState):
            raise ValueError('General EFT1 requires ordinary split scalar state')
    elif status != 'Success':
        raise ValueError('EFT1 RGBeta status is not Success')

    state = state.validated()
    _validate_intermediate_scalar_metadata(payload, state)
    betas = payload.get('report_betas')
    if not isinstance(betas, Mapping):
        raise ValueError('EFT1 report_betas missing')
    if fmt == GENERAL_FORMAT:
        forbidden = FERMION_PARAMETERS.intersection(betas)
        if forbidden:
            raise ValueError(
                f'General EFT1 still contains integrated-out fermion betas: {sorted(forbidden)}'
            )
    _, layout = pack_intermediate_scalar_state(state)
    validate_beta_keys({n: None for n in betas}, layout)
    env = intermediate_scalar_state_environment(state)
    result = {
        str(n): evaluate_inputform_expression(str(expr), env)
        for n, expr in betas.items()
    }
    validate_beta_keys(result, layout)
    return result


def intermediate_scalar_derivative_from_payload(
    payload: Mapping[str, Any], state: IntermediateScalarState,
):
    state = state.validated()
    _, layout = pack_intermediate_scalar_state(state)
    result = evaluate_intermediate_scalar_rgbeta_payload(payload, state)
    return beta_values_to_derivative(
        result, layout, divide_by_loop_factor=True,
    ), layout


EFT1State = IntermediateScalarState
eft1_state_environment = intermediate_scalar_state_environment
evaluate_eft1_rgbeta_payload = evaluate_intermediate_scalar_rgbeta_payload
eft1_derivative_from_payload = intermediate_scalar_derivative_from_payload
