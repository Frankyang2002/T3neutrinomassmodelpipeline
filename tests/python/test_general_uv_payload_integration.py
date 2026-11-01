"""Batch 2: native general-RGBeta payload acceptance by the existing UV evaluator."""
import numpy as np
import pytest
from Numerical.core.State import SMNumericalState, T3Representation, T3UVState
from Numerical.core.StateVector import pack_uv_state
from Numerical.core.RGBetaEvaluator import evaluate_rgbeta_payload, derivative_from_rgbeta_payload


def state():
    z = np.zeros((3, 3), complex)
    sm = SMNumericalState(.35, .65, 1.2, .13, z, z, z)
    return T3UVState(
        1e7, T3Representation(1, 3, 2, -2), sm,
        z, z, np.eye(3, dtype=complex) * 1e5,
        1e6, 4e6, .1, .2, .3, .4, .5, .01+.02j,
        lambdaH2Adj=.02, lambdaS2Adj=.03,
        general_couplings={
            'lambdaT3GeneralQ000004I1': .02+.03j,
            'muT3GeneralS2Sq': 1.5+2.5j,
            'y2Conj': np.eye(3, dtype=complex)*(.01+.04j),
        },
    )


def payload(s):
    _, layout = pack_uv_state(s)
    return {
        'format': 'general-t3-quadratic-rgbeta-experimental-v1',
        'status': 'SuccessCandidateMassBeta',
        'metadata': {
            'SharedScalar': False, 'dS1': 1, 'dS2': 3, 'dF': 2, 'alpha': -2,
            'GeneralScalarStage': 'UV',
            'ReportBetaConvention': '16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g',
        },
        'report_betas': {name: ('0*y2Conj' if name == 'y2Conj' else
                                '0*MF' if name == 'MF' else
                                '0*y1' if name == 'y1' else
                                '0*y2' if name == 'y2' else
                                '0*yu' if name == 'yu' else
                                '0*yd' if name == 'yd' else
                                '0*ye' if name == 'ye' else '0')
                         for name in layout.names},
    }


def test_general_payload_accepted_and_derivative_shape():
    s = state()
    p = payload(s)
    values = evaluate_rgbeta_payload(p, s)
    assert set(values) == set(pack_uv_state(s)[1].names)
    dy, layout = derivative_from_rgbeta_payload(p, s)
    assert dy.shape == (layout.size,)
    np.testing.assert_allclose(dy, 0.0)


def test_general_payload_rejects_missing_extra_beta():
    s = state()
    p = payload(s)
    del p['report_betas']['muT3GeneralS2Sq']
    with pytest.raises(ValueError, match='missing state variables'):
        evaluate_rgbeta_payload(p, s)


def test_general_payload_rejects_incorrect_stage():
    s = state()
    p = payload(s)
    p['metadata']['GeneralScalarStage'] = 'EFT1'
    with pytest.raises(ValueError, match='UV-stage'):
        evaluate_rgbeta_payload(p, s)


def test_restricted_payload_remains_strict():
    s = state()
    p = payload(s)
    p['format'] = 'legacy'
    with pytest.raises(ValueError, match='status is not Success'):
        evaluate_rgbeta_payload(p, s)
