"""Checks the real-export benchmark state builder (no fake beta formulas)."""
import numpy as np
import pytest
from Numerical.orchestration.GeneralRunningBenchmark import initial_state, _encode
from Numerical.core.StateVector import pack_uv_state


def payload():
    names = ["gY", "g2", "g3", "lambdaH", "yu", "yd", "ye", "y1", "y2", "MF",
             "mS1Sq", "mS2Sq", "lambdaS1", "lambdaS2", "lambdaH1", "lambdaH2",
             "lambda12", "lambdaT3", "lambdaH2Adj", "lambdaS2Adj",
             "lambdaT3GeneralQ000004I1", "muT3GeneralS2Sq"]
    return {"metadata": {"dS1": 1, "dS2": 3, "dF": 2, "alpha": -2,
                         "GeneralScalarModel": "A_-2"},
            "report_betas": {k: "0" for k in names}}


def test_initial_state_contains_all_rgbeta_parameters():
    p = payload()
    model, state = initial_state(p, mu_uv=1e7, mu_f=1e5, mu_s=1e3)
    _, layout = pack_uv_state(state)
    assert model == "A_-2"
    assert set(layout.names) == set(p["report_betas"])
    assert state.general_couplings["muT3GeneralS2Sq"].imag != 0


def test_user_override_complex_general_coupling():
    _, state = initial_state(payload(), mu_uv=1e7, mu_f=1e5, mu_s=1e3,
        overrides={"lambdaT3GeneralQ000004I1": {"real": 0.06, "imag": -0.015}})
    assert state.general_couplings["lambdaT3GeneralQ000004I1"] == 0.06-0.015j
    assert _encode(state.general_couplings["lambdaT3GeneralQ000004I1"])["imag"] == -0.015


def test_missing_required_beta_rejected():
    p = payload()
    del p["report_betas"]["lambdaH2Adj"]
    with pytest.raises(ValueError, match="required structural quartic"):
        initial_state(p, mu_uv=1e7, mu_f=1e5, mu_s=1e3)
