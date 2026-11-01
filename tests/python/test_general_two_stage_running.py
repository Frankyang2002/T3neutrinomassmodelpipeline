"""Two-stage ODE regression using the existing UV/EFT1 production runners.

Synthetic beta expressions test plumbing, not the physics of a Wolfram export.
"""
from __future__ import annotations

import copy
import numpy as np
import pytest

from Numerical.core.State import SMNumericalState, T3Representation, T3UVState
from Numerical.core.StateVector import pack_uv_state
from Numerical.running.UVRunner import run_uv_segment
from Numerical.running.FermionThresholdBoundary import build_intermediate_scalar_boundary
from Numerical.running.IntermediateScalarRunner import run_intermediate_scalar_segment
from Numerical.running.IntermediateScalarStateVector import pack_intermediate_scalar_state

CONVENTION = "16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g"
Q = "lambdaT3GeneralQ000004I1"
M = "muT3GeneralS2Sq"


def initial_state():
    ones = np.eye(3, dtype=complex)
    sm = SMNumericalState(
        gY=.35, g2=.64, g3=1.05, lambdaH=.13,
        yu=.3*ones, yd=.02*ones, ye=.01*ones,
    )
    return T3UVState(
        mu_gev=1e7,
        representation=T3Representation(1, 3, 2, -2),
        sm=sm, y1=.08*ones, y2=.09*ones, MF=1e5*ones,
        mS1Sq=1e6, mS2Sq=2e6,
        lambdaS1=.2, lambdaS2=.3, lambdaH1=.04,
        lambdaH2=.05, lambda12=.06, lambdaT3=.01+.02j,
        lambdaH2Adj=.01, lambdaS2Adj=.01,
        general_couplings={Q:.02+.03j, M:100+25j},
    ).validated()


def payload(state, stage):
    if stage == 'UV':
        _, layout = pack_uv_state(state)
    else:
        _, layout = pack_intermediate_scalar_state(state)
    # Every block has an explicit beta expression. Real matrices stay matrices.
    betas = {name: f"0*{name}" for name in layout.names}
    betas[Q] = f"0.5*{Q}"
    betas[M] = f"0.25*{M}"
    rep = state.representation
    meta = {
        "SharedScalar":False, "dS1":rep.d_s1, "dS2":rep.d_s2,
        "dF":rep.d_f,"alpha":rep.alpha,
        "ReportBetaConvention":CONVENTION,
        "GeneralScalarStage":stage,
    }
    if stage == 'EFT1':
        meta.update(IntegratedField="F",ActiveBSMFields=["S1","S2"])
    return {
        "format":"general-t3-quadratic-rgbeta-experimental-v1",
        "status":"SuccessCandidateMassBeta", "metadata":meta,
        "report_betas":betas,
    }


def test_general_couplings_survive_uv_eft1_numerical_running():
    uv0 = initial_state()
    uv = run_uv_segment(uv0, payload(uv0,'UV'), 1e5)
    assert uv.solver_success
    f_state = uv.final_state
    eft0 = build_intermediate_scalar_boundary(f_state)
    assert set(eft0.general_couplings) == {Q,M}
    assert eft0.general_couplings[Q] == f_state.general_couplings[Q]
    assert eft0.general_couplings[M] == f_state.general_couplings[M]
    eft = run_intermediate_scalar_segment(eft0, payload(eft0,'EFT1'), 1e3)
    assert eft.solver_success
    final = eft.final_state
    tspan = np.log(1e3 / 1e7)
    assert np.allclose(final.general_couplings[Q],
                       uv0.general_couplings[Q]*np.exp(.5*tspan/(16*np.pi**2)),rtol=1e-7)
    assert np.allclose(final.general_couplings[M],
                       uv0.general_couplings[M]*np.exp(.25*tspan/(16*np.pi**2)),rtol=1e-7)
    assert np.allclose(final.sm.yu,uv0.sm.yu)
    assert np.allclose(final.mS1Sq,uv0.mS1Sq)


def test_missing_general_beta_fails_closed():
    uv0 = initial_state()
    bad = payload(uv0,'UV')
    del bad['report_betas'][Q]
    with pytest.raises(ValueError,match='missing|Missing|mismatch|coverage'):
        run_uv_segment(uv0,bad,1e5)


def test_wrong_stage_is_rejected():
    uv0 = initial_state()
    wrong = payload(uv0,'EFT1')
    with pytest.raises(ValueError,match='UV|stage'):
        run_uv_segment(uv0,wrong,1e5)
