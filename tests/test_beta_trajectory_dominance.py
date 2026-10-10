"""No-physics tests for scale-dependent beta dominance formatting."""
from __future__ import annotations

from types import SimpleNamespace
import numpy as np

from Numerical.diagnostics.BetaTrajectoryDominance import analyse_states, COUPLINGS


def test_analyse_states_uses_complete_state_at_every_point(monkeypatch):
    from Numerical.diagnostics import BetaTrajectoryDominance as module
    class FakeResult:
        mu_gev = np.array([1e7, 1e5])
        n_points = 2
        def state_at_index(self, i):
            return SimpleNamespace(value=i)
    observed = []
    monkeypatch.setattr(module, "evaluate_rgbeta_payload",
                        lambda payload, state: {k: state.value for k in COUPLINGS})
    monkeypatch.setattr(module, "state_environment",
                        lambda state: {"value": state.value})
    def fake_analyse(expression, env, expected, *, target):
        observed.append((expression, env["value"], expected, target))
        return {"cancellation_ratio": 1.0, "categories": {
            "gauge": {"dominance": 1.0, "norm_16pi2_beta": 1.0}}}
    monkeypatch.setattr(module, "analyse_beta", fake_analyse)
    result = analyse_states(FakeResult(), {
        "report_betas": {k: k for k in COUPLINGS}})
    assert len(result) == 2
    assert result[0]["mu_gev"] == 1e7
    assert result[1]["mu_gev"] == 1e5
    assert len(observed) == 6
    assert observed[-1] == ("lambdaT3", 1, 1, "lambdaT3")


def test_representatives_cover_five_classes():
    from Numerical.diagnostics.BetaTrajectoryDominance import REPRESENTATIVE
    assert set(REPRESENTATIVE) == set("ABCDE")
    assert len(set(REPRESENTATIVE.values())) == 5
