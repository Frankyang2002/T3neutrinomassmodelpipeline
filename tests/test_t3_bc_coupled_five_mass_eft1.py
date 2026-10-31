"""Unit checks for the opt-in native-RGBeta B/C five-mass EFT1 adapter.

Run from repository root:
  python -X utf8 -m pytest -q tests/test_t3_bc_coupled_five_mass_eft1.py
"""
from pathlib import Path
import math

import numpy as np
import pytest

from scripts import RunT3BCFullScalarMassEFT1 as coupled
from scripts import VerifyT3BCFullScalarMassClosure as closure


@pytest.fixture
def example_betas():
    layout = ("gY", "mS1Sq", "mS2Sq", "lambdaH")
    old = {
        "gY": "gY^3", "lambdaH": "lambdaH^2",
        "mS1Sq": "WRONG_OLD_M1", "mS2Sq": "WRONG_OLD_M2",
        "muT3GeneralS1S2": "WRONG_OLD_MU",
    }
    mass = {
        "mHSq": "new_higgs", "mS1Sq": "new_m1",
        "mS2Sq": "new_m2", "muT3GeneralS1S2": "new_mu",
    }
    return layout, old, mass


def test_mass_betas_replace_all_three_old_betas(example_betas):
    layout, old, mass = example_betas
    merged = coupled.merge_beta_reports(old, mass, layout)
    assert set(merged) == set(layout) | {"mHSq", "muT3GeneralS1S2"}
    assert merged["gY"] == old["gY"]
    assert merged["lambdaH"] == old["lambdaH"]
    for key in closure.NAMES:
        assert merged[key] == mass[key]
    assert not any("WRONG_OLD" in s for s in merged.values())


@pytest.mark.parametrize("drop_from", ["old", "mass"])
def test_fail_closed_for_incomplete_beta_inventory(example_betas, drop_from):
    layout, old, mass = example_betas
    if drop_from == "old":
        old.pop("gY")
    else:
        mass.pop("mHSq")
    with pytest.raises(ValueError):
        coupled.merge_beta_reports(old, mass, layout)


def test_five_real_derivative_uses_one_loop_factor_once(monkeypatch):
    # Depend only on the adapter's vector conventions, not Stage33 or Wolfram.
    class Layout:
        size = 3
        names = ("gY", "mS1Sq", "mS2Sq")
    class Baseline:
        pass
    class State:
        baseline = Baseline()
        extra = {"newQuartic": 0.12 + 0.02j}

    def unpack(y, layout, registry, *, uv, mu_gev):
        assert not uv and y.shape == (3,) and mu_gev == 1e5
        return State()

    monkeypatch.setattr(coupled, "unpack_general_scalar_state", unpack)
    monkeypatch.setattr(coupled, "intermediate_scalar_state_environment",
                        lambda state: {"gY": 0.5})
    evaluations = []
    def evaluate(expr, env):
        assert env["mHSq"] == -7800
        assert env["muT3GeneralS1S2"] == 4200 + 3100j
        assert env["newQuartic"] == 0.12 + 0.02j
        evaluations.append(expr)
        return {"g": 16., "h": 48., "a": 32.,
                "b": 64., "mu": 8.+12.j}[expr]
    monkeypatch.setattr(coupled, "evaluate_inputform_expression", evaluate)
    def base_derivative(values, layout, *, divide_by_loop_factor):
        assert divide_by_loop_factor is True
        assert set(values) == set(layout.names)
        return np.array([values[k] for k in layout.names]) / coupled.LOOP
    monkeypatch.setattr(coupled, "beta_values_to_derivative", base_derivative)
    y = np.array([0.5, 1.44e6, 1.96e6, 4200, 3100, -7800.])
    beta = {"gY": "g", "mHSq": "h", "mS1Sq": "a",
            "mS2Sq": "b", "muT3GeneralS1S2": "mu"}
    dy = coupled.beta_derivative(y, mu_gev=1e5, layout=Layout(),
                                 registry={}, betas=beta)
    assert evaluations == list(beta.values())
    np.testing.assert_allclose(dy * coupled.LOOP, [16.,32.,64.,8.,12.,48.],
                               rtol=0, atol=1e-12)


@pytest.mark.parametrize("model,filename", [
    ("B_-1", "B_m1_eft1_full_mass_v1.json"),
    ("C_-1", "C_m1_eft1_full_mass_v1.json"),
])
def test_existing_live_mass_betas_are_evaluable(model, filename):
    path = coupled.MASS_DIR / filename
    if not path.is_file():
        pytest.skip(f"live RGBeta mass output not present locally: {path}")
    report = closure.check_payload(closure.read_payload(path), model, "EFT1")
    env = closure.benchmark_environment(model, "EFT1")
    values = closure.evaluate_mass_vector(report, env)
    assert values.shape == (5,)
    assert np.isfinite(values).all()
    # Verify Higgs feedback into both BSM masses and vice versa at the same
    # native benchmark point by changing only the quadratic masses.
    delta_h = 1234.5
    env_modified = dict(env, mHSq=env["mHSq"]+delta_h)
    values_modified = closure.evaluate_mass_vector(report, env_modified)
    np.testing.assert_allclose(values_modified[1]-values[1],
                               4*env["lambdaH1"]*delta_h,
                               rtol=1e-9, atol=1e-6)
    np.testing.assert_allclose(values_modified[2]-values[2],
                               4*env["lambdaH2"]*delta_h,
                               rtol=1e-9, atol=1e-6)
