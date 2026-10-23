"""Order counting and opt-in sidecar behaviour, with no fitted-physics claims."""
import numpy as np
import pytest

from Numerical.diagnostics.T3BHeavyLogAlignment import matrix_to_json
from Numerical.diagnostics.T3BExploratoryCoupledBenchmark import (
    calculate, fixed_order_log_from_common_boundary,
)
from Numerical.diagnostics.T3BPerturbativeAccuracy import (
    _solve_expansion, perturbative_accuracy,
)
from Numerical.orchestration.T3BCoupledLLSidecar import check_low_energy


def _samples():
    zero = matrix_to_json(np.zeros((3, 3)))
    return [
        dict(mu_gev=m, g2=0.6, lambdaH=0.13, lambda12=0.07,
             lambda12Adj=0.025, lambdaT3=[0.24, 0.03],
             ye=zero, yu=zero, yd=zero)
        for m in (1e5, 1e4, 1e3)
    ]


def _inputs():
    y1 = np.diag([0.1, 0.2, 0.3]).astype(complex)
    y2 = np.diag([0.5, 0.4, 0.6]).astype(complex)
    masses = np.array([9.8e4, 1.1e5, 1.2e5])
    return y1, y2, masses


def _coupled_data():
    samples = _samples()
    y1, y2, masses = _inputs()
    direct = fixed_order_log_from_common_boundary(samples, y1, y2, masses)
    baseline = direct + np.eye(3) * 1e-8
    ctl = {
        "status": "FixedOrderControlComplete", "mu_f_gev": 1e5,
        "mu_s_gev": 1e3, "intermediate_running": {"samples": samples},
        "fixed_order_c5_at_ms": matrix_to_json(baseline),
    }
    nrm = {
        "status": "UnfittedNormalizationPassed",
        "production_direct_c5": matrix_to_json(direct),
        "predicted_direct_c5": matrix_to_json(direct),
    }
    coupled = calculate(ctl, nrm, y1, y2, masses)
    return ctl, coupled, y1, y2, masses


def test_order3_is_first_reverse_influence_on_c5():
    s = _samples()
    a, b, m = _inputs()
    with_reverse = _solve_expansion(s, a, b, m, max_order=4, reverse=True)
    no_reverse = _solve_expansion(s, a, b, m, max_order=4, reverse=False)
    for n in (0, 1, 2):
        np.testing.assert_allclose(with_reverse[n], no_reverse[n], rtol=1e-10, atol=1e-23)
    assert np.linalg.norm(with_reverse[3] - no_reverse[3]) > 0


def test_expansion_reconstructs_coupled_solution_and_direct_anchor():
    ctl, coupled, y1, y2, masses = _coupled_data()
    result = perturbative_accuracy(ctl, coupled, y1, y2, masses)
    assert result["status"] == "PerturbativeLLBookkeepingComplete"
    assert result["one_loop_frozen_direct_production_residual"] < 1e-10
    assert result["resummed_c5_reconstruction_residual"] < 1e-7
    assert result["terms"]["order_2"]["norm_relative_to_order1"] > 0
    assert result["terms"]["order_3"]["reverse_contribution_relative_to_order1"] > 0
    assert result["production_replacement_allowed"] is False
    assert result["order2_complete_nlo_matching"] is False


def test_offdiagonal_constant_case_first_order_is_direct():
    s = _samples()
    a, b, masses = _inputs()
    # With frozen coupling and no diagonal beta, X[1] is exactly the production LL.
    first = _solve_expansion(s, a, b, masses, max_order=1, diagonal=False)[1]
    direct = fixed_order_log_from_common_boundary(s, a, b, masses)
    np.testing.assert_allclose(first, direct, rtol=1e-11, atol=1e-24)


def test_cannot_consider_incorrect_threshold_or_status():
    ctl, coupled, y1, y2, masses = _coupled_data()
    coupled["threshold_scheme"] = "endpoint_extrapolation"
    with pytest.raises(ValueError, match="common_scale_log_boundary"):
        perturbative_accuracy(ctl, coupled, y1, y2, masses)


def test_low_energy_linked_only_on_matching_c5_change():
    coupled = {"relative_difference_from_fixed_order": 0.025}
    low = {"status": "ExploratoryLowEnergyComparisonComplete",
           "threshold_scheme": "common_scale_log_boundary",
           "relative_c5_difference_at_ms": 0.025,
           "relative_c5_difference_at_low_scale": 0.0249,
           "mu_low_gev": 100.0}
    check = check_low_energy(coupled, low)
    assert check["low_energy_is_exploratory"] is True
    low["relative_c5_difference_at_ms"] = 0.04
    with pytest.raises(ValueError, match="inconsistent"):
        check_low_energy(coupled, low)


def test_max_order_validation():
    s = _samples()
    a, b, m = _inputs()
    with pytest.raises(ValueError, match="max_order"):
        _solve_expansion(s, a, b, m, max_order=0)


def test_fitted_scale_background_order_selection_on_shared_integration_mesh():
    """Stress regression: small complex neutrino Yukawas and changing couplings.

    Two separate adaptive integrations can disagree at the sub-1e-8 relative
    level even when Wilson-order 1 and 2 must be reverse-independent. The
    paired solver uses precisely the same numerical stages for both branches.
    """
    from Numerical.diagnostics.T3BPerturbativeAccuracy import _solve_expansion_pair

    mus = np.geomspace(1e5, 1e3, 9)
    eye = np.eye(3)
    ye = np.diag([3e-6, 6e-4, 1e-2]).astype(complex)
    yu = np.diag([5e-6, 3e-3, 0.8]).astype(complex)
    yd = np.diag([1e-5, 4e-4, 1.6e-2]).astype(complex)
    samples = []
    for j, mu in enumerate(mus):
        x = j / (len(mus) - 1)
        samples.append(dict(
            mu_gev=float(mu), g2=0.65 - 0.03*x,
            lambdaH=0.13 + 0.015*x, lambda12=0.08 + 0.02*x,
            lambda12Adj=0.025 - 0.01*x, lambdaT3=[0.10 + 0.012*x, 0],
            ye=matrix_to_json(ye*(1 + 0.03*x)),
            yu=matrix_to_json(yu*(1 + 0.025*x)),
            yd=matrix_to_json(yd*(1 + 0.01*x)),
        ))
    y1 = 1e-4 * np.array([[1.5, 0.2j, 0.4],
                          [0.3, -1.1, 0.8j],
                          [0.9j, 0.5, -0.6]], dtype=complex)
    y2 = 1e-4 * np.array([[0.2, 0.5, 0.3j],
                          [-0.7, 0.3j, 0.4],
                          [0.4j, -0.5, 0.9]], dtype=complex)
    masses = np.array([97529.9, 107033.9, 118015.6])
    with_reverse, no_reverse = _solve_expansion_pair(samples, y1, y2, masses)
    assert np.linalg.norm(with_reverse[1]) > 0
    assert np.linalg.norm(with_reverse[2]) > 0
    for n in (0, 1, 2):
        np.testing.assert_array_equal(with_reverse[n], no_reverse[n])
    assert np.linalg.norm(with_reverse[3] - no_reverse[3]) > 0
