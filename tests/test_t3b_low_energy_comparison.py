"""Regression tests for the low-energy fixed/coupled comparison."""
from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace
import importlib.util

import numpy as np
import pytest


try:
    from Numerical.diagnostics import T3BLowEnergyComparison as m
except ImportError:
    # Allow local pure-helper tests without the full user repository.
    source = (Path(__file__).resolve().parents[1] / "Numerical" /
              "diagnostics" / "T3BLowEnergyComparison.py")
    spec = importlib.util.spec_from_file_location("T3BLowEnergyComparison", source)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)


def encoded(matrix):
    return m._matrix_to_json(np.asarray(matrix, dtype=complex))


def example():
    fixed = np.array([[2, 0.3 + 0.1j, 0], [0.3 + 0.1j, 3, 0.2], [0, 0.2, 4]], complex)*1e-15
    direct = fixed * 0.4
    hard = fixed - direct
    evolved = direct * 0.98
    total = hard + evolved
    ctl = {"status": "FixedOrderControlComplete", "fixed_order_c5_at_ms": encoded(fixed)}
    cpl = {"status": "ExploratoryCoupledMatrixComplete",
           "threshold_scheme": "common_scale_log_boundary",
           "common_boundary_fixed_order_direct_residual": 1e-16,
           "production_direct_c5_at_ms": encoded(direct),
           "production_hard_c5_at_ms": encoded(hard),
           "coupled_eft1_c5_at_ms": encoded(evolved),
           "exploratory_total_c5_at_ms": encoded(total),
           "relative_difference_from_fixed_order": float(np.linalg.norm(total-fixed)/np.linalg.norm(fixed))}
    return ctl, cpl, fixed, total


def test_load_complex_symmetric_matrix():
    matrix = np.array([[1j, 2+3j, 4j], [2+3j, -2, 0], [4j, 0, 3]], complex)
    assert np.array_equal(m._matrix_from_json(encoded(matrix), name="unit"), matrix)
    m._require_majorana(matrix, "unit")


def test_reject_non_symmetric_and_incomplete_c5():
    a = np.eye(3, dtype=complex)
    a[0, 1] = 1
    with pytest.raises(ValueError, match="not symmetric"):
        m._require_majorana(a, "broken")
    with pytest.raises(ValueError, match="expected 3x3"):
        m._matrix_from_json({"real": [1], "imag": [1]}, name="broken")


def test_matching_decomposition_rejects_double_counting():
    ctl, cpl, fixed, total = example()
    a, b, checks = m._extract_threshold_matrices(ctl, cpl)
    assert np.array_equal(a, fixed)
    assert np.array_equal(b, total)
    assert checks["fixed_order_decomposition_residual"] < 1e-14
    assert checks["coupled_decomposition_residual"] < 1e-14
    cpl["exploratory_total_c5_at_ms"] = encoded(total + np.eye(3)*1e-16)
    with pytest.raises(ValueError, match="decomposition"):
        m._extract_threshold_matrices(ctl, cpl)


def test_reject_old_threshold_scheme():
    ctl, cpl, _, _ = example()
    cpl["threshold_scheme"] = "endpoint_extrapolation"
    with pytest.raises(ValueError, match="common_scale_log_boundary"):
        m._extract_threshold_matrices(ctl, cpl)


def test_pdg_pmns_mixing_angles():
    t12, t13, t23 = np.deg2rad([33.0, 8.5, 48.0])
    def r12(t):
        c, s = np.cos(t), np.sin(t)
        return np.array([[c,s,0],[-s,c,0],[0,0,1]])
    def r13(t):
        c, s = np.cos(t), np.sin(t)
        return np.array([[c,0,s],[0,1,0],[-s,0,c]])
    def r23(t):
        c, s = np.cos(t), np.sin(t)
        return np.array([[1,0,0],[0,c,s],[0,-s,c]])
    got = m._mixing_angles(np.abs(r23(t23)@r13(t13)@r12(t12)))
    np.testing.assert_allclose([got["theta12_deg"], got["theta13_deg"], got["theta23_deg"]],
                               [33, 8.5, 48], atol=1e-12)


def test_saved_ms_background_consistency():
    e = np.eye(3, dtype=complex)
    state = SimpleNamespace(mu_gev=1000., lambda12=0.1, lambda12Adj=0.2,
        lambdaT3=0.3+0.04j, sm=SimpleNamespace(g2=0.63, lambdaH=0.25,
        ye=0.01*e, yu=0.7*e, yd=0.05*e))
    ren = SimpleNamespace(eft1_threshold_state=state)
    sample = {"mu_gev":1000, "g2":0.63, "lambdaH":0.25,
              "lambda12":0.1,"lambda12Adj":0.2,"lambdaT3":[0.3,0.04],
              "ye":encoded(state.sm.ye),"yu":encoded(state.sm.yu),"yd":encoded(state.sm.yd)}
    ctl = {"intermediate_running":{"samples":[sample]}}
    assert max(m._check_saved_ms_background(ctl, ren).values()) < 1e-14
    sample["g2"] = 0.71
    with pytest.raises(ValueError, match="g2"):
        m._check_saved_ms_background(ctl, ren)


def test_comparison_deltas():
    a = {"masses_ev":[0.01,0.02,0.05], "delta_m21_sq_ev2":0.0003,
         "delta_m31_sq_ev2":0.0024,"delta_m32_sq_ev2":0.0021,
         "mixing_angles":{"theta12_deg":30,"theta13_deg":8,"theta23_deg":45}}
    b = {"masses_ev":[0.009,0.018,0.045], "delta_m21_sq_ev2":0.000243,
         "delta_m31_sq_ev2":0.001944,"delta_m32_sq_ev2":0.001701,
         "mixing_angles":{"theta12_deg":31,"theta13_deg":8.5,"theta23_deg":44}}
    result=m._comparisons(a,b)
    np.testing.assert_allclose(result["relative_mass_differences"], [-0.1]*3)
    assert result["mixing_angle_differences_deg"]["theta23_deg"] == -1


def test_real_sm_weinberg_running_is_linear_in_k():
    # This test runs against the actual repository package when installed.
    runner = pytest.importorskip("Numerical.running.SMWeinbergEvolution")
    matrix = np.diag([1e-6, 2e-6, 3e-6]).astype(complex)
    common = dict(gY=0.35, g2=0.64, g3=1.1, lambdaH=0.25,
        ye=np.diag([1e-5, 2e-3, 0.01]),
        yu=np.diag([1e-5, 0.005, 0.9]),
        yd=np.diag([1e-5, 0.001, 0.04]))
    s1 = runner.SMWeinbergInitialConditions(K=matrix, **common)
    s2 = runner.SMWeinbergInitialConditions(K=1.03*matrix, **common)
    a=runner.evolve_sm_weinberg(s1,1000.,100.)
    b=runner.evolve_sm_weinberg(s2,1000.,100.)
    assert a.solver_success and b.solver_success
    np.testing.assert_allclose(b.K,1.03*a.K,rtol=1e-7,atol=1e-12)
    np.testing.assert_allclose(a.ye,b.ye,rtol=1e-10,atol=1e-12)
