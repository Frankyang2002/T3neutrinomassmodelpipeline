"""Synthetic fixture tests for the shared cross-/within-model HTML renderer."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from Numerical.plotting.RunningTrajectoryData import (
    SCENARIOS, discover_saved_runs, extract_trajectory,
)
from Numerical.orchestration.BuildComparisonReports import build_reports

MODEL = "T3_dS1_3_dS2_3_dF_2_alpha_m2"


def fixture_payload(factor: float) -> dict:
    mu_i = np.geomspace(1e5, 1e3, 5)
    mu_f = np.geomspace(1e3, 1e2, 5)
    zero = np.zeros((3, 3))
    hard = np.eye(3) * factor * 1e-12
    direct = np.eye(3) * factor * 2e-13
    fractions = np.log(mu_i / mu_i[0]) / np.log(mu_i[-1] / mu_i[0])
    running_direct = np.stack([direct * f for f in fractions])
    combined = hard + direct
    final_c5 = np.stack([combined * (1 + .02 * k) for k in range(5)])
    mat = lambda m: {"real": m.tolist(), "imag": zero.tolist(), "abs": np.abs(m).tolist()}
    return {"status": "Success", "scales_gev": {"mu_uv": 1e7,
            "mu_fermion_threshold": 1e5, "mu_scalar_threshold": 1e3, "mu_low": 1e2},
        "uv_running": {"mu_gev": [1e7, 1e6, 1e5], "y1_frobenius_norm": [.2,.21,.22],
                       "y2_frobenius_norm": [.1,.11,.12], "lambdaT3_abs": [.1,.11,.12]},
        "intermediate_direct_weinberg": {"mu_gev": mu_i.tolist(),
                                          "delta_c5_abs": np.abs(running_direct).tolist()},
        "c5_threshold_contributions": {"hard": mat(hard), "direct_running": mat(direct),
                                       "combined": mat(combined)},
        "final_running": {"mu_gev": mu_f.tolist(), "c5_abs": np.abs(final_c5).tolist(),
                          "delta_m21_sq_ev2": (factor*np.linspace(7e-5,8e-5,5)).tolist(),
                          "delta_m3l_sq_ev2": (factor*np.linspace(2e-3,2.5e-3,5)).tolist(),
                          "masses_ev": (factor*np.stack([np.linspace(.01,.011,5),
                          np.linspace(.02,.021,5),np.linspace(.05,.051,5)],axis=1)).tolist()}}


def test_extraction_norm_includes_off_diagonals():
    payload = fixture_payload(1)
    payload["final_running"]["c5_abs"][0][0][1] = 3e-12
    mu, norm = extract_trajectory(payload, "c5")
    assert mu.size == 5
    assert norm[0] == pytest.approx(np.sqrt(3*(1.2e-12)**2 + (3e-12)**2))


def test_four_scenarios_generate_with_shared_dashboard(tmp_path: Path):
    input_root, out = tmp_path / "raw", tmp_path / "reports"
    for index, scenario in enumerate(SCENARIOS):
        path = input_root / scenario / MODEL / "T3_E_alpha_m2" / "data" / "running_diagnostics.json"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(fixture_payload(index + 1)), encoding="utf-8")
    found = discover_saved_runs(input_root)
    assert len(found[MODEL]) == 4
    result = build_reports(input_root, out, strict=True)
    assert len(result["models"]) == 1
    assert (out / MODEL / "comparison_report.html").is_file()
    interactive = out / MODEL / "interactive_comparison.html"
    assert interactive.is_file()
    page = interactive.read_text(encoding="utf-8")
    assert "smallY_largeL" in page
    # These controls belong to InteractiveModelComparison's canonical renderer.
    assert 'id="zoomInt"' in page
    assert 'id="zoomUV"' in page
    assert 'id="zoomFinal"' in page
    assert 'id="showComponents"' in page
    assert "Show all scenarios" in page
    assert (out / MODEL / "figures" / "weinberg_eft_comparison.png").is_file()
    assert not result["models"][MODEL]["reconstruction_warnings"]


def test_strict_rejects_missing_scenario(tmp_path: Path):
    path = tmp_path / "raw" / SCENARIOS[0] / MODEL / "T3_E_alpha_m2" / "data" / "running_diagnostics.json"
    path.parent.mkdir(parents=True)
    path.write_text(json.dumps(fixture_payload(1)), encoding="utf-8")
    with pytest.raises(ValueError, match="missing scenarios"):
        build_reports(tmp_path / "raw", tmp_path / "reports", strict=True)
