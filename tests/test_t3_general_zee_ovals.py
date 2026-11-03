"""Read-only post-processing tests. No Wolfram, Matchete, or RGEs required."""
import csv
import json
import math
from pathlib import Path

import numpy as np
import pytest

from scripts import PlotT3GeneralZeeOvals as plot


def _fake_point(model, direction, gamma):
    k = {"T1": 1., "T2": 3., "T3": -2.}[direction]
    return {
        "model": model, "direction": direction, "gamma": gamma, "status": "Success",
        "chi2": 0.01 + 2 * (gamma - 1)**2,
        "delta_m21_sq_ev2": 7.537e-5 + k * (gamma-1) * 1e-5,
        "delta_m3l_sq_ev2": 2.511e-3 + k*(gamma-1)*1e-4,
        "sin2_theta12": 0.3088 + k*(gamma-1)*0.1,
        "sin2_theta13": 0.02248 + k*(gamma-1)*0.01,
        "sin2_theta23": 0.47 + k*(gamma-1)*0.05,
        "error": "",
    }


def _fixture(tmp_path, *, status="Complete", bad_direction=False):
    src = tmp_path / "source"
    src.mkdir()
    rows = [_fake_point("A_-2", d, g) for d in plot.DIRECTIONS
            for g in (0.99, 1., 1.01)]
    if bad_direction:
        rows.pop()
    with (src / "points.csv").open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    (src / "summary.json").write_text(json.dumps({
        "schema": "t3-general-zee-compensated-sensitivity-v1",
        "models": [{"model": "A_-2", "status": status, "failed_points": 0,
                    "baseline_chi2_saved": 0.01}],
        "gamma_grid": [0.99, 1., 1.01]
    }), encoding="utf-8")
    nu = tmp_path / "nufit.json"
    nu.write_text(json.dumps({
        "ordering": "NO", "dataset": "NuFIT 6.1 (2025), IC24 NO",
        "central_values": {k: _fake_point("A_-2", "T1", 1.)[k] for k in plot.OBS},
        "one_sigma_errors": {
            "delta_m21_sq_ev2": 9.7e-7, "delta_m3l_sq_ev2": 2.05e-5,
            "sin2_theta12": 0.00665, "sin2_theta13": 0.00057,
            "sin2_theta23": 0.0155,
        }
    }), encoding="utf-8")
    return src, nu


def test_juno_latest_measured_not_projection():
    assert plot.JUNO_2026["central_values"]["sin2_theta12"] == 0.3036
    assert plot.JUNO_2026["one_sigma_errors"]["sin2_theta12"] == 0.0064
    assert plot.JUNO_2026["central_values"]["delta_m21_sq_ev2"] == 7.388e-5
    assert plot.JUNO_2026["one_sigma_errors"]["delta_m21_sq_ev2"] == 7.8e-7
    assert "2609.37895" in plot.JUNO_2026["url"]


def test_joint_gaussian_radii():
    assert math.isclose(plot.JOINT_RADII["68% joint"]**2, 2.30, rel_tol=.003)
    assert math.isclose(plot.JOINT_RADII["99.73% joint"]**2, 11.83, rel_tol=.003)


def test_read_complete_scan(tmp_path):
    src, _ = _fixture(tmp_path)
    data, meta = plot.load_points(src / "points.csv", src / "summary.json", models=["A_-2"])
    assert len(data["A_-2"]["T1"]) == 3
    assert plot.centre(data["A_-2"])["gamma"] == 1.
    assert meta["A_-2"]["baseline_chi2_saved"] == 0.01


def test_centered_slope(tmp_path):
    src, _ = _fixture(tmp_path)
    data, _ = plot.load_points(src / "points.csv", src / "summary.json", models=["A_-2"])
    deriv, delta_g = plot.central_slope(data["A_-2"]["T2"], "sin2_theta12")
    expected = 0.3*(0.02)/(math.log(1.01) - math.log(0.99))
    assert deriv == pytest.approx(expected, rel=1e-12)
    assert delta_g == pytest.approx(0.02)


def test_model_covariance_is_about_unmoved_centre(tmp_path):
    src, _ = _fixture(tmp_path)
    data, _ = plot.load_points(src / "points.csv", src / "summary.json", models=["A_-2"])
    cov = plot.centered_covariance(data["A_-2"], "sin2_theta12", "delta_m21_sq_ev2")
    assert cov.shape == (2,2)
    assert np.allclose(cov, cov.T)
    assert np.linalg.eigvalsh(cov)[-1] > 0
    assert np.linalg.eigvalsh(cov)[0] >= -1e-17


def test_incomplete_status_rejected(tmp_path):
    src, _ = _fixture(tmp_path, status="PartialSuccess")
    with pytest.raises(ValueError, match="incomplete"):
        plot.load_points(src / "points.csv", src / "summary.json", models=["A_-2"])


def test_bad_direction_rejected(tmp_path):
    src, _ = _fixture(tmp_path, bad_direction=True)
    with pytest.raises(ValueError, match="Bad gamma grid|Inconsistent gamma grids"):
        plot.load_points(src / "points.csv", src / "summary.json", models=["A_-2"])


def test_full_end_to_end_one_model(tmp_path):
    src, nu = _fixture(tmp_path)
    out = tmp_path / "out"
    assert plot.main(["--source", str(src), "--nufit", str(nu),
                      "--output-dir", str(out), "--models", "A_-2",
                      "--no-per-model"]) == 0
    assert (out / "model_grid_sin2_theta12_vs_delta_m21_sq_ev2.png").stat().st_size > 1000
    assert (out / "solar_sensitivity_ranking.csv").is_file()
    assert (out / "baseline_nufit_juno_comparison.csv").is_file()
    prov = json.loads((out / "provenance.json").read_text())
    assert prov["plots_use_saved_points_only"] is True
    assert prov["models"] == ["A_-2"]
    assert "NOT" in prov["physics_scope"]
