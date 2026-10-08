"""Tests for the threshold-anchored Weinberg dashboard and optional components."""
from __future__ import annotations

import json
import math
import re
from pathlib import Path

import numpy as np
import pytest

from Numerical.plotting import InteractiveModelComparison as dashboard
from Numerical.plotting.WeinbergMatchedContinuation import (
    reconstruct_matched_continuation,
)


def _matrix_block(matrix: np.ndarray) -> dict:
    return {"real": matrix.real.tolist(), "imag": matrix.imag.tolist(),
            "abs": np.abs(matrix).tolist()}


def _diagnostics(*, complex_case: bool = False, include_complex_samples: bool = False,
                 with_direct: bool = True, model: str = "T3_A_alpha_m2") -> dict:
    mu = np.array([1e5, 1e4, 1e3])
    t = np.log(mu / 1e3) / np.log(1e5 / 1e3)
    hard = np.array([[3, 0.5, 0.25], [0.5, 2, -0.5], [0.25, -0.5, 1.5]], dtype=complex) * 1e-7
    direct_s = np.array([[-2, 1, -0.5], [1, -1, -0.25], [-0.5, -0.25, 0.5]], dtype=complex) * 1e-7
    direct_f = np.array([[.15, -.05, .02], [-.05, -.1, .01], [.02, .01, -.08]], dtype=complex) * 1e-7
    if complex_case:
        hard += 1j * np.array([[.3, .1, 0], [.1, .2, 0], [0, 0, .3]]) * 1e-7
        direct_s += 1j * np.array([[-.6, .2, .05], [.2, .4, 0], [.05, 0, -.2]]) * 1e-7
        direct_f += 1j * np.array([[-.12, 0, .03], [0, .02, 0], [.03, 0, .04]]) * 1e-7
    intermediate = direct_s[None, :, :] + t[:, None, None] * (direct_f - direct_s)[None, :, :]
    combined = hard + direct_s
    payload = {
        "status": "Success",
        "benchmark_search": {"model": {"model_key": model}},
        "scales_gev": {"mu_uv": 1e7, "mu_fermion_threshold": 1e5,
                       "mu_scalar_threshold": 1e3, "mu_low": 100.0},
        "final_running": {"mu_gev": [1e3, 316.0, 100.0],
                          "c5_abs": [np.abs(combined).tolist(),
                                     np.abs(.99 * combined).tolist(),
                                     np.abs(.98 * combined).tolist()]},
        "c5_threshold_contributions": {"hard": _matrix_block(hard),
                                       "direct_running": _matrix_block(direct_s),
                                       "combined": _matrix_block(combined)},
    }
    if with_direct:
        payload["intermediate_direct_weinberg"] = {
            "mu_gev": mu.tolist(), "delta_c5_abs": np.abs(intermediate).tolist(),
        }
        if include_complex_samples:
            payload["intermediate_direct_weinberg"].update({
                "delta_c5_real": intermediate.real.tolist(),
                "delta_c5_imag": intermediate.imag.tolist(),
            })
    return payload


def _write(root: Path, *, model: str = "T3_A_alpha_m2", **kwargs) -> None:
    path = root / "T3_dS1_1_dS2_3_dF_2_alpha_m2" / "data" / "running_diagnostics.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(_diagnostics(model=model, **kwargs)), encoding="utf-8")


def _render(tmp_path: Path) -> tuple[str, dict]:
    _write(tmp_path / "input")
    target = dashboard.write_dashboard(
        tmp_path / "input", tmp_path / "reports" / "interactive_comparison.html", "Test T3"
    )
    page = target.read_text(encoding="utf-8")
    match = re.search(r"const DATA=(.*?);\nconst META=", page, re.S)
    assert match is not None
    return page, json.loads(match.group(1))


def test_real_legacy_continuation_complex_sum_and_exact_matching() -> None:
    payload = _diagnostics()
    result = reconstruct_matched_continuation(payload)
    assert result.method == "validated_real_log_affine_reconstruction"
    assert result.maximum_relative_magnitude_residual < 1e-12
    np.testing.assert_allclose(result.combined, result.direct + result.hard)
    np.testing.assert_allclose(np.abs(result.combined[-1]), payload["final_running"]["c5_abs"][0])
    # Magnitudes must be taken AFTER the complex sum: they are not additive.
    assert abs(result.combined[-1, 0, 0]) < abs(result.hard[0, 0]) + abs(result.direct[-1, 0, 0])
    np.testing.assert_allclose(result.direct[0, 0, 0], .15e-7)


def test_explicit_complex_trajectory_is_used_without_phase_inference() -> None:
    payload = _diagnostics(complex_case=True, include_complex_samples=True)
    result = reconstruct_matched_continuation(payload)
    assert result.method == "saved_complex_intermediate_trajectory"
    assert abs(result.direct[0, 0, 0].imag) > 0
    np.testing.assert_allclose(np.abs(result.combined[-1]), payload["final_running"]["c5_abs"][0])


def test_legacy_complex_magnitudes_are_rejected_instead_of_fabricated() -> None:
    with pytest.raises(ValueError, match="Complex legacy"):
        reconstruct_matched_continuation(_diagnostics(complex_case=True))


def test_legacy_curve_must_be_log_affine_and_match_all_saved_points() -> None:
    payload = _diagnostics()
    payload["intermediate_direct_weinberg"]["delta_c5_abs"][1][0][0] *= 1.3
    with pytest.raises(ValueError, match="log-affine"):
        reconstruct_matched_continuation(payload)


def test_threshold_matching_is_checked_before_continuity() -> None:
    payload = _diagnostics()
    payload["c5_threshold_contributions"]["combined"]["real"][0][0] += 1e-7
    with pytest.raises(ValueError, match=r"Hard \+ direct"):
        reconstruct_matched_continuation(payload)


def test_model_collect_keeps_raw_data_and_adds_checked_continuation(tmp_path: Path) -> None:
    _write(tmp_path)
    data = dashboard.collect(tmp_path)
    model = "T3_A_alpha_m2"
    assert data["models"] == [model]
    assert data["overlays"]["c5_11"][model]["y"][0] == pytest.approx(.15e-7)
    assert data["matched"]["c5_11"][model]["y"][0] == pytest.approx(3.15e-7)
    assert data["matched"]["c5_11"][model]["y"][-1] == pytest.approx(1e-7)
    assert data["components"]["c5_11"][model]["hard"]["y"] == pytest.approx([3e-7] * 3)
    assert data["components"]["c5_11"][model]["direct"]["y"][-1] == pytest.approx(2e-7)
    assert data["continuation_methods"][model] == "validated_real_log_affine_reconstruction"
    matrix = np.abs(np.asarray(_matrix_block(np.array([[3, .5, .25], [.5, 2, -.5], [.25, -.5, 1.5]]) * 1e-7)["real"]))
    assert data["components"]["c5_norm"][model]["hard"]["y"][0] == pytest.approx(np.linalg.norm(matrix))


def test_unusable_data_does_not_get_a_fake_continuous_curve(tmp_path: Path) -> None:
    _write(tmp_path, complex_case=True)
    data = dashboard.collect(tmp_path)
    assert data["models"] == ["T3_A_alpha_m2"]
    assert "T3_A_alpha_m2" not in data["matched"].get("c5_11", {})
    assert "Complex legacy" in data["continuation_issues"]["T3_A_alpha_m2"]


def test_selector_checkbox_and_solid_vs_dotted_implementation(tmp_path: Path) -> None:
    page, data = _render(tmp_path)
    assert '<option value="c5_norm">' in page
    assert '<option value="c5_11">' in page
    assert '<option value="llss_kernel">' in page
    assert '<option value="direct_c5_11">' not in page
    assert 'id="showComponents" type="checkbox"' in page
    assert 'showComponents.onchange=draw' in page
    assert 'componentControl.hidden=!weinberg' in page
    assert 'componentLegend.hidden=!showParts' in page
    assert "drawLine(parts.hard,col,1.5,[2,5],.45)" in page
    assert "drawLine(parts.direct,col,1.7,[2,4],.75)" in page
    assert "drawLine(continuation,col,2.6)" in page
    assert "drawLine(final,col,2.6)" in page
    assert "drawLine(direct,col,2.3)" not in page
    assert "ctx.lineTo(X(scale),Y(atFinal))" not in page
    assert "modelControls.hidden=llss" in page
    assert "canvas.dataset.yMin=String(ymin)" in page
    assert data["matched"]["c5_11"]["T3_A_alpha_m2"]["y"][-1] == pytest.approx(1e-7)


def test_refresh_four_scenarios_using_saved_diagnostics(tmp_path: Path) -> None:
    raw = tmp_path / "raw"
    report = tmp_path / "reports"
    display = tmp_path / "display"
    for scenario in dashboard.COMPARISON_SCENARIOS:
        _write(raw / scenario)
    paths = dashboard.refresh_full_comparison(raw, report, display)
    assert len(paths) == 4
    for scenario in dashboard.COMPARISON_SCENARIOS:
        one = report / scenario / "interactive_comparison.html"
        two = display / scenario / "interactive_comparison.html"
        assert one.is_file() and two.is_file()
        assert one.read_bytes() == two.read_bytes()


def test_invalid_intermediate_matrix_not_fabricated() -> None:
    payload = _diagnostics()
    payload["intermediate_direct_weinberg"]["delta_c5_abs"] = [[[1, 2, 3]]]
    with pytest.raises(ValueError):
        reconstruct_matched_continuation(payload)


def test_old_comparison_frobenius_still_correct(tmp_path: Path) -> None:
    _write(tmp_path)
    data = dashboard.collect(tmp_path)
    model = "T3_A_alpha_m2"
    original = np.asarray(_diagnostics()["final_running"]["c5_abs"])
    assert data["quantities"]["c5_norm"][model]["y"][0] == pytest.approx(np.linalg.norm(original[0]))
    reconstructed = data["matched"]["c5_norm"][model]["y"]
    assert reconstructed[-1] == pytest.approx(np.linalg.norm(original[0]))
