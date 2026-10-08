"""Regression tests for the clean combined Weinberg/LLSS dashboard."""
from __future__ import annotations

import json
import math
import re
from pathlib import Path

import pytest

from Numerical.plotting import InteractiveModelComparison as dashboard


def _matrix(value: float) -> list[list[float]]:
    return [[value, value / 4, 0.0], [value / 4, value / 2, 0.0], [0.0, 0.0, value / 3]]


def _diagnostics(*, with_direct: bool = True, model: str = "T3_A_alpha_m2") -> dict:
    payload = {
        "status": "Success",
        "benchmark_search": {"model": {"model_key": model}},
        "scales_gev": {
            "mu_uv": 1e7,
            "mu_fermion_threshold": 1e5,
            "mu_scalar_threshold": 1e3,
            "mu_low": 100.0,
        },
        "final_running": {
            "mu_gev": [1000.0, 316.0, 100.0],
            "c5_abs": [_matrix(2e-7), _matrix(1.9e-7), _matrix(1.8e-7)],
        },
    }
    if with_direct:
        payload["intermediate_direct_weinberg"] = {
            "mu_gev": [1e5, 1e4, 1e3],
            "delta_c5_abs": [_matrix(0), _matrix(5e-8), _matrix(1e-7)],
        }
    return payload


def _write(root: Path, *, with_direct: bool = True, model: str = "T3_A_alpha_m2") -> None:
    path = root / "T3_dS1_1_dS2_3_dF_2_alpha_m2" / "data" / "running_diagnostics.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(_diagnostics(with_direct=with_direct, model=model)), encoding="utf-8")


def _render(tmp_path: Path, *, with_direct: bool = True) -> tuple[str, dict]:
    _write(tmp_path / "input", with_direct=with_direct)
    output = dashboard.write_dashboard(
        tmp_path / "input", tmp_path / "reports" / "interactive_comparison.html", "Test T3"
    )
    page = output.read_text(encoding="utf-8")
    match = re.search(r"const DATA=(.*?);\nconst META=", page, re.S)
    assert match is not None
    return page, json.loads(match.group(1))


def test_saved_direct_and_final_are_distinct_and_preserved(tmp_path: Path) -> None:
    _write(tmp_path)
    data = dashboard.collect(tmp_path)
    model = "T3_A_alpha_m2"
    assert data["models"] == [model]
    assert data["quantities"]["c5_11"][model]["y"] == [2e-7, 1.9e-7, 1.8e-7]
    assert data["overlays"]["c5_11"][model]["y"] == [0, 5e-8, 1e-7]
    assert data["quantities"]["direct_c5_11"][model]["x"] == [1e5, 1e4, 1e3]
    expected_norm = 2e-7 * math.sqrt(1 + 2 * .25 ** 2 + .5 ** 2 + (1 / 3) ** 2)
    assert data["quantities"]["c5_norm"][model]["y"][0] == pytest.approx(expected_norm)


def test_no_unavailable_direct_coefficient_is_fabricated(tmp_path: Path) -> None:
    _write(tmp_path, with_direct=False)
    data = dashboard.collect(tmp_path)
    assert "c5_11" in data["quantities"]
    assert "c5_11" not in data["overlays"]


def test_clean_selector_and_separate_llss_view(tmp_path: Path) -> None:
    page, data = _render(tmp_path)
    assert '<option value="c5_norm">' in page
    assert '<option value="c5_11">' in page
    assert '<option value="llss_kernel">' in page
    assert '<option value="direct_c5_11">' not in page
    assert 'const kernel=llss?kernelCurve():null' in page
    assert 'modelControls.hidden=llss' in page
    assert '[hidden]{display:none!important}' in page
    assert 'id="axisTitle"' in page
    assert 'NOT a numerically evaluated LLSS Wilson coefficient' in page
    assert data["overlays"]["c5_11"]["T3_A_alpha_m2"]["y"][-1] == 1e-7


def test_matching_jump_is_explicit_and_no_dashed_operator_lines(tmp_path: Path) -> None:
    page, _ = _render(tmp_path)
    assert 'const atDirect=matchedEndpoint(direct,scale)' not in page  # declarations inline with prefix
    assert 'matchedEndpoint(direct,scale),atFinal=matchedEndpoint(final,scale)' in page
    assert 'ctx.lineTo(X(scale),Y(atFinal))' in page
    assert 'drawLine(direct,col,2.3)' in page
    assert 'drawLine(final,col,2.3)' in page
    assert 'hard-matching jump' in page
    assert 'ctx.setLineDash([7,4])' not in page
    assert 'LLSS one-loop kernel η' in page
    assert 'right axis' not in page.lower()


def test_refresh_uses_saved_diagnostics_for_four_scenarios(tmp_path: Path) -> None:
    raw = tmp_path / "raw"
    report = tmp_path / "report"
    display = tmp_path / "display"
    for scenario in dashboard.COMPARISON_SCENARIOS:
        _write(raw / scenario)
    paths = dashboard.refresh_full_comparison(raw, report, display)
    assert len(paths) == 4
    for scenario in dashboard.COMPARISON_SCENARIOS:
        source = report / scenario / "interactive_comparison.html"
        copy = display / scenario / "interactive_comparison.html"
        assert source.is_file() and copy.is_file()
        assert source.read_bytes() == copy.read_bytes()


def test_bad_matrix_shape_rejected() -> None:
    data = _diagnostics()
    data["intermediate_direct_weinberg"]["delta_c5_abs"] = [[[1.0, 2.0, 3.0]]]
    with pytest.raises((IndexError, ValueError)):
        dashboard._extract(data, dashboard.QUANTITIES["direct_c5_11"])


def test_y_axis_has_auto_and_fixed_absolute_modes(tmp_path: Path) -> None:
    page, _ = _render(tmp_path)
    assert '<select id="yScaleMode">' in page
    assert '<option value="auto">' in page
    assert '<option value="fixed">' in page
    assert "const fixedY=yScaleMode.value==='fixed'" in page
    assert 'if(!fixedY&&!enabled[m])return' in page
    assert 'if(fixedY||(Math.log10(x)' in page
    assert 'canvas.dataset.yMin=String(ymin)' in page
    assert 'canvas.dataset.yMax=String(ymax)' in page
    assert 'yScaleMode.onchange=draw' in page
