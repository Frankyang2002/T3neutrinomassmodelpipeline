"""Compatibility regression tests for the current presentation-only report API.

The former collect/build_within_model_comparisons interface has been retired.
The current API lives in RunningTrajectoryData and BuildComparisonReports.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from Numerical.plotting.RunningTrajectoryData import SCENARIOS, SavedRun, discover_saved_runs
from Numerical.plotting.WithinModelComparison import LABELS, generate_model_figures
from Numerical.orchestration.BuildComparisonReports import build_reports
from studies.FullT3Study import COMPARISON_SCENARIOS, LAMBDA_T3_NORMALISATION


def test_scenario_values_are_current() -> None:
    assert SCENARIOS == tuple(COMPARISON_SCENARIOS)
    assert COMPARISON_SCENARIOS == {
        "smallY_smallL": (0.005, 0.1),
        "smallY_largeL": (0.005, 1.0),
        "largeY_smallL": (0.5, 0.1),
        "largeY_largeL": (0.5, 1.0),
    }
    assert all(scenario in LABELS for scenario in SCENARIOS)


def test_multiplicity_only_benchmark_factors() -> None:
    assert LAMBDA_T3_NORMALISATION == {
        "A": 1.0, "B": 1.0, "C": 1.0, "D": 1.0, "E": 0.5,
    }


@pytest.mark.parametrize("model", [
    "T3_A_alpha_m2", "T3_B_alpha_m1", "T3_C_alpha_p1",
    "T3_D_alpha_p0", "T3_E_alpha_m2",
])
def test_legacy_class_keys_are_not_fabricated_as_current_model_keys(model: str, tmp_path: Path) -> None:
    """Current discovery accepts representation-encoded model directories only."""
    path = tmp_path / SCENARIOS[0] / model / "data" / "running_diagnostics.json"
    path.parent.mkdir(parents=True)
    path.write_text(json.dumps({"status": "Success"}), encoding="utf-8")
    assert discover_saved_runs(tmp_path) == {}


def test_no_raw_diagnostics_raises_instead_of_generating_fake_reports(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError, match="No saved full-comparison diagnostics"):
        build_reports(tmp_path / "raw", tmp_path / "reports", strict=True)


def test_four_scenarios_preserve_distinct_saved_runs(tmp_path: Path) -> None:
    """Use real directory layout without requiring generated physics fixtures."""
    model = "T3_dS1_3_dS2_3_dF_2_alpha_m2"
    root = tmp_path / "comparison"
    from runpy import run_path
    fixture_payload = run_path(str(Path(__file__).with_name("test_within_model_reports.py")))["fixture_payload"]
    for idx, scenario in enumerate(SCENARIOS):
        path = root / scenario / model / "T3_E_alpha_m2" / "data" / "running_diagnostics.json"
        path.parent.mkdir(parents=True)
        path.write_text(json.dumps(fixture_payload(idx + 1)), encoding="utf-8")
    runs = discover_saved_runs(root)[model]
    assert tuple(runs) == SCENARIOS
    assert len({str(run.path) for run in runs.values()}) == 4
    output = tmp_path / "reports"
    summary = build_reports(root, output, strict=True)
    assert (output / model / "interactive_comparison.html").is_file()
    assert (output / model / "figures" / "c5_comparison.png").is_file()
    assert set(summary["models"][model]["scenarios"]) == set(SCENARIOS)
