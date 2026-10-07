"""Tests for the comparison-only supervisor display builder."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from Numerical.plotting import BuildDisplayResults


def _diagnostic() -> dict:
    return {
        "status": "Success",
        "benchmark_search": {
            "model": {
                "model_key": "T3_dS1_1_dS2_3_dF_2_alpha_p0",
                "d_s1": 1,
                "d_s2": 3,
                "d_f": 2,
                "alpha": 0,
            }
        },
        "low_energy": {
            "masses_ev": [0.01, 0.02, 0.05],
            "delta_m21_sq_ev2": 7.5e-5,
            "delta_m31_sq_ev2": 2.5e-3,
            "takagi_residual": 1.0e-12,
        },
        "final_running": {
            "sin2_theta12": [0.31, 0.30],
            "sin2_theta13": [0.023, 0.022],
            "sin2_theta23": [0.57, 0.56],
        },
    }


def test_display_ignores_legacy_optimal_outputs(tmp_path: Path) -> None:
    data = tmp_path / "output" / "full"
    reports = tmp_path / "Reports" / "output" / "full"
    display = tmp_path / "Reports" / "output" / "display"

    scenario = "smallY_smallL"
    model = "T3_dS1_1_dS2_3_dF_2_alpha_p0"

    comparison_data = data / "comparison" / scenario / model / "data"
    comparison_data.mkdir(parents=True)
    (comparison_data / "running_diagnostics.json").write_text(
        json.dumps(_diagnostic()),
        encoding="utf-8",
    )

    # A stale successful optimal diagnostic must not enter the new display.
    optimal_data = data / "optimal" / model / "data"
    optimal_data.mkdir(parents=True)
    (optimal_data / "running_diagnostics.json").write_text(
        json.dumps(_diagnostic()),
        encoding="utf-8",
    )

    comparison_figures = reports / "comparison" / scenario / model / "figures"
    comparison_figures.mkdir(parents=True)
    for name in BuildDisplayResults.KEEP_FIGURES:
        (comparison_figures / name).write_bytes(b"figure")

    stale_optimal = reports / "optimal" / model / "figures"
    stale_optimal.mkdir(parents=True)
    (stale_optimal / "c5_running.png").write_bytes(b"stale")

    BuildDisplayResults.build_display(data, reports, display)

    assert not (display / "optimal").exists()

    summary = display / "model_summary.csv"
    assert summary.is_file()
    with summary.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 1
    assert rows[0]["scenario"] == scenario
    assert rows[0]["model"] == model

    copied = display / "comparison" / scenario / "models" / model
    assert (copied / "intermediate_direct_weinberg_running.png").is_file()
    assert (copied / "c5_threshold_contributions.png").is_file()
    assert (copied / "c5_running.png").is_file()

    readme = (display / "README.md").read_text(encoding="utf-8")
    assert "Legacy optimal/optimizer outputs are not read or copied." in readme
    assert "O(hbar^2)" in readme


def test_comparison_scenarios_match_full_study() -> None:
    assert BuildDisplayResults.SCENARIOS == {
        "smallY_smallL": (0.005, 0.1),
        "smallY_largeL": (0.005, 1.0),
        "largeY_smallL": (0.5, 0.1),
        "largeY_largeL": (0.5, 1.0),
    }
