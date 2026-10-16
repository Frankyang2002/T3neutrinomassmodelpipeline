"""Tests for the normalized comparison-only supervisor display builder."""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

from Numerical.plotting import BuildDisplayResults
from studies.FullT3Study import LAMBDA_T3_NORMALISATION


def _diagnostic() -> dict:
    return {
        "status": "Success",
        "benchmark_search": {"model": {
            "model_key": "T3_dS1_2_dS2_2_dF_3_alpha_m1",
            "d_s1": 2, "d_s2": 2, "d_f": 3, "alpha": -1,
        }},
        "low_energy": {
            "masses_ev": [.01, .02, .05],
            "delta_m21_sq_ev2": 7.5e-5,
            "delta_m31_sq_ev2": 2.5e-3,
            "takagi_residual": 1e-12,
        },
        "final_running": {
            "sin2_theta12": [.31, .30],
            "sin2_theta13": [.023, .022],
            "sin2_theta23": [.57, .56],
        },
    }


def test_display_records_reference_and_effective_lambda(tmp_path: Path, monkeypatch) -> None:
    """Exercise summary generation; plot contents are tested separately."""
    data = tmp_path / "output" / "full"
    reports = tmp_path / "Reports" / "output" / "full"
    display = tmp_path / "Reports" / "output" / "display"
    scenario = "smallY_smallL"
    model = "T3_dS1_2_dS2_2_dF_3_alpha_m1"
    source = data / "comparison" / scenario / model / "data"
    source.mkdir(parents=True)
    diagnostic = source / "running_diagnostics.json"
    diagnostic.write_text(json.dumps(_diagnostic()), encoding="utf-8")

    calls = []

    def fake_figures(path: Path, output: Path) -> Path:
        calls.append(Path(path))
        output.mkdir(parents=True, exist_ok=True)
        for name in BuildDisplayResults.KEEP_FIGURES:
            (output / name).write_bytes(b"figure")
        return output

    monkeypatch.setattr(BuildDisplayResults, "generate_running_result_figures", fake_figures)
    BuildDisplayResults.build_display(data, reports, display)
    assert calls == [diagnostic]
    with (display / "model_summary.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 1
    assert rows[0]["model_class"] == "C"
    assert math.isclose(float(rows[0]["lambdaT3_reference"]), 0.1)
    assert math.isclose(float(rows[0]["lambdaT3_effective"]), 0.1)
    assert math.isclose(float(rows[0]["lambdaT3_normalisation_factor"]), 1.0)
    readme = (display / "README.md").read_text(encoding="utf-8")
    assert "lambda_T3^ref = 0.1" in readme
    assert "lambda_T3^eff" in readme
    assert "physical matching group factors" in readme


def test_comparison_scenarios_match_full_study() -> None:
    assert BuildDisplayResults.SCENARIOS == {
        "smallY_smallL": (0.005, 0.1),
        "smallY_largeL": (0.005, 1.0),
        "largeY_smallL": (0.5, 0.1),
        "largeY_largeL": (0.5, 1.0),
    }
    assert BuildDisplayResults.LAMBDA_T3_NORMALISATION == LAMBDA_T3_NORMALISATION
