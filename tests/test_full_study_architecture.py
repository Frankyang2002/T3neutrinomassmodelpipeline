"""Regression tests for the multiplicity-only comparison ``--full`` orchestrator."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from studies import FullT3Study


def _args(**overrides) -> argparse.Namespace:
    values = {
        "debug_reports": False, "force": False, "numerical": None,
        "threshold": None, "threshold_scale": None,
        "reset_numerical_configs": False,
    }
    values.update(overrides)
    return argparse.Namespace(**values)


def test_full_study_grid_and_comparison_scenarios_are_current() -> None:
    assert len(FullT3Study.MODELS) == 16
    assert {item[0] for item in FullT3Study.MODELS} == {"A", "B", "C", "D", "E"}
    assert FullT3Study.COMPARISON_SCENARIOS == {
        "smallY_smallL": (0.005, 0.1),
        "smallY_largeL": (0.005, 1.0),
        "largeY_smallL": (0.5, 0.1),
        "largeY_largeL": (0.5, 1.0),
    }
    assert FullT3Study.LAMBDA_T3_NORMALISATION == {
        "A": 1.0, "B": 1.0, "C": 1.0, "D": 1.0, "E": 0.5,
    }


def test_common_cli_args_forward_only_shared_child_options() -> None:
    args = _args(
        debug_reports=True, force=True, numerical=Path("point.json"),
        threshold=[["F"], ["S1", "S2"]], threshold_scale=["MF", "MS"],
    )
    assert FullT3Study._common_args(args) == [
        "--debug-reports", "--force", "--threshold", "F",
        "--threshold", "S1", "S2", "--threshold-scale", "MF",
        "--threshold-scale", "MS",
    ]


def test_lambda_t3_compensation_is_applied_only_at_initialization() -> None:
    template = {"base_state": {"ordinary": {
        "y1": {"real": [[0.0] * 3 for _ in range(3)],
               "imag": [[0.0] * 3 for _ in range(3)]},
        "y2": {"real": [[0.0] * 3 for _ in range(3)],
               "imag": [[0.0] * 3 for _ in range(3)]},
        "lambdaT3": {"real": 0.0, "imag": 0.0},
    }}, "benchmark_search": {}, "sensitivity": {}}
    cases = {
        "A": (1, 3, 2, 1.0),
        "B": (2, 2, 1, 1.0),
        "C": (2, 2, 3, 1.0),
        "D": (3, 1, 2, 1.0),
        "E": (3, 3, 2, 0.5),
    }
    for cls, (d1, d2, df, factor) in cases.items():
        config = FullT3Study._fixed_comparison_config(
            template, cls, d1, d2, df, 0, 0.005, 0.1
        )
        point = config["comparison_point"]
        assert math.isclose(
            config["base_state"]["ordinary"]["lambdaT3"]["real"],
            0.1 * factor,
        )
        assert point["lambdaT3_reference"] == 0.1
        assert math.isclose(point["lambdaT3_effective"], 0.1 * factor)
        assert math.isclose(point["normalisation_factor"], factor)
        assert "multiplicity" in point["normalisation_basis"]
        assert "unbroken-SU(2)" in point["normalisation_basis"]


def test_full_study_runs_only_common_parameter_comparisons(
    tmp_path: Path, monkeypatch
) -> None:
    project_root = tmp_path / "project"
    output_root = project_root / "output" / "full"
    report_root = project_root / "Reports" / "output" / "full"
    config_root = project_root / "configs" / "generated_model_sets"
    active_root = project_root / "configs" / "generated_models"
    pipeline_script = project_root / "pipeline.py"
    template_path = project_root / "configs" / "template.json"
    template_path.parent.mkdir(parents=True)
    pipeline_script.write_text("# test\n", encoding="utf-8")
    template_path.write_text(json.dumps({
        "base_state": {"ordinary": {
            "y1": {"real": [[.1, 0, 0]] * 3, "imag": [[0, 0, 0]] * 3},
            "y2": {"real": [[.1, 0, 0]] * 3, "imag": [[0, 0, 0]] * 3},
            "lambdaT3": {"real": .1, "imag": 0},
        }},
        "benchmark_search": {}, "sensitivity": {},
    }), encoding="utf-8")
    for name, value in [
        ("PROJECT_ROOT", project_root),
        ("PIPELINE_SCRIPT", pipeline_script),
        ("OUTPUT_ROOT", output_root),
        ("REPORT_ROOT", report_root),
        ("CONFIG_SET_ROOT", config_root),
        ("ACTIVE_CONFIG_DIR", active_root),
    ]:
        monkeypatch.setattr(FullT3Study, name, value)
    monkeypatch.setattr(
        FullT3Study, "MODELS",
        (("A", 1, 3, 2, 0), ("C", 2, 2, 3, -1), ("E", 3, 3, 2, -2)),
    )
    monkeypatch.setattr(
        FullT3Study, "COMPARISON_SCENARIOS",
        {"smallY_smallL": (0.005, 0.1), "largeY_largeL": (0.5, 1.0)},
    )
    calls = []
    monkeypatch.setattr(
        FullT3Study, "_run", lambda command: calls.append(list(command)) or 0
    )
    monkeypatch.setattr(FullT3Study, "_dashboard", lambda *_a, **_k: None)
    assert FullT3Study.run_full_study(_args(numerical=template_path)) == 0
    assert len(calls) == 6
    assert all("/optimal/" not in c[c.index("--study") + 1] for c in calls)
    configs = [
        json.loads(Path(c[c.index("--numerical") + 1]).read_text())
        for c in calls[:3]
    ]
    effective = [c["base_state"]["ordinary"]["lambdaT3"]["real"] for c in configs]
    assert all(math.isclose(actual, expected) for actual, expected in zip(
        effective, (0.1, 0.1, 0.05), strict=True
    ))
    assert all(c["benchmark_search"]["enabled"] is False for c in configs)
    summary = json.loads((output_root / "full_numerical_summary.json").read_text())
    assert summary["study_mode"] == "comparison_only"
    assert summary["optimizer_in_full_study"] is False
    normalisation = summary["comparison_definition"]["lambdaT3_normalisation"]
    assert normalisation["physical_group_factors_modified"] is False
    assert normalisation["factors_by_class"] == FullT3Study.LAMBDA_T3_NORMALISATION
