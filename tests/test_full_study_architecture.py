"""Regression tests for the comparison-only ``--full`` orchestrator."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from studies import FullT3Study


def _args(**overrides) -> argparse.Namespace:
    values = {
        "debug_reports": False,
        "force": False,
        "numerical": None,
        "threshold": None,
        "threshold_scale": None,
        "reset_numerical_configs": False,
    }
    values.update(overrides)
    return argparse.Namespace(**values)


def test_full_study_grid_and_comparison_scenarios_are_current() -> None:
    assert len(FullT3Study.MODELS) == 16
    assert {item[0] for item in FullT3Study.MODELS} == {
        "A", "B", "C", "D", "E",
    }

    assert FullT3Study.COMPARISON_SCENARIOS == {
        "smallY_smallL": (0.005, 0.1),
        "smallY_largeL": (0.005, 1.0),
        "largeY_smallL": (0.5, 0.1),
        "largeY_largeL": (0.5, 1.0),
    }


def test_common_cli_args_forward_only_shared_child_options() -> None:
    args = _args(
        debug_reports=True,
        force=True,
        numerical=Path("point.json"),
        threshold=[["F"], ["S1", "S2"]],
        threshold_scale=["MF", "MS"],
    )

    assert FullT3Study._common_args(args) == [
        "--debug-reports",
        "--force",
        "--threshold",
        "F",
        "--threshold",
        "S1",
        "S2",
        "--threshold-scale",
        "MF",
        "--threshold-scale",
        "MS",
    ]


def test_full_study_runs_only_common_parameter_comparisons(
    tmp_path: Path,
    monkeypatch,
) -> None:
    project_root = tmp_path / "project"
    output_root = project_root / "output" / "full"
    report_root = project_root / "Reports" / "output" / "full"
    config_root = project_root / "configs" / "generated_model_sets"
    active_root = project_root / "configs" / "generated_models"
    pipeline_script = project_root / "pipeline.py"
    template_path = project_root / "configs" / "template.json"

    template_path.parent.mkdir(parents=True)
    pipeline_script.write_text("# test pipeline\n", encoding="utf-8")
    template_path.write_text(
        json.dumps(
            {
                "base_state": {
                    "ordinary": {
                        "y1": {
                            "real": [[0.1, 0.0, 0.0]] * 3,
                            "imag": [[0.0, 0.0, 0.0]] * 3,
                        },
                        "y2": {
                            "real": [[0.1, 0.0, 0.0]] * 3,
                            "imag": [[0.0, 0.0, 0.0]] * 3,
                        },
                        "lambdaT3": {
                            "real": 0.1,
                            "imag": 0.0,
                        },
                    }
                },
                "benchmark_search": {},
                "sensitivity": {},
            }
        ),
        encoding="utf-8",
    )

    monkeypatch.setattr(
        FullT3Study,
        "PROJECT_ROOT",
        project_root,
    )
    monkeypatch.setattr(
        FullT3Study,
        "PIPELINE_SCRIPT",
        pipeline_script,
    )
    monkeypatch.setattr(
        FullT3Study,
        "OUTPUT_ROOT",
        output_root,
    )
    monkeypatch.setattr(
        FullT3Study,
        "REPORT_ROOT",
        report_root,
    )
    monkeypatch.setattr(
        FullT3Study,
        "CONFIG_SET_ROOT",
        config_root,
    )
    monkeypatch.setattr(
        FullT3Study,
        "ACTIVE_CONFIG_DIR",
        active_root,
    )
    monkeypatch.setattr(
        FullT3Study,
        "MODELS",
        (
            ("A", 1, 3, 2, 0),
            ("B", 2, 2, 1, -1),
        ),
    )
    monkeypatch.setattr(
        FullT3Study,
        "COMPARISON_SCENARIOS",
        {
            "smallY_smallL": (0.005, 0.1),
            "largeY_largeL": (0.5, 1.0),
        },
    )

    calls: list[list[str]] = []
    monkeypatch.setattr(
        FullT3Study,
        "_run",
        lambda command: calls.append(list(command)) or 0,
    )
    monkeypatch.setattr(
        FullT3Study,
        "_dashboard",
        lambda *_args, **_kwargs: None,
    )

    status = FullT3Study.run_full_study(
        _args(numerical=template_path)
    )

    assert status == 0
    assert len(calls) == 4

    studies = [
        call[call.index("--study") + 1]
        for call in calls
    ]
    assert studies == [
        "full/comparison/smallY_smallL/T3_dS1_1_dS2_3_dF_2_alpha_p0",
        "full/comparison/smallY_smallL/T3_dS1_2_dS2_2_dF_1_alpha_m1",
        "full/comparison/largeY_largeL/T3_dS1_1_dS2_3_dF_2_alpha_p0",
        "full/comparison/largeY_largeL/T3_dS1_2_dS2_2_dF_1_alpha_m1",
    ]
    assert all("/optimal/" not in study for study in studies)

    for call in calls:
        config_path = Path(
            call[call.index("--numerical") + 1]
        )
        config = json.loads(
            config_path.read_text(encoding="utf-8")
        )
        assert config["benchmark_search"]["enabled"] is False
        assert config["sensitivity"]["enabled"] is False

    summary_path = output_root / "full_numerical_summary.json"
    summary = json.loads(
        summary_path.read_text(encoding="utf-8")
    )
    assert summary["status"] == "Success"
    assert summary["study_mode"] == "comparison_only"
    assert summary["optimizer_in_full_study"] is False
    assert summary["ordinary_model_count"] == 2
    assert summary["shared_scalar_models_included"] is False
    assert [
        item["mode"]
        for item in summary["results"]
    ] == [
        "comparison",
        "comparison",
    ]
