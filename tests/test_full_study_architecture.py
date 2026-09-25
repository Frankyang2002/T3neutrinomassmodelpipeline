"""Regression tests for the ``--full`` study orchestrator."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from types import SimpleNamespace

from studies import FullT3Study


def _args(**overrides) -> argparse.Namespace:
    values = {
        "debug_reports": False,
        "force": False,
        "numerical": None,
        "threshold": None,
        "threshold_scale": None,
    }
    values.update(overrides)
    return argparse.Namespace(**values)


def test_forwarded_full_study_cli_contract_is_preserved() -> None:
    args = _args(
        debug_reports=True,
        force=True,
        numerical=Path("point.json"),
        threshold=[["F"], ["S1", "S2"]],
        threshold_scale=["MF", "MS"],
    )

    assert FullT3Study._forward_common_cli_args(args) == [
        "--debug-reports",
        "--force",
        "--numerical",
        "point.json",
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


def test_full_study_runs_the_three_historical_child_modes(
    tmp_path: Path,
    monkeypatch,
) -> None:
    project_root = tmp_path / "project"
    output_root = project_root / "output"
    report_root = project_root / "Reports" / "output"
    pipeline_script = project_root / "pipeline.py"
    project_root.mkdir()
    pipeline_script.write_text("# test pipeline\n", encoding="utf-8")

    monkeypatch.setattr(FullT3Study, "PROJECT_ROOT", project_root)
    monkeypatch.setattr(FullT3Study, "PIPELINE_SCRIPT", pipeline_script)
    monkeypatch.setattr(FullT3Study, "OUTPUT_DIR", output_root)
    monkeypatch.setattr(FullT3Study, "REPORT_OUTPUT_DIR", report_root)

    calls: list[list[str]] = []

    def fake_run(command, cwd):
        assert cwd == project_root
        calls.append(list(command))
        study = command[command.index("--study") + 1].split("/")[-1]
        aggregate = output_root / "full" / study / "t3_model_comparison.json"
        aggregate.parent.mkdir(parents=True, exist_ok=True)
        aggregate.write_text(
            json.dumps(
                [
                    {"BuildStatus": "Success", "MatchingStatus": "Success"},
                    {"BuildStatus": "Success", "MatchingStatus": "Failed"},
                ]
            ),
            encoding="utf-8",
        )
        return SimpleNamespace(returncode=0)

    clock = iter((0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0))
    monkeypatch.setattr(FullT3Study.subprocess, "run", fake_run)
    monkeypatch.setattr(FullT3Study.time, "time", lambda: next(clock))

    status = FullT3Study.run_full_study(_args())

    assert status == 0
    assert [call[2] for call in calls] == [
        "--smoke",
        "--hypercharge-comparison",
        "--dimension-comparison",
    ]
    assert [call[4] for call in calls] == [
        "full/smoke",
        "full/hypercharge",
        "full/dimensions",
    ]

    summary = json.loads(
        (output_root / "full" / "full_run_summary.json").read_text(encoding="utf-8")
    )
    assert summary["Status"] == "Success"
    assert [item["Study"] for item in summary["Studies"]] == [
        "smoke",
        "hypercharge",
        "dimensions",
    ]
    assert [item["SuccessfulModels"] for item in summary["Studies"]] == [1, 1, 1]
    assert [item["TotalModels"] for item in summary["Studies"]] == [2, 2, 2]
