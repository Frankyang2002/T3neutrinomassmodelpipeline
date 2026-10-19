"""Regression tests for the current 16-model frozen fitted-reference full study.

All runner tests use temporary inputs and mocked subprocess execution: they do
not launch Matchete, RGBeta, fitting, or the numerical T3 pipeline.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from studies import FullT3Study


REPRESENTATIONS = (
    (1, 3, 2, -4), (1, 3, 2, -2), (1, 3, 2, 0),
    (2, 2, 1, -3), (2, 2, 1, -1), (2, 2, 1, 1),
    (2, 2, 3, -3), (2, 2, 3, -1), (2, 2, 3, 1),
    (3, 1, 2, -2), (3, 1, 2, 0), (3, 1, 2, 2),
    (3, 3, 2, -4), (3, 3, 2, -2), (3, 3, 2, 0), (3, 3, 2, 2),
)


def _args(**overrides) -> argparse.Namespace:
    fields = {
        "debug_reports": False, "force": False, "numerical": None,
        "threshold": None, "threshold_scale": None,
        "reset_numerical_configs": False,
    }
    fields.update(overrides)
    return argparse.Namespace(**fields)


def _make_frozen_manifest(tmp_path: Path, monkeypatch) -> list[tuple[str, tuple[int, int, int, int], Path]]:
    root = tmp_path / "project"
    benchmark_root = root / "Numerical" / "benchmarks" / "fitted_reference"
    benchmark_root.mkdir(parents=True)
    rows = []
    expected = []
    for idx, dims in enumerate(REPRESENTATIONS):
        name = f"T3_fixture_{idx:02d}"
        relative = Path("Numerical/benchmarks/fitted_reference") / f"{name}.json"
        config = root / relative
        representation = dict(zip(("d_s1", "d_s2", "d_f", "alpha"), dims))
        payload = {
            "representation": representation,
            "benchmark_search": {"enabled": False},
            "sensitivity": {"enabled": False},
            "base_state": {"ordinary": {"lambdaT3": {"real": 0.1, "imag": 0.0}}},
        }
        config.write_text(json.dumps(payload), encoding="utf-8")
        rows.append({
            "model": name,
            "representation": {**representation, "shared_scalar": False},
            "config_path": relative.as_posix(),
            "sha256": hashlib.sha256(config.read_bytes()).hexdigest(),
        })
        expected.append((name, dims, config))
    manifest_path = benchmark_root / "benchmark_manifest.json"
    manifest_path.write_text(json.dumps({"models": rows}), encoding="utf-8")
    monkeypatch.setattr(FullT3Study, "PROJECT_ROOT", root)
    monkeypatch.setattr(FullT3Study, "PIPELINE_SCRIPT", root / "pipeline.py")
    monkeypatch.setattr(FullT3Study, "MANIFEST_PATH", manifest_path)
    monkeypatch.setattr(FullT3Study, "OUTPUT_ROOT", root / "output/full/fitted_reference")
    monkeypatch.setattr(FullT3Study, "REPORT_ROOT", root / "Reports/output/full/fitted_reference")
    return expected


def test_common_cli_args_and_fermion_first_defaults() -> None:
    assert FullT3Study._common_args(_args()) == [
        "--threshold", "F", "--threshold", "S1", "S2"
    ]
    args = _args(
        debug_reports=True, force=True,
        threshold=[["F"], ["S1", "S2"]], threshold_scale=["MF", "MS"],
    )
    assert FullT3Study._common_args(args) == [
        "--debug-reports", "--force", "--threshold", "F",
        "--threshold", "S1", "S2", "--threshold-scale", "MF",
        "--threshold-scale", "MS",
    ]


def test_manifest_loads_exactly_sixteen_frozen_models(tmp_path: Path, monkeypatch) -> None:
    expected = _make_frozen_manifest(tmp_path, monkeypatch)
    assert FullT3Study._load_models() == expected
    assert len(expected) == 16
    assert len(set(dims[:3] for _, dims, _ in expected)) == 5


def test_manifest_rejects_modified_frozen_configuration(tmp_path: Path, monkeypatch) -> None:
    expected = _make_frozen_manifest(tmp_path, monkeypatch)
    config = expected[0][2]
    payload = json.loads(config.read_text(encoding="utf-8"))
    payload["base_state"]["ordinary"]["lambdaT3"]["real"] = 0.2
    config.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="hash mismatch"):
        FullT3Study._load_models()


def test_manifest_rejects_enabled_benchmark_search(tmp_path: Path, monkeypatch) -> None:
    expected = _make_frozen_manifest(tmp_path, monkeypatch)
    config = expected[0][2]
    payload = json.loads(config.read_text(encoding="utf-8"))
    payload["benchmark_search"]["enabled"] = True
    config.write_text(json.dumps(payload), encoding="utf-8")
    manifest = json.loads(FullT3Study.MANIFEST_PATH.read_text(encoding="utf-8"))
    manifest["models"][0]["sha256"] = hashlib.sha256(config.read_bytes()).hexdigest()
    FullT3Study.MANIFEST_PATH.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(ValueError, match="enables a search"):
        FullT3Study._load_models()


def test_full_study_rejects_custom_numerical_config(tmp_path: Path, monkeypatch) -> None:
    _make_frozen_manifest(tmp_path, monkeypatch)
    with pytest.raises(ValueError, match="do not also supply --numerical"):
        FullT3Study.run_full_study(_args(numerical=Path("custom.json")))


def test_full_study_runs_one_frozen_point_per_model(tmp_path: Path, monkeypatch) -> None:
    expected = _make_frozen_manifest(tmp_path, monkeypatch)
    commands: list[list[str]] = []
    def fake_subprocess(command, *, cwd, check):
        assert cwd == FullT3Study.PROJECT_ROOT
        assert check is False
        commands.append(list(command))
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(FullT3Study.subprocess, "run", fake_subprocess)
    # Grouping is tested separately; never create grouping files during this test.
    from Numerical.orchestration import RunYukawaGrouping
    monkeypatch.setattr(RunYukawaGrouping, "run_grouping", lambda: {"groups": []})

    assert FullT3Study.run_full_study(_args()) == 0
    assert len(commands) == 16
    for command, (name, dims, config) in zip(commands, expected, strict=True):
        assert command[command.index("--dims") + 1:command.index("--dims") + 4] == [str(v) for v in dims[:3]]
        assert command[command.index("--alpha") + 1] == str(dims[3])
        assert command[command.index("--numerical") + 1] == str(config)
        assert command[command.index("--study") + 1] == f"full/fitted_reference/{name}"
        assert "--reset-numerical-configs" in command
    summary = json.loads((FullT3Study.OUTPUT_ROOT / "full_fitted_reference_summary.json").read_text(encoding="utf-8"))
    assert summary["study_mode"] == "fitted_reference"
    assert summary["model_count"] == 16
    assert summary["optimizer_in_full_study"] is False
    assert summary["status"] == "Success"
    assert len(summary["results"]) == 16
    assert all(row["config_sha256"] for row in summary["results"])
