"""Architecture contracts for the Python/Wolfram Lagrangian boundary."""

from __future__ import annotations

import ast
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RUNNER = PROJECT_ROOT / "Lagrangian" / "Runner.py"
VALIDATION = PROJECT_ROOT / "Lagrangian" / "ModelValidation.py"
WOLFRAM = PROJECT_ROOT / "Lagrangian" / "WolframRunner.py"
RESULTS = PROJECT_ROOT / "Lagrangian" / "MatchingResults.py"


def _source(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def test_runner_remains_the_stable_public_boundary() -> None:
    source = _source(RUNNER)

    for name in (
        "def run_model(",
        "def validate_dimensions(",
        "def validate_shared_dimensions(",
        "def obtain_class_dimensions(",
    ):
        assert name in source

    assert "Lagrangian.ModelValidation" in source
    assert "Lagrangian.WolframRunner" in source
    assert "Lagrangian.MatchingResults" in source


def test_external_process_execution_is_isolated_in_wolfram_runner() -> None:
    runner_source = _source(RUNNER)
    wolfram_source = _source(WOLFRAM)

    assert "subprocess" not in runner_source
    assert "shutil.rmtree" not in runner_source
    assert "subprocess.Popen" in wolfram_source
    assert "wolframscript" in wolfram_source
    assert "RunModel.wl" in wolfram_source


def test_model_validation_has_no_process_or_result_io() -> None:
    source = _source(VALIDATION)

    assert "valid_t3_topology_dimensions" in source
    assert "t3_has_neutral_bsm_component" in source
    assert "valid_shared_scalar_topology_dimensions" in source

    assert "subprocess" not in source
    assert "comparison_summary.json" not in source
    assert "RunRecord" not in source


def test_matching_result_ingestion_is_separate_from_process_execution() -> None:
    source = _source(RESULTS)

    assert "comparison_summary.json" in source
    assert "RunRecord" in source
    assert "SequentialStageExportFailed" in source
    assert "subprocess" not in source


def test_runner_contains_no_matching_summary_json_logic() -> None:
    source = _source(RUNNER)

    assert "json.loads" not in source
    assert "comparison_summary.json" not in source
    assert "SequentialStageExportFailed" not in source
