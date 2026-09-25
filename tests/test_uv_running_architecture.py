"""Architecture and output-contract checks for ultraviolet RGE running."""

from __future__ import annotations

import ast
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
PIPELINE_PATH = PROJECT_ROOT / "pipeline.py"
UV_PATH = PROJECT_ROOT / "RGE" / "running" / "UVRunning.py"


def _source(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _tree(path: Path) -> ast.Module:
    return ast.parse(_source(path))


def _function(path: Path, name: str) -> ast.FunctionDef:
    matches = [
        node
        for node in _tree(path).body
        if isinstance(node, ast.FunctionDef) and node.name == name
    ]
    assert len(matches) == 1
    return matches[0]


def test_pipeline_imports_uv_running_as_one_named_stage() -> None:
    imports = [
        node
        for node in _tree(PIPELINE_PATH).body
        if isinstance(node, ast.ImportFrom) and node.module == "RGE.running.UVRunning"
    ]

    assert len(imports) == 1
    assert {(item.name, item.asname) for item in imports[0].names} == {
        ("run_uv_rge", "run_uv_rgbeta_stage")
    }
    assert "RGBetaT3Running" not in _source(PIPELINE_PATH)


def test_uv_running_keeps_existing_machine_output_contract() -> None:
    source = _source(UV_PATH)
    constants = {
        node.value
        for node in ast.walk(_tree(UV_PATH))
        if isinstance(node, ast.Constant) and isinstance(node.value, str)
    }

    assert "uv_rgbeta_rge.json" in constants
    assert {"UVRGEStatus", "UVRGEError", "UVRGEFile", "UVRGEBetaCount"} <= constants
    assert "json.dumps(result.raw, indent=2)" in source


def test_uv_stage_still_runs_before_intermediate_intervals() -> None:
    node = _function(PIPELINE_PATH, "_run_intermediate_eft_stages")
    source = ast.get_source_segment(_source(PIPELINE_PATH), node)
    assert source is not None

    uv_position = source.index("run_uv_rgbeta_stage")
    interval_position = source.index("pipeline_plan.running_intervals")
    assert uv_position < interval_position
