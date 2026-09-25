"""Architecture contracts for T3 UV construction and EFT matching."""

from __future__ import annotations

import ast
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
PIPELINE_PATH = PROJECT_ROOT / "pipeline.py"
MATCHING_PATH = PROJECT_ROOT / "Lagrangian" / "T3ModelMatching.py"


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


def test_pipeline_imports_explicit_matching_stage() -> None:
    imports = [
        node
        for node in _tree(PIPELINE_PATH).body
        if isinstance(node, ast.ImportFrom)
        and node.module == "Lagrangian.T3ModelMatching"
    ]
    assert len(imports) == 1
    assert {(item.name, item.asname) for item in imports[0].names} == {
        ("build_and_match_t3_models", "_build_and_match_t3_models")
    }


def test_matching_module_is_the_only_boundary_to_runner_entry_points() -> None:
    source = _source(MATCHING_PATH)
    assert "from Lagrangian.Runner import" in source
    assert "obtain_class_dimensions" in source
    assert "validate_dimensions" in source
    assert "validate_shared_dimensions" in source

    pipeline_source = _source(PIPELINE_PATH)
    assert "Lagrangian.Runner" not in pipeline_source


def test_matching_uses_pipeline_plan_only_at_compatibility_boundary() -> None:
    node = _function(MATCHING_PATH, "build_and_match_t3_model")
    source = ast.get_source_segment(_source(MATCHING_PATH), node)
    assert source is not None

    assert "threshold_plan = pipeline_plan.threshold_plan" in source
    assert "request.kind == \"class\"" in source
    assert "request.kind == \"shared_dimensions\"" in source
    assert "validate_dimensions" in source


def test_matching_keeps_rge_tensor_export_disabled() -> None:
    source = _source(MATCHING_PATH)
    assert "export_rge_tensors=False" in source
    assert "debug_reports=debug_reports" in source
    assert "force=force" in source
