"""Architecture contracts for T3 study/model selection."""

from __future__ import annotations

import ast
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
PIPELINE_PATH = PROJECT_ROOT / "pipeline.py"
STUDY_PATH = PROJECT_ROOT / "model" / "T3Study.py"


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


def test_pipeline_separates_model_selection_from_matching() -> None:
    source = _source(PIPELINE_PATH)
    imports = [
        node
        for node in _tree(PIPELINE_PATH).body
        if isinstance(node, ast.ImportFrom) and node.module == "model.T3Study"
    ]

    assert len(imports) == 1
    assert {(item.name, item.asname) for item in imports[0].names} == {
        ("select_study_models", "_select_study_models")
    }

    assert "common.T3Model" not in source
    assert "Lagrangian.Runner" not in source
    assert "HYPERCHARGE_COMPARISON" not in source
    assert "DIMENSION_COMPARISON" not in source


def test_study_module_preserves_comparison_scan_definitions() -> None:
    source = _source(STUDY_PATH)

    assert 'HYPERCHARGE_ALPHAS = tuple(range(-4, 3))' in source
    assert '("A", "B", "C", "D", "E")' in source
    assert 'label, candidates = "smoke", tuple(SMOKE)' in source
    assert 'label, candidates = "hypercharge comparison", HYPERCHARGE_COMPARISON' in source
    assert 'label, candidates = "dimension comparison", DIMENSION_COMPARISON' in source


def test_neutrality_filter_still_uses_t3_model_predicate() -> None:
    node = _function(STUDY_PATH, "_neutral_scan_points")
    source = ast.get_source_segment(_source(STUDY_PATH), node)
    assert source is not None

    assert "T3_CLASSES[model_class]" in source
    assert "t3_has_neutral_bsm_component" in source


def test_study_selection_no_longer_calls_lagrangian_runner() -> None:
    source = _source(STUDY_PATH)
    selection = _function(STUDY_PATH, "select_study_models")
    selection_source = ast.get_source_segment(source, selection)
    assert selection_source is not None

    assert "from Lagrangian.Runner" not in source
    assert "validate_shared_dimensions" not in selection_source
    assert "validate_dimensions" not in selection_source
    assert "obtain_class_dimensions" not in selection_source
    assert "T3ModelRequest" in selection_source


def test_legacy_build_run_records_wrapper_is_retained() -> None:
    wrapper = _function(STUDY_PATH, "build_run_records")
    source = ast.get_source_segment(_source(STUDY_PATH), wrapper)
    assert source is not None

    assert "select_study_models" in source
    assert "build_and_match_t3_models" in source
    assert "parser.error" in source
