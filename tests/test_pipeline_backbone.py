"""Architecture regression checks for the central pipeline backbone.

``pipeline.py`` is intentionally the project map: it owns ordering and data
flow, while CLI parsing, model selection, detailed RGEs, low-energy physics,
validation, and reports live in specialised modules.
"""

from __future__ import annotations

import ast
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
PIPELINE_PATH = PROJECT_ROOT / "pipeline.py"


def _pipeline_source() -> str:
    return PIPELINE_PATH.read_text(encoding="utf-8")


def _pipeline_tree() -> ast.Module:
    return ast.parse(_pipeline_source())


def _function(tree: ast.Module, name: str) -> ast.FunctionDef:
    matches = [
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name == name
    ]
    assert len(matches) == 1, f"expected exactly one {name}()"
    return matches[0]


def _called_names(node: ast.AST) -> set[str]:
    names: set[str] = set()
    for child in ast.walk(node):
        if isinstance(child, ast.Call) and isinstance(child.func, ast.Name):
            names.add(child.func.id)
    return names


def test_pipeline_remains_compact_central_backbone() -> None:
    source = _pipeline_source()
    assert len(source.splitlines()) < 500
    assert "Central backbone for the T3 neutrino-mass calculation" in source
    assert "this file owns their order" in source


def test_pipeline_imports_configuration_helpers_instead_of_implementing_them() -> None:
    tree = _pipeline_tree()
    imports = [
        node
        for node in tree.body
        if isinstance(node, ast.ImportFrom) and node.module == "common.PipelineCLI"
    ]

    assert len(imports) == 1
    assert {(item.name, item.asname) for item in imports[0].names} == {
        ("build_argument_parser", "_build_argument_parser"),
        ("resolve_model_mode", "_resolve_model_mode"),
        ("resolve_pipeline_plan", "_resolve_pipeline_plan"),
        ("study_name", "_study_name"),
    }

    function_names = {
        node.name for node in tree.body if isinstance(node, ast.FunctionDef)
    }
    assert "_build_argument_parser" not in function_names
    assert "_resolve_model_mode" not in function_names
    assert "_resolve_pipeline_plan" not in function_names
    assert "_study_name" not in function_names
    assert "_build_run_records" not in function_names
    assert "_select_study_models" not in function_names
    assert "_build_and_match_t3_models" not in function_names
    assert "print_summary" not in function_names
    assert "finish_runs" not in function_names


def test_main_keeps_complete_calculation_order_explicit() -> None:
    tree = _pipeline_tree()
    main = _function(tree, "main")
    source = ast.get_source_segment(_pipeline_source(), main)
    assert source is not None

    ordered_calls = (
        "_build_argument_parser",
        "_resolve_model_mode",
        "_resolve_pipeline_plan",
        "_select_study_models",
        "_build_and_match_t3_models",
        "_attach_threshold_metadata",
        "organise_matched_weinberg_coefficient",
        "_run_intermediate_eft_stages",
        "_run_low_energy_neutrino_stages",
        "finish_pipeline_run",
    )
    positions = [source.index(name) for name in ordered_calls]
    assert positions == sorted(positions)


def test_main_uses_one_pipeline_plan_for_threshold_configuration() -> None:
    tree = _pipeline_tree()
    main = _function(tree, "main")
    called = _called_names(main)

    assert "_resolve_pipeline_plan" in called
    assert "_select_study_models" in called
    assert "_build_and_match_t3_models" in called
    assert "_attach_threshold_metadata" in called
    assert "_run_intermediate_eft_stages" in called

    source = ast.get_source_segment(_pipeline_source(), main)
    assert source is not None
    assert "pipeline_plan=pipeline_plan" in source


def test_metadata_attachment_delegates_to_pipeline_plan() -> None:
    tree = _pipeline_tree()
    attach = _function(tree, "_attach_threshold_metadata")
    source = ast.get_source_segment(_pipeline_source(), attach)
    assert source is not None

    assert "pipeline_plan.summary_metadata()" in source
    assert "pipeline_plan.build_stage_records(record.output_dir)" in source


def test_intermediate_dispatch_uses_physical_running_intervals() -> None:
    tree = _pipeline_tree()
    runner = _function(tree, "_run_intermediate_eft_stages")
    source = ast.get_source_segment(_pipeline_source(), runner)
    assert source is not None

    assert "pipeline_plan.running_intervals" in source
    assert "_stage_for_running_interval" in source
    assert "run_intermediate_eft_interval" in source
    assert "run_intermediate_eft_validation" in source


def test_pipeline_has_no_eft1_specific_dispatch() -> None:
    source = _pipeline_source()

    assert "from RGE.running.eft1" not in source
    assert "run_eft1_" not in source
    assert "record.first_eft_stage" not in source
    assert "_has_f_first_eft_stage" not in source
    assert "_is_f_then_s1_s2_plan" not in source
    assert "_is_f_then_scalar_threshold_plan" not in source


def test_full_study_orchestration_is_outside_the_single_study_backbone() -> None:
    source = _pipeline_source()
    tree = _pipeline_tree()

    imports = [
        node
        for node in tree.body
        if isinstance(node, ast.ImportFrom) and node.module == "studies.FullT3Study"
    ]
    assert len(imports) == 1
    assert {(item.name, item.asname) for item in imports[0].names} == {
        ("run_full_study", None),
    }

    function_names = {
        node.name for node in tree.body if isinstance(node, ast.FunctionDef)
    }
    assert "run_full_study" not in function_names
    assert "_forward_common_cli_args" not in function_names
    assert "import subprocess" not in source
    assert "import json" not in source
    assert "import time" not in source


def test_unsupported_intermediate_content_blocks_authoritative_low_energy_path() -> None:
    source = _pipeline_source()
    runner = _function(_pipeline_tree(), "_run_intermediate_eft_stages")
    runner_source = ast.get_source_segment(source, runner)
    assert runner_source is not None

    assert "IncompleteUnsupportedFieldContent" in runner_source
    assert "IntermediateEFTUnsupportedIntervals" in runner_source
    assert "if not outcome.attempted" in runner_source
    assert "physics_failed = True" in runner_source

    readiness = _function(_pipeline_tree(), "_ready_for_low_energy_neutrino_stages")
    readiness_source = ast.get_source_segment(source, readiness)
    assert readiness_source is not None
    assert "IntermediateEFTProductionStatus" in readiness_source
    assert "IncompleteUnsupportedFieldContent" in readiness_source
