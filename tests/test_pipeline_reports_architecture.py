"""Architecture contracts for pipeline-level report generation.

The central ``pipeline.py`` should show that reports are generated after the
physics calculation, but detailed LaTeX/report orchestration belongs under
``Reports/``.  These checks protect the existing report and aggregate-output
contract while keeping presentation details out of the physics backbone.
"""

from __future__ import annotations

import ast
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
PIPELINE_PATH = PROJECT_ROOT / "pipeline.py"
REPORTS_PATH = PROJECT_ROOT / "Reports" / "PipelineReports.py"


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
    assert len(matches) == 1, f"expected exactly one {name}() in {path.name}"
    return matches[0]


def test_pipeline_imports_only_high_level_report_orchestration() -> None:
    imports = [
        node
        for node in _tree(PIPELINE_PATH).body
        if isinstance(node, ast.ImportFrom)
        and node.module == "Reports.PipelineReports"
    ]

    assert len(imports) == 1
    assert {alias.name for alias in imports[0].names} == {
        "finish_pipeline_run",
        "print_pipeline_summary",
    }

    source = _source(PIPELINE_PATH)
    assert "Reports.ReportGeneration" not in source
    assert "Reports.RGEComparison" not in source
    assert "Reports.GroupFactorReports" not in source


def test_pipeline_retains_compatibility_report_entry_points() -> None:
    print_wrapper = _function(PIPELINE_PATH, "print_summary")
    finish_wrapper = _function(PIPELINE_PATH, "finish_runs")

    print_source = ast.get_source_segment(_source(PIPELINE_PATH), print_wrapper)
    finish_source = ast.get_source_segment(_source(PIPELINE_PATH), finish_wrapper)
    assert print_source is not None
    assert finish_source is not None

    assert "print_pipeline_summary" in print_source
    assert "finish_pipeline_run" in finish_source


def test_report_module_preserves_aggregate_json_contract() -> None:
    source = _source(REPORTS_PATH)
    constants = {
        node.value
        for node in ast.walk(_tree(REPORTS_PATH))
        if isinstance(node, ast.Constant) and isinstance(node.value, str)
    }

    assert "t3_model_comparison.json" in constants
    assert "T3 MODEL SUMMARY" in source
    assert "BuildStatus" in constants
    assert "MatchingStatus" in constants
    assert "UVRGEStatus" in constants
    assert "json.dumps" in source
    assert "[record.summary for record in records]" in source


def test_report_generation_order_is_preserved() -> None:
    node = _function(REPORTS_PATH, "finish_pipeline_run")
    source = ast.get_source_segment(_source(REPORTS_PATH), node)
    assert source is not None

    ordered_calls = (
        "print_pipeline_summary",
        "write_bsm_uv_field_table",
        "write_bsm_matched_field_table",
        "write_c5_coefficient_report",
        "write_and_compile_rge_comparison",
        "write_and_compile_eft1_rge_comparison",
        "write_and_compile_final_eft_rge_comparison",
        "write_and_compile_stage_group_factor_reports",
    )
    positions = [source.index(name) for name in ordered_calls]
    assert positions == sorted(positions)


def test_report_summary_uses_generated_stage_content_not_hard_coded_eft1() -> None:
    source = _source(REPORTS_PATH)
    summary = _function(REPORTS_PATH, "_print_generated_report_summary")
    summary_source = ast.get_source_segment(source, summary)
    assert summary_source is not None

    assert "Stage-aware reports:" in summary_source
    assert "_unique_eft_stage_labels(records)" in summary_source
    assert "Lagrangian/{stage_label}" in summary_source
    assert "intermediate_rge_tex is not None" in summary_source
    assert "group_factor_tex" in summary_source

    # Historical report writers may still emit this compatibility path, but
    # the terminal report summary must not invent it for studies where the
    # corresponding physical stage does not exist.
    assert '"Lagrangian/EFT_1_after_F"' not in summary_source
    assert '"RGE/EFT_1_after_F"' not in summary_source


def test_main_still_finishes_only_after_all_physics_stages() -> None:
    node = _function(PIPELINE_PATH, "main")
    source = ast.get_source_segment(_source(PIPELINE_PATH), node)
    assert source is not None

    intermediate = source.index("_run_intermediate_eft_stages")
    low_energy = source.index("_run_low_energy_neutrino_stages")
    reports = source.index("finish_runs")

    assert intermediate < low_energy < reports
