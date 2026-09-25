"""Structural checks for generic intermediate-EFT production dispatch."""

from __future__ import annotations

import ast
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DISPATCH_PATH = PROJECT_ROOT / "RGE" / "running" / "IntermediateEFTRunning.py"
BACKEND_PATH = (
    PROJECT_ROOT
    / "RGE"
    / "running"
    / "backends"
    / "ScalarOnlyAfterFermion.py"
)


def _source(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _function(path: Path, name: str) -> ast.FunctionDef:
    tree = ast.parse(_source(path))
    matches = [
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name == name
    ]
    assert len(matches) == 1, f"expected exactly one {name}() in {path}"
    return matches[0]


def test_dispatch_module_is_small_and_contains_no_detailed_eft1_physics() -> None:
    dispatch = _source(DISPATCH_PATH)

    assert len(dispatch.splitlines()) < 160
    assert "from RGE.running.eft1" not in dispatch
    assert "run_rgbeta_t3_eft1" not in dispatch
    assert "build_final_weinberg_coefficient" not in dispatch
    assert "rerun_threshold2_with_running" not in dispatch


def test_public_runner_dispatches_by_physical_backend() -> None:
    node = _function(DISPATCH_PATH, "run_intermediate_eft_interval")
    source = ast.get_source_segment(_source(DISPATCH_PATH), node)
    assert source is not None

    assert "intermediate_running_backend(record, interval)" in source
    assert "run_scalar_only_after_fermion" in source
    assert "stage.level" not in source
    assert "first_eft_stage" not in source
    assert '"NotImplementedForFieldContent"' in source


def test_scalar_only_support_condition_lives_with_its_backend() -> None:
    dispatch = _source(DISPATCH_PATH)
    backend = _source(BACKEND_PATH)
    node = _function(BACKEND_PATH, "supports")
    support_source = ast.get_source_segment(backend, node)
    assert support_source is not None

    assert "_physical_scalar_fields" not in dispatch
    assert "record.physical_scalar_fields" in support_source
    assert 'frozenset({"S"}' not in support_source
    assert "interval.active_heavy_fields" in support_source
    assert "interval.entered_by.fields_to_integrate" in support_source
    assert "interval.exited_by.fields_to_integrate" in support_source
    assert "interval.exited_by.after.is_fully_decoupled" in support_source
    assert ".level" not in support_source


def test_verified_backend_keeps_existing_output_contract() -> None:
    backend = _source(BACKEND_PATH)

    for key in (
        "EFT1RenormalisableRGEStatus",
        "EFT1WilsonRGEStatus",
        "EFT1FullFlavorBridgeStatus",
        "WeinbergCoefficientFile",
        "FinalWeinbergCoefficientFile",
    ):
        assert key in backend

    for filename in (
        "eft1_rgbeta_rge.json",
        "eft1_wilson_rge.json",
        "eft1_wilson_flavor_seed.json",
        "final_weinberg_coefficient.json",
    ):
        assert filename in backend


def test_backend_run_function_reads_like_the_physical_sequence() -> None:
    backend = _source(BACKEND_PATH)
    node = _function(BACKEND_PATH, "run")
    source = ast.get_source_segment(backend, node)
    assert source is not None

    expected_order = [
        "_run_renormalisable_rge",
        "_run_dimension_five_wilson_rge",
        "_run_full_flavor_weinberg_matching",
    ]
    positions = [source.index(name) for name in expected_order]
    assert positions == sorted(positions)


def test_production_backend_contains_no_non_authoritative_component_transport() -> None:
    backend = _source(BACKEND_PATH)
    assert "run_eft1_wilson_transport" not in backend
    assert "eft1_wilson_at_S_threshold.json" not in backend
    assert "normalize_pole_rge_consistency" not in backend


def test_production_backend_uses_descriptive_intermediate_interfaces() -> None:
    backend = _source(BACKEND_PATH)

    assert "from RGE.running.eft1" not in backend
    assert "from RGE.running.intermediate" in backend
    assert "run_scalar_only_renormalisable_rge" in backend
    assert "run_scalar_only_dimension_five_wilson_rge" in backend
    assert "build_direct_weinberg_running" in backend
    assert "resume_scalar_threshold_with_running" in backend


def test_legacy_not_applicable_keys_remain_compatibility_only() -> None:
    dispatch = _source(DISPATCH_PATH)
    backend_node = _function(DISPATCH_PATH, "intermediate_running_backend")
    compatibility_node = _function(
        DISPATCH_PATH,
        "mark_intermediate_running_not_applicable",
    )
    backend_source = ast.get_source_segment(dispatch, backend_node)
    compatibility_source = ast.get_source_segment(dispatch, compatibility_node)
    assert backend_source is not None
    assert compatibility_source is not None

    assert "EFT1" not in backend_source
    assert "EFT1RenormalisableRGEStatus" in compatibility_source
    assert "EFT1WilsonRGEStatus" in compatibility_source
