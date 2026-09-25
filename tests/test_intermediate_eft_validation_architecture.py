"""Architecture checks for independent intermediate-EFT validation backends."""

from __future__ import annotations

import ast
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DISPATCH_PATH = PROJECT_ROOT / "validation" / "IntermediateEFTValidation.py"
BACKEND_PATH = (
    PROJECT_ROOT
    / "validation"
    / "backends"
    / "ScalarOnlyAfterFermionValidation.py"
)
PIPELINE_PATH = PROJECT_ROOT / "pipeline.py"
PRODUCTION_PATH = PROJECT_ROOT / "RGE" / "running" / "backends" / "ScalarOnlyAfterFermion.py"


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


def test_validation_dispatch_is_small_and_backend_agnostic() -> None:
    dispatch = _source(DISPATCH_PATH)

    assert len(dispatch.splitlines()) < 140
    assert "run_eft1_wilson_transport" not in dispatch
    assert "normalize_pole_rge_consistency" not in dispatch
    assert "build_final_weinberg_coefficient" not in dispatch


def test_validation_backend_owns_component_transport_regression() -> None:
    validation = _source(BACKEND_PATH)
    production = _source(PRODUCTION_PATH)

    assert "legacy_component_wilson_transport" in validation
    assert "run_eft1_wilson_transport" not in validation
    assert "run_eft1_wilson_transport" not in production
    assert "eft1_wilson_at_S_threshold.json" in validation
    assert "eft1_wilson_at_S_threshold.json" not in production


def test_validation_backend_owns_independent_pole_rge_check() -> None:
    validation = _source(BACKEND_PATH)
    production = _source(PRODUCTION_PATH)

    assert "normalize_pole_rge_consistency" in validation
    assert "normalize_pole_rge_consistency" not in production


def test_validation_dispatch_follows_production_field_content() -> None:
    node = _function(DISPATCH_PATH, "intermediate_validation_backend")
    source = ast.get_source_segment(_source(DISPATCH_PATH), node)
    assert source is not None

    assert "intermediate_running_backend(record, interval)" in source
    assert "SCALAR_ONLY_AFTER_FERMION_BACKEND" in source
    assert ".level" not in source


def test_validation_backend_run_is_checks_only() -> None:
    node = _function(BACKEND_PATH, "run")
    source = ast.get_source_segment(_source(BACKEND_PATH), node)
    assert source is not None

    assert "_run_component_transport_regression" in source
    assert "_validate_consistency_metadata" in source
    assert "build_final_weinberg_coefficient" not in source
    assert "rerun_threshold2_with_running" not in source


def test_pipeline_still_shows_production_then_validation_as_separate_steps() -> None:
    pipeline = _source(PIPELINE_PATH)
    node = _function(PIPELINE_PATH, "_run_intermediate_eft_stages")
    source = ast.get_source_segment(pipeline, node)
    assert source is not None

    production_pos = source.index("run_intermediate_eft_interval")
    validation_pos = source.index("run_intermediate_eft_validation")
    assert production_pos < validation_pos


def test_legacy_validation_report_keys_are_preserved_in_validation_backend() -> None:
    validation = _source(BACKEND_PATH)

    for key in (
        "EFT1WilsonTransportStatus",
        "EFT1WilsonTransportFile",
        "EFT1WilsonTransportEqualScaleCheck",
        "EFT1WilsonRGEWeinbergSubspaceValid",
        "EFT1FlavorSeedOneGenerationCheck",
        "EFT1FlavorTransportEqualScaleCheck",
    ):
        assert key in validation


def test_validation_backend_owns_interpretation_of_production_diagnostics() -> None:
    validation = _source(BACKEND_PATH)
    production = _source(PRODUCTION_PATH)

    for key in (
        "EFT1WilsonRGEWeinbergSubspaceValid",
        "EFT1WilsonRGEWeinbergBeta",
        "EFT1FlavorSeedOneGenerationCheck",
        "EFT1FlavorTransportEqualScaleCheck",
    ):
        assert key in validation
        assert key not in production

    for filename in (
        "eft1_wilson_rge.json",
        "eft1_wilson_flavor_seed.json",
        "eft1_wilson_flavor_at_S_threshold.json",
    ):
        assert filename in validation


def test_validation_backend_does_not_import_historical_eft1_package_directly() -> None:
    validation = _source(BACKEND_PATH)
    assert "from RGE.running.eft1" not in validation
    assert "LegacyEFT1Compatibility" in validation
