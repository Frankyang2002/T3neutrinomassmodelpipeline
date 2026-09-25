"""Architecture checks for descriptive intermediate-EFT running interfaces."""

from __future__ import annotations

import ast
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
INTERFACE_DIR = PROJECT_ROOT / "RGE" / "running" / "intermediate"
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


def _assert_thin_delegate(
    path: Path,
    public_name: str,
    compatibility_name: str,
) -> None:
    source = _source(path)
    node = _function(path, public_name)
    function_source = ast.get_source_segment(source, node)
    assert function_source is not None

    assert len(source.splitlines()) < 90
    assert compatibility_name in function_source
    assert "RGE.running.eft1" not in source
    assert "subprocess" not in source
    assert "json.loads" not in source
    assert "sympy" not in source


def test_scalar_only_running_interfaces_use_physical_names() -> None:
    expected = {
        "ScalarOnlyRenormalisableRunning.py": (
            "run_scalar_only_renormalisable_rge",
            "legacy_scalar_only_renormalisable_rge",
        ),
        "ScalarOnlyWilsonRunning.py": (
            "run_scalar_only_dimension_five_wilson_rge",
            "legacy_scalar_only_dimension_five_wilson_rge",
        ),
        "FlavorWilsonBoundary.py": (
            "export_full_flavor_wilson_boundary",
            "legacy_export_full_flavor_wilson_boundary",
        ),
        "DirectWeinbergRunning.py": (
            "build_direct_weinberg_running",
            "legacy_build_direct_weinberg_running",
        ),
        "ScalarThresholdMatching.py": (
            "resume_scalar_threshold_with_running",
            "legacy_resume_scalar_threshold_with_running",
        ),
    }

    for filename, (public_name, compatibility_name) in expected.items():
        path = INTERFACE_DIR / filename
        assert path.is_file()
        _assert_thin_delegate(path, public_name, compatibility_name)


def test_direct_weinberg_interface_exposes_insertion_export() -> None:
    path = INTERFACE_DIR / "DirectWeinbergRunning.py"
    _assert_thin_delegate(
        path,
        "export_direct_weinberg_insertion",
        "legacy_export_direct_weinberg_insertion",
    )


def test_production_backend_no_longer_imports_ordinal_eft1_modules() -> None:
    backend = _source(BACKEND_PATH)

    assert "from RGE.running.eft1" not in backend
    assert "from RGE.running.intermediate" in backend

    for public_name in (
        "run_scalar_only_renormalisable_rge",
        "run_scalar_only_dimension_five_wilson_rge",
        "export_full_flavor_wilson_boundary",
        "build_direct_weinberg_running",
        "export_direct_weinberg_insertion",
        "resume_scalar_threshold_with_running",
    ):
        assert public_name in backend


def test_legacy_implementation_names_are_absent_from_production_backend() -> None:
    backend = _source(BACKEND_PATH)

    for implementation_name in (
        "run_rgbeta_t3_eft1",
        "run_eft1_wilson_rge",
        "run_flavor_seed_export",
        "build_direct_weinberg_flavor_transport",
        "export_direct_weinberg_matchete",
        "rerun_threshold2_with_running",
    ):
        assert implementation_name not in backend


def test_existing_output_contract_remains_historical_for_compatibility() -> None:
    backend = _source(BACKEND_PATH)

    for key in (
        "EFT1RenormalisableRGEStatus",
        "EFT1WilsonRGEStatus",
        "EFT1FlavorSeedStatus",
        "EFT1FullFlavorBridgeStatus",
    ):
        assert key in backend

    for filename in (
        "eft1_rgbeta_rge.json",
        "eft1_wilson_rge.json",
        "eft1_wilson_flavor_seed.json",
        "eft1_wilson_flavor_at_S_threshold.json",
    ):
        assert filename in backend


def test_terminal_messages_use_physical_intermediate_eft_language() -> None:
    backend = _source(BACKEND_PATH)

    assert "starting RGBeta EFT1" not in backend
    assert "starting EFT1 dimension-five" not in backend
    assert "full-flavor EFT1 bridge failed" not in backend
    assert "scalar-only intermediate-EFT" in backend
