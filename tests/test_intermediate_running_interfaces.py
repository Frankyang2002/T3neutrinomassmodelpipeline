"""Architecture checks for the physical intermediate-EFT running implementation."""

from __future__ import annotations

import ast
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
INTERMEDIATE_DIR = PROJECT_ROOT / "RGE" / "running" / "intermediate"
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


def test_scalar_only_public_interfaces_use_physical_names() -> None:
    expected = {
        "ScalarOnlyRenormalisableRunning.py": "run_scalar_only_renormalisable_rge",
        "ScalarOnlyWilsonRunning.py": "run_scalar_only_dimension_five_wilson_rge",
        "FlavorWilsonBoundary.py": "export_full_flavor_wilson_boundary",
        "DirectWeinbergRunning.py": "build_direct_weinberg_running",
        "ScalarThresholdMatching.py": "resume_scalar_threshold_with_running",
    }

    for filename, public_name in expected.items():
        path = INTERMEDIATE_DIR / filename
        assert path.is_file()
        _function(path, public_name)
        assert "LegacyEFT1Compatibility" not in _source(path)
        assert "from RGE.running.eft1" not in _source(path)


def test_direct_weinberg_module_owns_insertion_export() -> None:
    path = INTERMEDIATE_DIR / "DirectWeinbergRunning.py"
    _function(path, "export_direct_weinberg_insertion")
    source = _source(path)
    assert "mathematica_code" in source
    assert "LegacyEFT1Compatibility" not in source


def test_scalar_threshold_module_owns_resume_implementation() -> None:
    path = INTERMEDIATE_DIR / "ScalarThresholdMatching.py"
    source = _source(path)
    _function(path, "resume_scalar_threshold_with_running")
    assert "subprocess.Popen" in source
    assert "ScalarThresholdContinuation" in source


def test_dimension_five_rge_module_owns_tensor_calculation() -> None:
    path = INTERMEDIATE_DIR / "ScalarOnlyWilsonTensorRGE.py"
    source = _source(path)
    _function(path, "run_scalar_only_wilson_rge")
    assert "calculate_complete_wilson_tensor_rge" in source
    assert "ScalarOnlyRGEContext" in source


def test_production_backend_uses_physical_intermediate_interfaces() -> None:
    backend = _source(BACKEND_PATH)

    assert "from RGE.running.eft1" not in backend
    assert "LegacyEFT1Compatibility" not in backend

    for public_name in (
        "run_scalar_only_renormalisable_rge",
        "run_scalar_only_dimension_five_wilson_rge",
        "export_full_flavor_wilson_boundary",
        "build_direct_weinberg_running",
        "export_direct_weinberg_insertion",
        "resume_scalar_threshold_with_running",
    ):
        assert public_name in backend


def test_historical_output_contract_remains_serialized_compatibility_only() -> None:
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
