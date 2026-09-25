"""Keep historical ordinal-EFT imports behind one explicit compatibility boundary."""

from __future__ import annotations

import ast
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
INTERMEDIATE_DIR = PROJECT_ROOT / "RGE" / "running" / "intermediate"
COMPATIBILITY_PATH = INTERMEDIATE_DIR / "LegacyEFT1Compatibility.py"


def _source(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def test_only_compatibility_module_imports_historical_eft1_package() -> None:
    offenders: list[str] = []
    for path in sorted(INTERMEDIATE_DIR.glob("*.py")):
        if path == COMPATIBILITY_PATH:
            continue
        if "RGE.running.eft1" in _source(path):
            offenders.append(path.name)

    assert offenders == []


def test_compatibility_boundary_contains_all_required_legacy_entry_points() -> None:
    source = _source(COMPATIBILITY_PATH)

    required_imports = (
        "run_rgbeta_t3_eft1",
        "run_eft1_wilson_rge",
        "run_eft1_wilson_transport",
        "run_flavor_seed_export",
        "build_direct_weinberg_flavor_transport",
        "export_direct_weinberg_matchete",
        "Threshold2Continuation",
        "rerun_threshold2_with_running",
    )
    for name in required_imports:
        assert name in source

    required_aliases = (
        "legacy_scalar_only_renormalisable_rge",
        "legacy_scalar_only_dimension_five_wilson_rge",
        "legacy_component_wilson_transport",
        "legacy_export_full_flavor_wilson_boundary",
        "legacy_build_direct_weinberg_running",
        "legacy_export_direct_weinberg_insertion",
        "LegacyScalarThresholdContinuation",
        "legacy_resume_scalar_threshold_with_running",
    )
    for name in required_aliases:
        assert name in source


def test_compatibility_boundary_has_no_physics_logic() -> None:
    tree = ast.parse(_source(COMPATIBILITY_PATH))

    functions = [node for node in tree.body if isinstance(node, ast.FunctionDef)]
    classes = [node for node in tree.body if isinstance(node, ast.ClassDef)]

    assert functions == []
    assert classes == []
    assert len(_source(COMPATIBILITY_PATH).splitlines()) < 65


def test_legacy_boundary_is_explicitly_exported_and_documented() -> None:
    source = _source(COMPATIBILITY_PATH)
    tree = ast.parse(source)

    assert ast.get_docstring(tree)
    assert "serialized output contract" in ast.get_docstring(tree)
    assert "__all__" in source
