"""Contracts for retiring the historical ordinal EFT1 source architecture."""

from __future__ import annotations

import ast
import importlib
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
INTERMEDIATE_DIR = PROJECT_ROOT / "RGE" / "running" / "intermediate"
OLD_DIR = PROJECT_ROOT / "RGE" / "running" / "eft1"

SCAN_ROOTS = (
    PROJECT_ROOT / "pipeline.py",
    PROJECT_ROOT / "common",
    PROJECT_ROOT / "Lagrangian",
    PROJECT_ROOT / "model",
    PROJECT_ROOT / "Numerical",
    PROJECT_ROOT / "physics",
    PROJECT_ROOT / "Reports",
    PROJECT_ROOT / "RGE",
    PROJECT_ROOT / "studies",
    PROJECT_ROOT / "validation",
    PROJECT_ROOT / "tests",
)

LEGACY_PREFIXES = (
    "RGE.running.eft1",
    "Numerical.EFT1",
)


def _python_files(root: Path):
    if root.is_file():
        if root.suffix == ".py":
            yield root
        return
    if not root.exists():
        return
    yield from (
        path
        for path in root.rglob("*.py")
        if "__pycache__" not in path.parts
    )


def _legacy_imports() -> list[str]:
    offenders: list[str] = []
    for root in SCAN_ROOTS:
        for path in _python_files(root):
            tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
            for node in ast.walk(tree):
                if isinstance(node, ast.ImportFrom):
                    module = node.module or ""
                    if module.startswith(LEGACY_PREFIXES):
                        offenders.append(
                            f"{path.relative_to(PROJECT_ROOT).as_posix()}:{node.lineno}:{module}"
                        )
                elif isinstance(node, ast.Import):
                    for alias in node.names:
                        if alias.name.startswith(LEGACY_PREFIXES):
                            offenders.append(
                                f"{path.relative_to(PROJECT_ROOT).as_posix()}:{node.lineno}:{alias.name}"
                            )
    return offenders


def test_historical_eft1_source_package_is_gone() -> None:
    assert not OLD_DIR.exists()
    assert not (
        INTERMEDIATE_DIR / "LegacyEFT1Compatibility.py"
    ).exists()


def test_physical_intermediate_modules_own_the_implementations() -> None:
    required = (
        "IntermediateMatcheteParsing.py",
        "ScalarOnlyTensorAdapters.py",
        "ScalarOnlyWilsonTensorRGE.py",
        "ScalarOnlyWilsonFlow.py",
        "DirectWeinbergRunning.py",
        "ScalarThresholdMatching.py",
    )
    for name in required:
        assert (INTERMEDIATE_DIR / name).is_file(), name


def test_repository_has_no_python_imports_of_retired_eft1_modules() -> None:
    assert _legacy_imports() == []


def test_group_factor_stack_imports_after_cleanup() -> None:
    # These imports caught stale dependencies that the earlier cleanup scan
    # missed because the smoke pipeline imports the report stack lazily through
    # Reports.PipelineReports.
    modules = (
        "RGE.group_factors.core.MixingQuarticTensorAlgebra",
        "RGE.group_factors.core.ScalarQuarticBasis",
        "RGE.group_factors.core.YukawaGroupFactors",
        "RGE.group_factors.recoupling.MixingQuarticRecoupling",
        "RGE.group_factors.recoupling.RGBetaMixingQuarticRecoupling",
        "Reports.GroupFactorReports",
    )
    for module in modules:
        importlib.import_module(module)


def test_pipeline_imports_after_cleanup() -> None:
    # Importing the pipeline must be possible without invoking Wolfram/Matchete.
    # This catches stale import paths before an expensive smoke run starts.
    importlib.import_module("pipeline")


def test_thesis_figure_generator_is_headless() -> None:
    path = PROJECT_ROOT / "Numerical" / "ThesisResultFigures.py"
    source = path.read_text(encoding="utf-8")
    use_index = source.index('matplotlib.use("Agg")')
    pyplot_index = source.index("import matplotlib.pyplot as plt")
    assert use_index < pyplot_index
