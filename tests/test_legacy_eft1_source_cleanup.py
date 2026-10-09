"""Regression contracts for retiring historical ordinal EFT1 source imports."""
from __future__ import annotations

import ast
import importlib
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
INTERMEDIATE_DIR = PROJECT_ROOT / "RGE" / "running" / "intermediate"
OLD_DIR = PROJECT_ROOT / "RGE" / "running" / "eft1"
SCAN_ROOTS = tuple(PROJECT_ROOT / p for p in (
    "pipeline.py", "common", "Lagrangian", "model", "Numerical", "physics",
    "Reports", "RGE", "studies", "validation", "tests",
))
LEGACY_PREFIXES = ("RGE.running.eft1", "Numerical.EFT1")


def _python_files(root: Path):
    if root.is_file():
        if root.suffix == ".py":
            yield root
        return
    if root.exists():
        yield from (path for path in root.rglob("*.py") if "__pycache__" not in path.parts)


def _legacy_imports() -> list[str]:
    offenders: list[str] = []
    for root in SCAN_ROOTS:
        for path in _python_files(root):
            tree = ast.parse(path.read_text(encoding="utf-8-sig"), filename=str(path))
            for node in ast.walk(tree):
                if isinstance(node, ast.ImportFrom):
                    module = node.module or ""
                    if module.startswith(LEGACY_PREFIXES):
                        offenders.append(f"{path.relative_to(PROJECT_ROOT).as_posix()}:{node.lineno}:{module}")
                elif isinstance(node, ast.Import):
                    for alias in node.names:
                        if alias.name.startswith(LEGACY_PREFIXES):
                            offenders.append(f"{path.relative_to(PROJECT_ROOT).as_posix()}:{node.lineno}:{alias.name}")
    return offenders


def test_historical_eft1_source_package_is_gone() -> None:
    assert not OLD_DIR.exists()
    assert not (INTERMEDIATE_DIR / "LegacyEFT1Compatibility.py").exists()


def test_physical_intermediate_modules_own_the_implementations() -> None:
    for name in (
        "IntermediateMatcheteParsing.py", "ScalarOnlyTensorAdapters.py",
        "ScalarOnlyWilsonTensorRGE.py", "ScalarOnlyWilsonFlow.py",
        "DirectWeinbergRunning.py", "ScalarThresholdMatching.py",
    ):
        assert (INTERMEDIATE_DIR / name).is_file(), name


def test_repository_has_no_python_imports_of_retired_eft1_modules() -> None:
    assert _legacy_imports() == []


def test_group_factor_stack_imports_after_cleanup() -> None:
    for module in (
        "RGE.group_factors.core.MixingQuarticTensorAlgebra",
        "RGE.group_factors.core.ScalarQuarticBasis",
        "RGE.group_factors.core.YukawaGroupFactors",
        "RGE.group_factors.recoupling.MixingQuarticRecoupling",
        "RGE.group_factors.recoupling.RGBetaMixingQuarticRecoupling",
        "Reports.GroupFactorReports",
    ):
        importlib.import_module(module)


def test_pipeline_imports_after_cleanup() -> None:
    importlib.import_module("pipeline")


def test_thesis_figure_generator_is_headless() -> None:
    source = (PROJECT_ROOT / "Numerical/plotting/ThesisResultFigures.py").read_text(encoding="utf-8-sig")
    assert source.index('matplotlib.use("Agg")') < source.index("import matplotlib.pyplot as plt")
