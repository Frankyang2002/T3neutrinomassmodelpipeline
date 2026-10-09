"""Architecture and documentation contracts for the current physical-stage refactor."""
from __future__ import annotations

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def _text(relative: str) -> str:
    return (PROJECT_ROOT / relative).read_text(encoding="utf-8-sig")


def test_central_documentation_matches_the_refactored_pipeline_boundaries() -> None:
    pipeline_doc = _text("INFO_PIPELINE.md")
    assert "pipeline.py" in pipeline_doc
    assert "scalar-only intermediate EFT" in pipeline_doc
    assert "Weinberg" in pipeline_doc
    assert "numerical" in pipeline_doc.lower()
    # Ownership is checked against code, rather than requiring an exhaustive
    # path listing in a human-facing overview document.
    for relative in (
        "common/PipelinePlan.py", "common/T3Fields.py", "model/T3Study.py",
        "Lagrangian/T3ModelMatching.py", "RGE/running/IntermediateEFTRunning.py",
        "RGE/running/backends/ScalarOnlyAfterFermion.py",
        "RGE/running/intermediate/ScalarOnlyWilsonTensorRGE.py",
        "Numerical/running/IntermediateScalarState.py", "physics/LowEnergyNeutrino.py",
        "Reports/PipelineReports.py",
    ):
        assert (PROJECT_ROOT / relative).is_file(), relative


def test_numerical_running_root_forwarders_are_retired() -> None:
    retired = (
        "UVRunner.py", "T3Trajectory.py", "IntermediateScalarState.py",
        "IntermediateScalarStateVector.py", "IntermediateScalarRGBetaEvaluator.py",
        "IntermediateScalarRunner.py", "FermionThresholdBoundary.py",
        "ScalarThresholdBoundary.py", "SMWeinbergEvolution.py", "SMWeinbergStage.py",
        "WeinbergTrajectory.py", "FinalC5TrajectoryAdapter.py",
    )
    for name in retired:
        assert not (PROJECT_ROOT / "Numerical" / name).exists()
        assert (PROJECT_ROOT / "Numerical" / "running" / name).is_file()


def test_numerical_fitting_root_forwarders_are_retired() -> None:
    for name in (
        "OscillationFit.py", "ParameterScan.py", "OscillationOptimizer.py",
        "SobolBenchmarkSearch.py", "BenchmarkSensitivity.py", "ScanCLI.py",
    ):
        assert not (PROJECT_ROOT / "Numerical" / name).exists()
        assert (PROJECT_ROOT / "Numerical" / "fitting" / name).is_file()


def test_numerical_diagnostics_root_forwarders_are_retired() -> None:
    for name in (
        "BestFitDiagnostics.py", "RunningDiagnostics.py", "IntermediateWeinbergDiagnostics.py",
        "FinalC5ContributionDiagnostics.py", "NumericalConfigCheck.py",
    ):
        assert not (PROJECT_ROOT / "Numerical" / name).exists()
        assert (PROJECT_ROOT / "Numerical" / "diagnostics" / name).is_file()


def test_numerical_plotting_root_forwarders_are_retired() -> None:
    for name in (
        "RunningResultFigures.py", "PlotOptimizerResults.py", "PlotScanResults.py",
        "ThesisResultFigures.py", "InteractiveModelComparison.py", "BuildDisplayResults.py",
    ):
        assert not (PROJECT_ROOT / "Numerical" / name).exists()
        assert (PROJECT_ROOT / "Numerical" / "plotting" / name).is_file()


def test_numerical_orchestration_root_forwarders_are_retired() -> None:
    for name in ("NumericalConfig.py", "PipelineNumericalResults.py"):
        assert not (PROJECT_ROOT / "Numerical" / name).exists()
        assert (PROJECT_ROOT / "Numerical" / "orchestration" / name).is_file()


def test_historical_final_c5_bridge_is_retired() -> None:
    assert not (PROJECT_ROOT / "Numerical" / "FinalC5Bridge.py").exists()
    assert (PROJECT_ROOT / "Numerical/running/FinalC5TrajectoryAdapter.py").is_file()


def test_documentation_states_the_project_hypercharge_conversion_explicitly() -> None:
    for relative in ("INFO_PIPELINE.md", "INFO_LAGRANGIAN.md"):
        source = _text(relative).replace(" ", "")
        assert "Q=T_3+Y" in source
        assert "Y_{\\rmRZY}=2Y" in source
        assert "Y(S_1)=\\frac{\\alpha}{2}" in source
        assert "Y(F)=\\frac{\\alpha+1}{2}" in source


def test_scalar_first_is_documented_as_outside_production_scope() -> None:
    pipeline = _text("INFO_PIPELINE.md").lower()
    assert "scalar-first" in pipeline
    assert "dimension-six" in pipeline
    assert "production" in pipeline
    for relative in ("INFO_PIPELINE.md", "INFO_RGE.md", "INFO_LAGRANGIAN.md"):
        source = _text(relative)
        assert "--allow-truncated-scalar-first" not in source
        assert "--eft-max-dimension" not in source


def test_historical_eft1_implementation_package_is_retired() -> None:
    assert not (PROJECT_ROOT / "RGE/running/intermediate/LegacyEFT1Compatibility.py").exists()
    source_roots = (
        PROJECT_ROOT / "pipeline.py", PROJECT_ROOT / "RGE/running/IntermediateEFTRunning.py",
        PROJECT_ROOT / "RGE/running/backends", PROJECT_ROOT / "RGE/running/intermediate",
        PROJECT_ROOT / "validation", PROJECT_ROOT / "physics", PROJECT_ROOT / "Numerical",
    )
    offenders = []
    for root in source_roots:
        paths = [root] if root.is_file() else root.rglob("*.py")
        for path in paths:
            source = path.read_text(encoding="utf-8-sig")
            if "from RGE.running.eft1" in source or "from Numerical.EFT1" in source:
                offenders.append(path.relative_to(PROJECT_ROOT).as_posix())
    assert offenders == []


def test_shared_scalar_identity_is_documented_as_physical_s_not_parallel_architecture() -> None:
    fields = _text("common/T3Fields.py")
    pipeline_doc = _text("INFO_PIPELINE.md")
    assert 'SHARED_SCALAR_FIELDS: tuple[str, ...] = ("S",)' in fields
    assert "F, S" in pipeline_doc
    assert "shared-scalar" in pipeline_doc.lower()
    assert "one physical scalar" in pipeline_doc.lower()


def test_project_helper_scripts_are_not_root_level() -> None:
    assert not (PROJECT_ROOT / "comparison_test.py").exists()
    assert not (PROJECT_ROOT / "regression.py").exists()
    assert (PROJECT_ROOT / "scripts/run_comparison_study.py").is_file()
    assert (PROJECT_ROOT / "scripts/run_regression.py").is_file()


def test_docs_require_external_smoke_regression_after_python_tests() -> None:
    """The actual smoke CLI must exist; docs may link to INFO_TESTS instead."""
    from common.PipelineCLI import build_argument_parser
    parser = build_argument_parser()
    option_strings = {o for a in parser._actions for o in a.option_strings}
    assert "--smoke" in option_strings
    info = _text("INFO_TESTS.md")
    assert "pytest" in info.lower()
    assert "smoke" in info.lower()
