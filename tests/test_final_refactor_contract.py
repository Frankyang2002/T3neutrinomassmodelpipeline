"""Final architecture/documentation contracts for the structural refactor."""

from __future__ import annotations

from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def _text(relative: str) -> str:
    return (PROJECT_ROOT / relative).read_text(encoding="utf-8")


def test_central_documentation_matches_the_refactored_pipeline_boundaries() -> None:
    pipeline_doc = _text("INFO_PIPELINE.md")

    required = (
        "pipeline.py",
        "common/PipelinePlan.py",
        "common/T3Fields.py",
        "model/T3Study.py",
        "Lagrangian/T3ModelMatching.py",
        "RGE/running/IntermediateEFTRunning.py",
        "RGE/running/backends/ScalarOnlyAfterFermion.py",
        "RGE/running/intermediate/ScalarOnlyWilsonTensorRGE.py",
        "Numerical/IntermediateScalarState.py",
        "validation/",
        "physics/LowEnergyNeutrino.py",
        "Reports/PipelineReports.py",
    )
    for name in required:
        assert name in pipeline_doc


def test_documentation_states_the_project_hypercharge_conversion_explicitly() -> None:
    for relative in ("INFO_PIPELINE.md", "INFO_LAGRANGIAN.md"):
        source = _text(relative)
        assert "Q=T_3+Y" in source.replace(" ", "")
        assert "Y_{\\rmRZY}=2Y" in source.replace(" ", "")
        assert "Y(S_1)=\\frac{\\alpha}{2}" in source.replace(" ", "")
        assert "Y(F)=\\frac{\\alpha+1}{2}" in source.replace(" ", "")


def test_scalar_first_dimension_six_limitation_is_documented() -> None:
    for relative in ("INFO_PIPELINE.md", "INFO_RGE.md", "INFO_LAGRANGIAN.md"):
        source = _text(relative)
        assert "--allow-truncated-scalar-first" in source
        assert "dimension-six" in source.lower()
        assert "non-authoritative" in source.lower() or "not a complete" in source.lower()


def test_historical_eft1_implementation_package_is_retired() -> None:
    assert not (
        PROJECT_ROOT
        / "RGE"
        / "running"
        / "intermediate"
        / "LegacyEFT1Compatibility.py"
    ).exists()

    source_roots = (
        PROJECT_ROOT / "pipeline.py",
        PROJECT_ROOT / "RGE" / "running" / "IntermediateEFTRunning.py",
        PROJECT_ROOT / "RGE" / "running" / "backends",
        PROJECT_ROOT / "RGE" / "running" / "intermediate",
        PROJECT_ROOT / "validation",
        PROJECT_ROOT / "physics",
        PROJECT_ROOT / "Numerical",
    )

    offenders: list[str] = []
    for root in source_roots:
        paths = [root] if root.is_file() else root.rglob("*.py")
        for path in paths:
            source = path.read_text(encoding="utf-8")
            if "from RGE.running.eft1" in source or "from Numerical.EFT1" in source:
                offenders.append(path.relative_to(PROJECT_ROOT).as_posix())

    assert offenders == []


def test_shared_scalar_identity_is_documented_as_physical_s_not_parallel_architecture() -> None:
    fields = _text("common/T3Fields.py")
    pipeline_doc = _text("INFO_PIPELINE.md")

    assert 'SHARED_SCALAR_FIELDS: tuple[str, ...] = ("S",)' in fields
    assert "F, S" in pipeline_doc
    assert "formal topology roles `S1` and `S2`" in pipeline_doc
    assert "only at the matching boundary" in pipeline_doc


def test_docs_require_external_smoke_regression_after_python_tests() -> None:
    source = _text("INFO_PIPELINE.md")
    pytest_position = source.index("python -m pytest -q")
    smoke_position = source.index("python pipeline.py --smoke")
    assert pytest_position < smoke_position
    assert "Wolfram/Matchete" in source
