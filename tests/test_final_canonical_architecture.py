"""Final canonical-caller and documentation architecture checks."""

from __future__ import annotations

from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def _text(relative: str) -> str:
    return (PROJECT_ROOT / relative).read_text(encoding="utf-8")


def test_internal_numerical_callers_use_physics_ownership_directly() -> None:
    fit = _text("Numerical/OscillationFit.py")
    scan = _text("Numerical/ParameterScan.py")

    assert "from physics.NeutrinoObservables import" in fit
    assert "RGE.phenomenology.NeutrinoObservables" not in fit

    assert "from physics.NeutrinoTrajectory import" in scan
    assert "from Numerical.WeinbergTrajectory import" not in scan


def test_scan_cli_uses_canonical_final_c5_trajectory_adapter() -> None:
    source = _text("Numerical/ScanCLI.py")

    assert "from Numerical.FinalC5TrajectoryAdapter import FinalC5TrajectoryEvaluator" in source
    assert "FinalC5TrajectoryEvaluator(" in source
    assert "from Numerical.FinalC5Bridge import" not in source


def test_t3_trajectory_exposes_descriptive_intermediate_aliases() -> None:
    source = _text("Numerical/T3Trajectory.py")

    assert "def intermediate_initial_state(" in source
    assert "def intermediate(" in source
    assert "def intermediate_threshold_state(" in source
    assert "def eft1_threshold_state(" in source


def test_docs_use_current_canonical_ownership() -> None:
    pipeline = _text("INFO_PIPELINE.md")
    rge = _text("INFO_RGE.md")
    lagrangian = _text("INFO_LAGRANGIAN.md")

    for token in (
        "RGE/matching/FinalWeinbergCoefficient.py",
        "RGE/running/weinberg/WeinbergRGE.py",
        "RGE/running/weinberg/FullFlavorWeinbergStage.py",
        "Numerical/SMWeinbergEvolution.py",
        "Numerical/FinalC5TrajectoryAdapter.py",
        "physics/NeutrinoObservables.py",
    ):
        assert token in pipeline
        assert token in rge

    for token in (
        "Lagrangian/ModelValidation.py",
        "Lagrangian/WolframRunner.py",
        "Lagrangian/MatchingResults.py",
    ):
        assert token in pipeline
        assert token in lagrangian
