"""Canonical ownership assertions based on live module paths and consolidated docs."""
from __future__ import annotations

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def _text(relative: str) -> str:
    return (PROJECT_ROOT / relative).read_text(encoding="utf-8-sig")


def test_internal_numerical_callers_use_physics_ownership_directly() -> None:
    fit = _text("Numerical/fitting/OscillationFit.py")
    scan = _text("Numerical/fitting/ParameterScan.py")
    assert "from physics.NeutrinoObservables import" in fit
    assert "RGE.phenomenology.NeutrinoObservables" not in fit
    assert "from physics.NeutrinoTrajectory import" in scan
    assert "from Numerical.WeinbergTrajectory import" not in scan
    assert "from Numerical.running.T3Trajectory import" in scan


def test_scan_cli_uses_canonical_final_c5_trajectory_adapter() -> None:
    source = _text("Numerical/fitting/ScanCLI.py")
    assert ("from Numerical.running.FinalC5TrajectoryAdapter "
            "import FinalC5TrajectoryEvaluator") in source
    assert "FinalC5TrajectoryEvaluator(" in source
    assert "from Numerical.FinalC5Bridge import" not in source


def test_t3_trajectory_exposes_descriptive_intermediate_aliases() -> None:
    source = _text("Numerical/running/T3Trajectory.py")
    for name in ("intermediate_initial_state", "intermediate", "intermediate_threshold_state", "eft1_threshold_state"):
        assert f"def {name}(" in source


def test_docs_use_current_canonical_ownership() -> None:
    overview = _text("README.md")
    physics = _text("docs/PHYSICS.md")
    numerical = _text("docs/NUMERICAL.md")
    assert "pipeline.py" in overview
    assert "Weinberg" in overview and "Weinberg" in physics
    assert "Matchete" in physics
    assert "Numerical/" in numerical
    for relative in (
        "RGE/matching/FinalWeinbergCoefficient.py",
        "RGE/running/weinberg/WeinbergRGE.py",
        "RGE/running/weinberg/FullFlavorWeinbergStage.py",
        "Numerical/running/SMWeinbergEvolution.py",
        "Numerical/running/FinalC5TrajectoryAdapter.py",
        "physics/NeutrinoObservables.py",
        "Lagrangian/ModelValidation.py",
        "Lagrangian/WolframRunner.py",
        "Lagrangian/MatchingResults.py",
    ):
        assert (PROJECT_ROOT / relative).is_file(), relative
