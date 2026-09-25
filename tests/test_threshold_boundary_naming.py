"""Architecture contracts for numerical threshold-boundary naming."""

from __future__ import annotations

from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
FERMION = PROJECT_ROOT / "Numerical" / "FermionThresholdBoundary.py"
FERMION_COMPAT = PROJECT_ROOT / "Numerical" / "ThresholdMatching.py"
SCALAR = PROJECT_ROOT / "Numerical" / "ScalarThresholdBoundary.py"
SCALAR_COMPAT = PROJECT_ROOT / "Numerical" / "FinalSMBoundary.py"
TRAJECTORY = PROJECT_ROOT / "Numerical" / "T3Trajectory.py"


def _source(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def test_fermion_threshold_module_describes_boundary_not_matching() -> None:
    source = _source(FERMION)

    assert "def fermion_threshold_masses(" in source
    assert "def build_intermediate_scalar_boundary(" in source
    assert "Matchete/RGE matching pipeline" in source


def test_scalar_threshold_module_describes_final_sm_boundary() -> None:
    source = _source(SCALAR)

    assert "class FinalSMBoundaryState:" in source
    assert "def build_final_sm_boundary(" in source
    assert "def build_sm_weinberg_initial_conditions(" in source


def test_historical_threshold_modules_are_compatibility_only() -> None:
    fermion = _source(FERMION_COMPAT)
    scalar = _source(SCALAR_COMPAT)

    assert "from Numerical.FermionThresholdBoundary import" in fermion
    assert "def project_after_fermion_threshold(" not in fermion

    assert "from Numerical.ScalarThresholdBoundary import" in scalar
    assert "def project_after_scalar_threshold(" not in scalar
    assert "def build_weinberg_initial_conditions(" not in scalar


def test_t3_trajectory_uses_physical_boundary_names() -> None:
    source = _source(TRAJECTORY)

    assert "from Numerical.FermionThresholdBoundary import" in source
    assert "from Numerical.ScalarThresholdBoundary import" in source
    assert "build_intermediate_scalar_boundary(" in source
    assert "build_final_sm_boundary(" in source
    assert "build_sm_weinberg_initial_conditions(" in source

    assert "Numerical.ThresholdMatching" not in source
    assert "Numerical.FinalSMBoundary" not in source
