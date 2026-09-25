"""Architecture contracts for Weinberg trajectory ownership."""

from __future__ import annotations

from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
NUMERICAL = PROJECT_ROOT / "Numerical" / "WeinbergTrajectory.py"
PHYSICS = PROJECT_ROOT / "physics" / "NeutrinoTrajectory.py"


def _source(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def test_numerical_trajectory_owns_integration_only() -> None:
    source = _source(NUMERICAL)

    assert "solve_ivp" in source
    assert "def run_weinberg_trajectory(" in source
    assert "class WeinbergTrajectoryResult" in source

    assert "def neutrino_mass_from_c5(" not in source
    assert "def charged_lepton_mass_basis_matrix(" not in source
    assert "def scale_dependent_neutrino_observables(" not in source
    assert "calculate_neutrino_observables" not in source


def test_physics_trajectory_owns_physical_interpretation() -> None:
    source = _source(PHYSICS)

    assert "def neutrino_mass_from_c5(" in source
    assert "def charged_lepton_mass_basis_matrix(" in source
    assert "def scale_dependent_neutrino_observables(" in source
    assert "calculate_neutrino_observables" in source

    assert "solve_ivp" not in source
    assert "Numerical.WeinbergTrajectory" not in source


def test_numerical_module_retains_legacy_import_surface() -> None:
    source = _source(NUMERICAL)

    assert "from physics.NeutrinoTrajectory import" in source
    assert "ScaleDependentNeutrinoPoint" in source
    assert "charged_lepton_mass_basis_matrix" in source
    assert "neutrino_mass_from_c5" in source
    assert "scale_dependent_neutrino_observables" in source
