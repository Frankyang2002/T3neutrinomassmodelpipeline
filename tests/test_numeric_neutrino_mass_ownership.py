"""Architecture contracts for numerical C5 -> m_nu ownership."""

from __future__ import annotations

from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
CANONICAL_RUNNING = PROJECT_ROOT / "Numerical" / "SMWeinbergEvolution.py"
HISTORICAL_RUNNING = PROJECT_ROOT / "Numerical" / "WeinbergRunning.py"
CANONICAL_STAGE = PROJECT_ROOT / "Numerical" / "SMWeinbergStage.py"
HISTORICAL_STAGE = PROJECT_ROOT / "Numerical" / "WeinbergStage.py"
MASS = PROJECT_ROOT / "physics" / "NeutrinoMass.py"


def _source(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def test_numeric_mass_conversion_is_implemented_under_physics() -> None:
    source = _source(MASS)

    assert "def numerical_neutrino_mass_matrix(" in source
    assert "return -(vev**2 / 2.0) * c5" in source


def test_canonical_numerical_evolution_reexports_mass_helper_from_physics() -> None:
    source = _source(CANONICAL_RUNNING)

    assert (
        "numerical_neutrino_mass_matrix as neutrino_mass_matrix"
        in source
    )
    assert "def neutrino_mass_matrix(" not in source


def test_historical_numerical_running_is_compatibility_only() -> None:
    source = _source(HISTORICAL_RUNNING)

    assert "from Numerical.SMWeinbergEvolution import" in source
    assert "def neutrino_mass_matrix(" not in source


def test_canonical_numerical_stage_calls_physics_mass_conversion_directly() -> None:
    source = _source(CANONICAL_STAGE)

    assert "from physics.NeutrinoMass import numerical_neutrino_mass_matrix" in source
    assert "mass_low = numerical_neutrino_mass_matrix(" in source
    assert "Numerical.WeinbergRunning" not in source


def test_historical_numerical_stage_is_compatibility_only() -> None:
    source = _source(HISTORICAL_STAGE)

    assert "from Numerical.SMWeinbergStage import" in source
    assert "def run_numerical_weinberg_stage(" not in source
