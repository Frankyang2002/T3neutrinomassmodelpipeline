"""Architecture contracts for final SM+Weinberg numerical naming."""

from __future__ import annotations

from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RUNNING = PROJECT_ROOT / "Numerical" / "running"
EVOLUTION = RUNNING / "SMWeinbergEvolution.py"
STAGE = RUNNING / "SMWeinbergStage.py"
TRAJECTORY = RUNNING / "WeinbergTrajectory.py"
LOW_ENERGY = PROJECT_ROOT / "physics" / "LowEnergyNeutrino.py"


def _source(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def test_canonical_evolution_names_describe_final_sm_weinberg_eft() -> None:
    source = _source(EVOLUTION)

    assert "class SMWeinbergInitialConditions:" in source
    assert "class SMWeinbergEvolutionResult:" in source
    assert "def evolve_sm_weinberg(" in source
    assert "solve_ivp" in source


def test_historical_numerical_weinberg_modules_are_removed() -> None:
    assert not (PROJECT_ROOT / "Numerical" / "WeinbergRunning.py").exists()
    assert not (PROJECT_ROOT / "Numerical" / "WeinbergStage.py").exists()


def test_canonical_stage_uses_canonical_evolution_names() -> None:
    source = _source(STAGE)

    assert "from Numerical.running.SMWeinbergEvolution import" in source
    assert "def run_sm_weinberg_numerical_stage(" in source
    assert "SMWeinbergInitialConditions(" in source
    assert "evolve_sm_weinberg(" in source
    assert "Numerical.WeinbergRunning" not in source


def test_trajectory_and_low_energy_facade_use_canonical_modules() -> None:
    trajectory = _source(TRAJECTORY)
    low_energy = _source(LOW_ENERGY)

    assert "from Numerical.running.SMWeinbergEvolution import" in trajectory
    assert "Numerical.WeinbergRunning" not in trajectory

    assert "from Numerical.running.SMWeinbergStage import" in low_energy
    assert "_run_sm_weinberg_numerical_stage(" in low_energy
    assert "Numerical.WeinbergStage" not in low_energy
