"""Architecture contracts for final SM+Weinberg numerical naming."""

from __future__ import annotations

from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
EVOLUTION = PROJECT_ROOT / "Numerical" / "SMWeinbergEvolution.py"
EVOLUTION_COMPAT = PROJECT_ROOT / "Numerical" / "WeinbergRunning.py"
STAGE = PROJECT_ROOT / "Numerical" / "SMWeinbergStage.py"
STAGE_COMPAT = PROJECT_ROOT / "Numerical" / "WeinbergStage.py"
TRAJECTORY = PROJECT_ROOT / "Numerical" / "WeinbergTrajectory.py"
LOW_ENERGY = PROJECT_ROOT / "physics" / "LowEnergyNeutrino.py"


def _source(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def test_canonical_evolution_names_describe_final_sm_weinberg_eft() -> None:
    source = _source(EVOLUTION)

    assert "class SMWeinbergInitialConditions:" in source
    assert "class SMWeinbergEvolutionResult:" in source
    assert "def evolve_sm_weinberg(" in source
    assert "solve_ivp" in source


def test_historical_runner_is_compatibility_only() -> None:
    source = _source(EVOLUTION_COMPAT)

    assert "from Numerical.SMWeinbergEvolution import" in source
    assert "def evolve_weinberg(" not in source
    assert "solve_ivp" not in source


def test_canonical_stage_uses_canonical_evolution_names() -> None:
    source = _source(STAGE)

    assert "from Numerical.SMWeinbergEvolution import" in source
    assert "def run_sm_weinberg_numerical_stage(" in source
    assert "SMWeinbergInitialConditions(" in source
    assert "evolve_sm_weinberg(" in source
    assert "Numerical.WeinbergRunning" not in source


def test_historical_stage_is_compatibility_only() -> None:
    source = _source(STAGE_COMPAT)

    assert "from Numerical.SMWeinbergStage import" in source
    assert "run_numerical_weinberg_stage = run_sm_weinberg_numerical_stage" in source
    assert "def run_numerical_weinberg_stage(" not in source


def test_trajectory_and_low_energy_facade_use_canonical_modules() -> None:
    trajectory = _source(TRAJECTORY)
    low_energy = _source(LOW_ENERGY)

    assert "from Numerical.SMWeinbergEvolution import" in trajectory
    assert "Numerical.WeinbergRunning" not in trajectory

    assert "from Numerical.SMWeinbergStage import" in low_energy
    assert "_run_sm_weinberg_numerical_stage(" in low_energy
    assert "Numerical.WeinbergStage" not in low_energy
