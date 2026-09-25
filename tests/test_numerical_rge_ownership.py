"""Architecture contracts for RGE-model versus numerical-runner ownership."""

from __future__ import annotations

from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RGE_WEINBERG = PROJECT_ROOT / "RGE" / "running" / "weinberg"
NUMERICAL = PROJECT_ROOT / "Numerical"


def test_final_eft_rge_model_stays_under_rge() -> None:
    path = RGE_WEINBERG / "WeinbergRGE.py"
    source = path.read_text(encoding="utf-8")

    assert path.is_file()
    assert "beta_weinberg_matrix" in source
    assert "sm_weinberg_beta" in source
    assert "solve_ivp" not in source


def test_numerical_weinberg_integration_lives_under_numerical() -> None:
    runner = NUMERICAL / "WeinbergRunning.py"
    stage = NUMERICAL / "WeinbergStage.py"

    assert runner.is_file()
    assert stage.is_file()
    assert "solve_ivp" in runner.read_text(encoding="utf-8")
    assert "run_numerical_weinberg_stage" in stage.read_text(encoding="utf-8")


def test_old_mixed_ownership_modules_are_removed() -> None:
    assert not (RGE_WEINBERG / "WeinbergRunning.py").exists()
    assert not (RGE_WEINBERG / "NumericalWeinbergStage.py").exists()


def test_symbolic_flavor_stage_depends_on_rge_model_not_numerical_runner() -> None:
    source = (RGE_WEINBERG / "FlavorMatchedWeinbergStage.py").read_text(
        encoding="utf-8"
    )

    assert "RGE.running.weinberg.WeinbergRGE" in source
    assert "Numerical.WeinbergRunning" not in source


def test_low_energy_facade_calls_numerical_stage_from_numerical_package() -> None:
    source = (PROJECT_ROOT / "physics" / "LowEnergyNeutrino.py").read_text(
        encoding="utf-8"
    )

    assert "from Numerical.WeinbergStage import" in source
    assert "RGE.running.weinberg.NumericalWeinbergStage" not in source
