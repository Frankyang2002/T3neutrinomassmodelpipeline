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
    runner = NUMERICAL / "SMWeinbergEvolution.py"
    stage = NUMERICAL / "SMWeinbergStage.py"

    assert runner.is_file()
    assert stage.is_file()
    assert "solve_ivp" in runner.read_text(encoding="utf-8")
    assert "run_sm_weinberg_numerical_stage" in stage.read_text(encoding="utf-8")


def test_historical_numerical_module_names_are_compatibility_only() -> None:
    runner = (NUMERICAL / "WeinbergRunning.py").read_text(encoding="utf-8")
    stage = (NUMERICAL / "WeinbergStage.py").read_text(encoding="utf-8")

    assert "from Numerical.SMWeinbergEvolution import" in runner
    assert "solve_ivp" not in runner

    assert "from Numerical.SMWeinbergStage import" in stage
    assert "def run_numerical_weinberg_stage(" not in stage


def test_old_mixed_rge_numerical_modules_are_removed() -> None:
    assert not (RGE_WEINBERG / "WeinbergRunning.py").exists()
    assert not (RGE_WEINBERG / "NumericalWeinbergStage.py").exists()


def test_symbolic_flavor_stage_depends_on_rge_model_not_numerical_runner() -> None:
    source = (RGE_WEINBERG / "FullFlavorWeinbergStage.py").read_text(
        encoding="utf-8"
    )

    assert "RGE.running.weinberg.WeinbergRGE" in source
    assert "Numerical.SMWeinbergEvolution" not in source
    assert "Numerical.WeinbergRunning" not in source


def test_historical_flavor_stage_is_compatibility_only() -> None:
    source = (RGE_WEINBERG / "FlavorMatchedWeinbergStage.py").read_text(
        encoding="utf-8"
    )

    assert "from RGE.running.weinberg.FullFlavorWeinbergStage import" in source
    assert "beta_weinberg_matrix" not in source


def test_low_energy_facade_calls_canonical_numerical_stage() -> None:
    source = (PROJECT_ROOT / "physics" / "LowEnergyNeutrino.py").read_text(
        encoding="utf-8"
    )

    assert "from Numerical.SMWeinbergStage import" in source
    assert "_run_sm_weinberg_numerical_stage(" in source
    assert "RGE.running.weinberg.NumericalWeinbergStage" not in source
