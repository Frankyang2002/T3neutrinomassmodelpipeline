"""Architecture contracts for the full-flavor Weinberg stage rename."""

from __future__ import annotations

from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
CANONICAL = (
    PROJECT_ROOT
    / "RGE"
    / "running"
    / "weinberg"
    / "FullFlavorWeinbergStage.py"
)
HISTORICAL = (
    PROJECT_ROOT
    / "RGE"
    / "running"
    / "weinberg"
    / "FlavorMatchedWeinbergStage.py"
)
LOW_ENERGY = PROJECT_ROOT / "physics" / "LowEnergyNeutrino.py"


def _source(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def test_canonical_full_flavor_stage_has_descriptive_name() -> None:
    source = _source(CANONICAL)

    assert "def run_full_flavor_weinberg_stage(" in source
    assert "def run_flavor_matched_weinberg_rge(" not in source
    assert "beta_weinberg_matrix" in source


def test_historical_stage_is_compatibility_only() -> None:
    source = _source(HISTORICAL)

    assert "from RGE.running.weinberg.FullFlavorWeinbergStage import" in source
    assert "run_flavor_matched_weinberg_rge = run_full_flavor_weinberg_stage" in source
    assert "beta_weinberg_matrix" not in source


def test_low_energy_orchestration_uses_canonical_stage_name() -> None:
    source = _source(LOW_ENERGY)

    assert "from RGE.running.weinberg.FullFlavorWeinbergStage import" in source
    assert "run_full_flavor_weinberg_stage(" in source
    assert "RGE.running.weinberg.FlavorMatchedWeinbergStage" not in source
