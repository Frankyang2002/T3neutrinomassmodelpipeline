"""Architecture contracts for the canonical full-flavor Weinberg stage."""

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
LOW_ENERGY = PROJECT_ROOT / "physics" / "LowEnergyNeutrino.py"


def _source(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def test_canonical_full_flavor_stage_has_descriptive_name() -> None:
    source = _source(CANONICAL)

    assert "def run_full_flavor_weinberg_stage(" in source
    assert "def run_flavor_matched_weinberg_rge(" not in source
    assert "beta_weinberg_matrix" in source


def test_historical_full_flavor_stage_is_removed() -> None:
    historical = (
        PROJECT_ROOT
        / "RGE"
        / "running"
        / "weinberg"
        / "FlavorMatchedWeinbergStage.py"
    )
    assert not historical.exists()


def test_low_energy_orchestration_uses_canonical_stage_name() -> None:
    source = _source(LOW_ENERGY)

    assert "from RGE.running.weinberg.FullFlavorWeinbergStage import" in source
    assert "run_full_flavor_weinberg_stage(" in source
    assert "RGE.running.weinberg.FlavorMatchedWeinbergStage" not in source
