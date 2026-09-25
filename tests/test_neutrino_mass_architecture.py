"""Architecture contracts for symbolic C5 -> m_nu ownership."""

from __future__ import annotations

from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RGE_STAGE = (
    PROJECT_ROOT
    / "RGE"
    / "running"
    / "weinberg"
    / "FullFlavorWeinbergStage.py"
)
RGE_COMPAT = (
    PROJECT_ROOT
    / "RGE"
    / "running"
    / "weinberg"
    / "FlavorMatchedWeinbergStage.py"
)
MASS_PHYSICS = PROJECT_ROOT / "physics" / "NeutrinoMass.py"
LOW_ENERGY = PROJECT_ROOT / "physics" / "LowEnergyNeutrino.py"


def _source(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def test_full_flavor_rge_stage_contains_only_c5_rge_implementation() -> None:
    source = _source(RGE_STAGE)

    assert "def run_full_flavor_weinberg_stage(" in source
    assert "beta_weinberg_matrix" in source
    assert "physics.NeutrinoMass" not in source
    assert "def build_neutrino_mass_matrix(" not in source
    assert "def run_symbolic_neutrino_mass_stage(" not in source
    assert "neutrino_mass_matrix.txt" not in source


def test_historical_full_flavor_stage_preserves_mass_compatibility_imports() -> None:
    source = _source(RGE_COMPAT)

    assert "from physics.NeutrinoMass import" in source
    assert "from RGE.running.weinberg.FullFlavorWeinbergStage import" in source
    assert "run_flavor_matched_weinberg_rge = run_full_flavor_weinberg_stage" in source
    assert "def build_neutrino_mass_matrix(" not in source
    assert "def run_symbolic_neutrino_mass_stage(" not in source


def test_symbolic_neutrino_mass_physics_owns_c5_to_mass_conversion() -> None:
    source = _source(MASS_PHYSICS)

    assert "def build_neutrino_mass_matrix(" in source
    assert "def run_symbolic_neutrino_mass_stage(" in source
    assert "m_nu = -(v^2/2) C5" in source
    assert "neutrino_mass_matrix.txt" in source
    assert "NeutrinoMassStatus" in source
    assert "beta_weinberg_matrix" not in source


def test_low_energy_facade_uses_physics_neutrino_mass_stage() -> None:
    source = _source(LOW_ENERGY)

    assert "from physics.NeutrinoMass import" in source
    assert (
        "run_symbolic_neutrino_mass_stage as _run_symbolic_neutrino_mass_stage"
        in source
    )
