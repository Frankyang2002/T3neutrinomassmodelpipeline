"""Architecture contracts for Weinberg flavor-matching names and ownership."""

from __future__ import annotations

from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
MATCHING = PROJECT_ROOT / "RGE" / "matching"
FLAVOR = MATCHING / "WeinbergFlavorMatching.py"
FINAL_ADAPTER = MATCHING / "FinalWeinbergAdapter.py"


def _source(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def test_one_generation_flavor_lifting_has_descriptive_names() -> None:
    source = _source(FLAVOR)

    assert "def extract_one_generation_loop_kernel(" in source
    assert "def build_majorana_c5_flavor_matrix(" in source
    assert "def build_flavor_c5_from_matchete(" in source
    assert "final_weinberg_coefficient.json" not in source


def test_hierarchical_final_c5_adapter_is_separate() -> None:
    source = _source(FINAL_ADAPTER)

    assert "def load_hierarchical_majorana_c5(" in source
    assert "def is_hierarchical_final_c5(" in source
    assert "ready_for_physical_majorana_numerics" in source
    assert "def build_flavor_c5_from_matchete(" not in source


def test_historical_mixed_flavor_module_is_removed() -> None:
    assert not (MATCHING / "FlavorC5Matching.py").exists()


def test_core_callers_use_descriptive_interfaces() -> None:
    callers = (
        PROJECT_ROOT / "RGE" / "running" / "weinberg" / "FullFlavorWeinbergStage.py",
        PROJECT_ROOT / "physics" / "NeutrinoMass.py",
        PROJECT_ROOT / "physics" / "LowEnergyNeutrino.py",
    )

    for path in callers:
        source = _source(path)
        assert "RGE.matching.FlavorC5Matching" not in source
