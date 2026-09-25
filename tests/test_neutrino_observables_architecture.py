"""Architecture contracts for neutrino-observable ownership."""

from __future__ import annotations

import ast
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
PHYSICS_PATH = PROJECT_ROOT / "physics" / "NeutrinoObservables.py"
LOW_ENERGY_PATH = PROJECT_ROOT / "physics" / "LowEnergyNeutrino.py"
COMPAT_PATH = (
    PROJECT_ROOT
    / "RGE"
    / "phenomenology"
    / "NeutrinoObservables.py"
)


def _source(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def test_neutrino_observables_live_under_physics() -> None:
    source = _source(PHYSICS_PATH)

    assert "class NeutrinoObservables" in source
    assert "def takagi_factorization(" in source
    assert "def calculate_neutrino_observables(" in source
    assert "def run_neutrino_observables_stage(" in source

    assert "neutrino_observables.json" in source
    assert "NeutrinoObservableStatus" in source
    assert "NeutrinoObservablesFile" in source


def test_low_energy_facade_imports_observables_from_physics() -> None:
    tree = ast.parse(_source(LOW_ENERGY_PATH))
    imports = [
        node
        for node in tree.body
        if isinstance(node, ast.ImportFrom)
        and node.module == "physics.NeutrinoObservables"
    ]

    assert len(imports) == 1
    assert {alias.name for alias in imports[0].names} == {
        "run_neutrino_observables_stage"
    }

    assert "RGE.phenomenology.NeutrinoObservables" not in _source(
        LOW_ENERGY_PATH
    )


def test_historical_rge_phenomenology_path_is_compatibility_only() -> None:
    source = _source(COMPAT_PATH)

    assert "from physics.NeutrinoObservables import" in source
    assert "def takagi_factorization(" not in source
    assert "def calculate_neutrino_observables(" not in source
    assert "def run_neutrino_observables_stage(" not in source
