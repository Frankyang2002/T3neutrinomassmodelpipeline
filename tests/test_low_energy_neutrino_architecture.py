"""Architecture contracts for the post-matching neutrino calculation.

The central ``pipeline.py`` must show the physical order of the low-energy
calculation, while detailed C5/RGE/file bookkeeping belongs in
``physics/LowEnergyNeutrino.py``.  These tests intentionally avoid numerical
physics changes; they protect the refactor boundary and existing output keys.
"""

from __future__ import annotations

import ast
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
PIPELINE_PATH = PROJECT_ROOT / "pipeline.py"
LOW_ENERGY_PATH = PROJECT_ROOT / "physics" / "LowEnergyNeutrino.py"


def _source(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _tree(path: Path) -> ast.Module:
    return ast.parse(_source(path))


def _function(path: Path, name: str) -> ast.FunctionDef:
    matches = [
        node
        for node in _tree(path).body
        if isinstance(node, ast.FunctionDef) and node.name == name
    ]
    assert len(matches) == 1, f"expected exactly one {name}() in {path.name}"
    return matches[0]


def test_pipeline_imports_named_low_energy_physics_stages() -> None:
    imports = [
        node
        for node in _tree(PIPELINE_PATH).body
        if isinstance(node, ast.ImportFrom)
        and node.module == "physics.LowEnergyNeutrino"
    ]

    assert len(imports) == 1
    assert {alias.name for alias in imports[0].names} == {
        "build_symbolic_neutrino_mass",
        "organise_matched_weinberg_coefficient",
        "run_full_flavor_weinberg_rge",
        "run_numerical_neutrino_observables",
        "run_smeft_weinberg_rge",
    }


def test_pipeline_no_longer_contains_low_energy_implementation_details() -> None:
    pipeline_source = _source(PIPELINE_PATH)
    function_names = {
        node.name
        for node in _tree(PIPELINE_PATH).body
        if isinstance(node, ast.FunctionDef)
    }

    assert "organise_c5_input" not in function_names
    assert "matched_c5_path" not in function_names
    assert "run_matched_weinberg_rge_stage" not in function_names
    assert "run_symbolic_flavor_weinberg_stage" not in function_names
    assert "run_symbolic_neutrino_mass_stage" not in function_names
    assert "run_numerical_weinberg_stage" not in function_names

    assert "RGE.matching.MatchedWeinbergRGE" not in pipeline_source
    assert "RGE.running.weinberg.FlavorMatchedWeinbergStage" not in pipeline_source
    assert "RGE.running.weinberg.NumericalWeinbergStage" not in pipeline_source
    assert "RGE.phenomenology.NeutrinoObservables" not in pipeline_source


def test_pipeline_keeps_low_energy_physics_order_explicit() -> None:
    node = _function(PIPELINE_PATH, "_run_low_energy_neutrino_stages")
    source = ast.get_source_segment(_source(PIPELINE_PATH), node)
    assert source is not None

    ordered_calls = (
        "run_smeft_weinberg_rge",
        "run_full_flavor_weinberg_rge",
        "build_symbolic_neutrino_mass",
        "run_numerical_neutrino_observables",
    )
    positions = [source.index(name) for name in ordered_calls]
    assert positions == sorted(positions)


def test_main_organises_c5_before_intermediate_and_low_energy_stages() -> None:
    node = _function(PIPELINE_PATH, "main")
    source = ast.get_source_segment(_source(PIPELINE_PATH), node)
    assert source is not None

    organise = source.index("organise_matched_weinberg_coefficient")
    intermediate = source.index("_run_intermediate_eft_stages")
    low_energy = source.index("_run_low_energy_neutrino_stages")

    assert organise < intermediate < low_energy


def test_low_energy_module_retains_existing_machine_output_contract() -> None:
    constants = {
        node.value
        for node in ast.walk(_tree(LOW_ENERGY_PATH))
        if isinstance(node, ast.Constant) and isinstance(node.value, str)
    }

    # These summary keys and payload keys are consumed by existing reports and
    # regression tooling.  Refactoring module boundaries must not rename them.
    required_tokens = {
        "WeinbergCoefficientFile",
        "RGEStatus",
        "C5BetaFile",
        "FlavorRGEStatus",
        "C5FlavorBetaMatrixFile",
        "NeutrinoMassStatus",
        "NeutrinoMassMatrixFile",
        "NumericalRGEStatus",
        "NeutrinoMassMatrixLowScaleFile",
        "NeutrinoObservableStatus",
        "NeutrinoObservablesFile",
    }

    assert required_tokens <= constants


def test_low_energy_module_keeps_existing_c5_data_directory_contract() -> None:
    source = _source(LOW_ENERGY_PATH)
    node = _function(LOW_ENERGY_PATH, "organise_matched_weinberg_coefficient")
    function_source = ast.get_source_segment(source, node)
    assert function_source is not None

    assert 'record.output_dir / "data"' in function_source
    assert "Refusing to overwrite existing RGE input" in function_source
    assert 'summary["WeinbergCoefficientFile"]' in function_source
