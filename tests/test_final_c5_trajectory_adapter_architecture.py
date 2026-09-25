"""Architecture contracts for the final-C5 numerical trajectory adapter."""

from __future__ import annotations

from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
CANONICAL = PROJECT_ROOT / "Numerical" / "FinalC5TrajectoryAdapter.py"
HISTORICAL = PROJECT_ROOT / "Numerical" / "FinalC5Bridge.py"


def _source(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def test_canonical_module_describes_trajectory_adapter_role() -> None:
    source = _source(CANONICAL)

    assert "class FinalC5EvaluationInputs:" in source
    assert "def extract_final_c5_inputs(" in source
    assert "def evaluate_final_c5_on_trajectory(" in source
    assert "class FinalC5TrajectoryEvaluator:" in source
    assert "from Numerical.SMWeinbergStage import evaluate_final_weinberg_json" in source


def test_historical_bridge_contains_only_compatibility_wrappers() -> None:
    source = _source(HISTORICAL)

    assert "from Numerical.FinalC5TrajectoryAdapter import" in source
    assert "class FinalC5EvaluationInputs:" not in source
    assert "def _majorana_takagi_mass_basis(" not in source
    assert "def extract_final_c5_inputs(" not in source

    # These thin wrappers are intentional because older callers patch symbols
    # on the historical module path.
    assert "def evaluate_final_c5_from_trajectory(" in source
    assert "class FinalC5TrajectoryBuilder:" in source
    assert "from Numerical.SMWeinbergStage import evaluate_final_weinberg_json" in source


def test_historical_symbols_remain_available_from_canonical_module() -> None:
    source = _source(CANONICAL)

    assert "FinalC5NumericalInputs = FinalC5EvaluationInputs" in source
    assert "FinalC5TrajectoryBuilder = FinalC5TrajectoryEvaluator" in source
    assert "final_c5_inputs_from_trajectory = extract_final_c5_inputs" in source
    assert "evaluate_final_c5_from_trajectory = evaluate_final_c5_on_trajectory" in source
