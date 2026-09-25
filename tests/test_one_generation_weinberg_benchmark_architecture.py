"""Architecture contracts for the one-generation Weinberg benchmark."""

from __future__ import annotations

from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
MODEL = (
    PROJECT_ROOT
    / "RGE"
    / "running"
    / "weinberg"
    / "OneGenerationWeinbergBenchmark.py"
)
STAGE = (
    PROJECT_ROOT
    / "RGE"
    / "running"
    / "weinberg"
    / "OneGenerationWeinbergBenchmarkStage.py"
)
OLD_MODEL = (
    PROJECT_ROOT
    / "RGE"
    / "running"
    / "weinberg"
    / "MatchedWeinbergRGE.py"
)
OLD_STAGE = (
    PROJECT_ROOT
    / "RGE"
    / "running"
    / "weinberg"
    / "MatchedWeinbergStage.py"
)
LOW_ENERGY = PROJECT_ROOT / "physics" / "LowEnergyNeutrino.py"


def _source(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def test_benchmark_model_has_descriptive_name_and_no_io() -> None:
    source = _source(MODEL)

    assert "def calculate_one_generation_weinberg_benchmark(" in source
    assert "calculate_complete_wilson_tensor_rge" in source
    assert "json.dumps" not in source
    assert "c5_beta.txt" not in source


def test_benchmark_stage_owns_output_contract() -> None:
    source = _source(STAGE)

    assert "def run_one_generation_weinberg_benchmark(" in source
    assert "c5_beta.txt" in source
    assert "rge_summary.json" in source
    assert "calculate_complete_wilson_tensor_rge" not in source


def test_historical_running_names_are_compatibility_only() -> None:
    model_source = _source(OLD_MODEL)
    stage_source = _source(OLD_STAGE)

    assert "from RGE.running.weinberg.OneGenerationWeinbergBenchmark import" in model_source
    assert "def calculate_matched_weinberg_rge(" not in model_source

    assert (
        "from RGE.running.weinberg.OneGenerationWeinbergBenchmarkStage import"
        in stage_source
    )
    assert "def run_matched_weinberg_rge(" not in stage_source


def test_low_energy_orchestration_uses_benchmark_name_directly() -> None:
    source = _source(LOW_ENERGY)

    assert (
        "from RGE.running.weinberg.OneGenerationWeinbergBenchmarkStage import"
        in source
    )
    assert "run_one_generation_weinberg_benchmark(" in source
    assert "RGE.running.weinberg.MatchedWeinbergStage" not in source
