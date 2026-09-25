"""Architecture contracts for matched-Weinberg parsing and benchmark ownership."""

from __future__ import annotations

from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
PARSER = PROJECT_ROOT / "RGE" / "matching" / "MatcheteC5Parsing.py"
BENCHMARK = (
    PROJECT_ROOT
    / "RGE"
    / "running"
    / "weinberg"
    / "OneGenerationWeinbergBenchmark.py"
)
BENCHMARK_STAGE = (
    PROJECT_ROOT
    / "RGE"
    / "running"
    / "weinberg"
    / "OneGenerationWeinbergBenchmarkStage.py"
)
COMPAT = PROJECT_ROOT / "RGE" / "matching" / "MatchedWeinbergRGE.py"
LOW_ENERGY = PROJECT_ROOT / "physics" / "LowEnergyNeutrino.py"


def _source(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def test_matchete_parser_contains_no_rge_or_output_stage() -> None:
    source = _source(PARSER)
    assert "def parse_matchete_c5(" in source
    assert "calculate_complete_wilson_tensor_rge" not in source
    assert "rge_summary.json" not in source


def test_one_generation_benchmark_owns_tensor_cross_check() -> None:
    source = _source(BENCHMARK)
    assert "def higgs_quartic(" in source
    assert "def calculate_one_generation_weinberg_benchmark(" in source
    assert "calculate_complete_wilson_tensor_rge" in source
    assert "json.dumps" not in source
    assert "c5_beta.txt" not in source


def test_benchmark_stage_owns_file_output() -> None:
    source = _source(BENCHMARK_STAGE)
    assert "def run_one_generation_weinberg_benchmark(" in source
    assert "c5_beta.txt" in source
    assert "rge_summary.json" in source
    assert "calculate_complete_wilson_tensor_rge" not in source


def test_low_energy_facade_uses_canonical_benchmark_stage() -> None:
    source = _source(LOW_ENERGY)
    assert (
        "from RGE.running.weinberg.OneGenerationWeinbergBenchmarkStage import"
        in source
    )
    assert "RGE.matching.MatchedWeinbergRGE" not in source


def test_historical_matching_module_uses_lazy_stage_wrapper() -> None:
    source = _source(COMPAT)

    assert "from RGE.matching.MatcheteC5Parsing import parse_matchete_c5" in source
    assert "def run_matched_weinberg_rge(" in source

    prefix = source.split("def run_matched_weinberg_rge(", 1)[0]
    assert "OneGenerationWeinbergBenchmarkStage" not in prefix
    assert (
        "from RGE.running.weinberg.OneGenerationWeinbergBenchmarkStage import"
        in source
    )
    assert "def calculate_matched_weinberg_rge(" not in source
