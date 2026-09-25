"""Architecture contracts for final-Weinberg coefficient ownership."""

from __future__ import annotations

from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
MATCHING = (
    PROJECT_ROOT
    / "RGE"
    / "matching"
    / "FinalWeinbergCoefficient.py"
)
HISTORICAL = (
    PROJECT_ROOT
    / "RGE"
    / "running"
    / "weinberg"
    / "FinalWeinbergCoefficient.py"
)


def _source(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def test_final_weinberg_implementation_lives_under_matching() -> None:
    source = _source(MATCHING)

    assert "def build_final_weinberg_coefficient(" in source
    assert "def normalize_pole_rge_consistency(" in source
    assert "def _build_msbar_hard_at_mf(" in source
    assert "def _build_physical_majorana_hard(" in source
    assert "fixed_order_one_loop_MSbar_hierarchical" in source


def test_historical_running_module_is_compatibility_only() -> None:
    source = _source(HISTORICAL)

    assert "from RGE.matching.FinalWeinbergCoefficient import" in source
    assert "def build_final_weinberg_coefficient(" not in source
    assert "def normalize_pole_rge_consistency(" not in source
    assert "def _build_msbar_hard_at_mf(" not in source
    assert "def _build_physical_majorana_hard(" not in source


def test_compatibility_path_preserves_cli_entrypoint() -> None:
    source = _source(HISTORICAL)

    assert '"main"' in source
    assert 'if __name__ == "__main__":' in source
    assert "raise SystemExit(main())" in source
