"""Compatibility API for the historical matched-Weinberg module path.

Parsing belongs to ``RGE.matching.MatcheteC5Parsing``. The one-generation
SMEFT cross-check is now explicitly named as a benchmark under
``RGE.running.weinberg``.
"""

from __future__ import annotations

from pathlib import Path

from RGE.matching.MatcheteC5Parsing import parse_matchete_c5
from RGE.running.weinberg.OneGenerationWeinbergBenchmark import (
    calculate_one_generation_weinberg_benchmark,
    higgs_quartic,
    lambdaH,
)


# Historical public function name.
calculate_matched_weinberg_rge = calculate_one_generation_weinberg_benchmark


def run_matched_weinberg_rge(
    c5_path: Path,
    output_dir: Path,
    *,
    debug_outputs: bool = False,
) -> dict:
    """Compatibility wrapper around the renamed benchmark stage."""

    from RGE.running.weinberg.OneGenerationWeinbergBenchmarkStage import (
        run_one_generation_weinberg_benchmark,
    )

    return run_one_generation_weinberg_benchmark(
        c5_path=c5_path,
        output_dir=output_dir,
        debug_outputs=debug_outputs,
    )


__all__ = [
    "calculate_matched_weinberg_rge",
    "calculate_one_generation_weinberg_benchmark",
    "higgs_quartic",
    "lambdaH",
    "parse_matchete_c5",
    "run_matched_weinberg_rge",
]
