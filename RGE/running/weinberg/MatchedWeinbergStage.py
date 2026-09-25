"""Compatibility API for the historical one-generation Weinberg stage name.

The canonical stage now lives in
``RGE.running.weinberg.OneGenerationWeinbergBenchmarkStage``.
"""

from RGE.running.weinberg.OneGenerationWeinbergBenchmarkStage import (
    run_one_generation_weinberg_benchmark,
)

# Historical public function name.
run_matched_weinberg_rge = run_one_generation_weinberg_benchmark

__all__ = [
    "run_matched_weinberg_rge",
    "run_one_generation_weinberg_benchmark",
]
