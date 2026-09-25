"""Compatibility API for the historical one-generation Weinberg RGE name.

The canonical implementation now lives in
``RGE.running.weinberg.OneGenerationWeinbergBenchmark``.
"""

from RGE.running.weinberg.OneGenerationWeinbergBenchmark import (
    calculate_one_generation_weinberg_benchmark,
    higgs_quartic,
    lambdaH,
)

# Historical public function name.
calculate_matched_weinberg_rge = calculate_one_generation_weinberg_benchmark

__all__ = [
    "calculate_matched_weinberg_rge",
    "calculate_one_generation_weinberg_benchmark",
    "higgs_quartic",
    "lambdaH",
]
