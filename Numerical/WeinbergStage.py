"""Compatibility surface for the historical numerical Weinberg stage name.

The canonical implementation now lives in ``Numerical.SMWeinbergStage``.
"""

from Numerical.SMWeinbergStage import (
    evaluate_final_weinberg_json,
    evaluate_symbolic_c5,
    run_sm_weinberg_numerical_stage,
)

# Historical public function name.
run_numerical_weinberg_stage = run_sm_weinberg_numerical_stage

__all__ = [
    "evaluate_final_weinberg_json",
    "evaluate_symbolic_c5",
    "run_numerical_weinberg_stage",
    "run_sm_weinberg_numerical_stage",
]
