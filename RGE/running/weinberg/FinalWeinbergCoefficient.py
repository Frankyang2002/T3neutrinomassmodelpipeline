"""Compatibility API for the historical final-Weinberg running path.

The final coefficient is matching/bookkeeping, not an RGE evolution model.
The implementation now lives in ``RGE.matching.FinalWeinbergCoefficient``.
"""

from RGE.matching.FinalWeinbergCoefficient import (
    WEINBERG_KEYS,
    build_final_weinberg_coefficient,
    main,
    normalize_pole_rge_consistency,
)

__all__ = [
    "WEINBERG_KEYS",
    "build_final_weinberg_coefficient",
    "main",
    "normalize_pole_rge_consistency",
]


if __name__ == "__main__":
    raise SystemExit(main())
