"""Explicit alternative 1-loop beta for the Fonseca--Olgoso--Santiago limit.

JHEP 07 (2026) 094 Eq.(22), arXiv:2512.15866v1: the six channels
excluding the extra crossed-Yukawa contraction. In the gaugeless limit,
these six channels agree numerically with an independently coded Eq.(22).

The full gauge-complete extension used here is the seven-term Ward-covariant
candidate with crossed_yukawa identically zero, retaining the transposed
fermion gauge Casimir. The general gauge-sector formula is NOT separately
certified against the independent publication; nor is the conflicting
September 2026 version of record for arXiv:2512.16890 verified. Thus this
mode is explicitly unverified as a general physical A(-4) RGE.
"""
from __future__ import annotations

import numpy as np

from .GeneralEFT1CrossedWardCandidate import candidate_beta
from .GeneralEFT1SevenTermFlow import Couplings, TERMS, total_beta

MODE = 'fonseca_eq22_six_term'


def six_term_beta(c: np.ndarray, couplings: Couplings) -> dict[str, np.ndarray]:
    """Reuse tested six non-crossed channels; crossed contribution exactly zero."""
    parts = candidate_beta(c, couplings)
    parts['crossed_yukawa'] = np.zeros_like(c, dtype=np.complex128)
    if set(parts) != set(TERMS):
        raise RuntimeError('Wilson beta-term inventory changed')
    return parts


def six_term_beta_total(c: np.ndarray, couplings: Couplings) -> np.ndarray:
    return total_beta(six_term_beta(c, couplings))
