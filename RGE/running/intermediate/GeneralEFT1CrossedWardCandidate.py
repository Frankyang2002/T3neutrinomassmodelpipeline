"""Opt-in Ward- and unitary-basis-covariant RGE *candidate*; NOT loop certified.

Use alongside GeneralEFT1SevenTermFlow, without modifying the upstream
symbolic (arXiv:2512.16890v1 Eq. 4.85) or fixed-order evaluator.

Two distinct tensor-index convention problems are exposed by covariance tests:

1. ``y_jla y_ikc C_lkbc`` in the printed crossed term is not covariant under
   arbitrary complex unitary changes of left-Weyl basis. ``C*`` is required
   for this index contraction to have the right transformation law. The
   numerical candidate below applies the *four* separate pair exchanges.
2. With coefficient indices defined by L5=C_ijab psi_i psi_j phi_a phi_b,
   the external-fermion contraction of the gauge wavefunction contribution
   requires ``C2(F).T``. The repository's ``gamma_f = -3g^2 C2(F) + ...``
   agrees in the canonical block basis where C2 is real/diagonal, but is not
   generally covariant under complex unitary rotations of reducible bases.

These corrections are mathematically compelled as covariance constraints for
this index convention; covariance does NOT independently fix loop coefficients
or certify the revised diagrammatic expression. No physical A(-4) result is
certified until a loop-level source/counterterm calculation resolves the
unstarred-C discrepancy in arXiv:2512.16890v1.
"""
from __future__ import annotations

import numpy as np
from RGE.running.intermediate.GeneralEFT1SevenTermFlow import (
    Couplings, TERMS, master_beta, total_beta,
)

CROSSED = "crossed_yukawa"
FERMION = "fermion_anomalous_dimension"
CANDIDATE_STATUS = "UNVERIFIED_WARD_AND_UNITARY_COVARIANT_RGE_CANDIDATE"


def crossed_candidate(c: np.ndarray, inputs: Couplings) -> np.ndarray:
    """Anti-linear crossed term: four independent external pair exchanges."""
    y = inputs.yukawa
    if c.shape != (y.shape[0], y.shape[0], y.shape[2], y.shape[2]):
        raise ValueError("Wilson/coupling dimensions do not agree")
    if not np.isfinite(c).all():
        raise ValueError("Non-finite Wilson entries")
    # Two successive contractions avoid constructing a dense anomalous matrix.
    base = np.einsum("jla,ikc,lkbc->ijab", y, y, c.conj(), optimize=True)
    return (base + base.swapaxes(0, 1) + base.swapaxes(2, 3)
            + base.transpose(1, 0, 3, 2))


def fermion_casimir_transpose_delta(c: np.ndarray, inputs: Couplings) -> np.ndarray:
    """Correct only the gauge Casimir orientation in the fermion gamma term.

    delta_gamma_f[j,k] = -3 sum_g g^2 [C2(F)^T - C2(F)]_{jk}.
    This vanishes on canonical irrep blocks with real diagonal Casimirs.
    """
    nf = c.shape[0]
    delta = np.zeros((nf, nf), dtype=complex)
    for g, pairs in inputs.gauge:
        for tf, _ in pairs:
            casimir = tf @ tf
            delta -= 3 * g*g * (casimir.T-casimir)
    base = np.einsum("jk,ikab->ijab", delta, c, optimize=True)
    return base + base.swapaxes(0, 1)


def candidate_beta(c: np.ndarray, inputs: Couplings) -> dict[str, np.ndarray]:
    """Full seven-term beta with two covariance-only candidate replacements."""
    parts = master_beta(c, inputs, corrected_crossed=True)
    parts[CROSSED] = crossed_candidate(c, inputs)
    parts[FERMION] = parts[FERMION] + fermion_casimir_transpose_delta(c, inputs)
    if set(parts) != set(TERMS):
        raise RuntimeError("Seven-term inventory changed unexpectedly")
    return parts


def candidate_beta_total(c: np.ndarray, inputs: Couplings) -> np.ndarray:
    return total_beta(candidate_beta(c, inputs))
