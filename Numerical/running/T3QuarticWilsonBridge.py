"""Opt-in bridge between explicit running quartic responses and the seven-term EFT1 RHS.

Does not alter the production solver or infer Matchete/RGBeta normalization.
Only source-backed full response data may be supplied for physics runs.

Project path: Numerical/running/T3QuarticWilsonBridge.py
"""
from __future__ import annotations

from dataclasses import replace
from typing import Mapping
import numpy as np

from Numerical.running.T3RunningQuarticTensor import (
    RunningQuarticBasis, QuarticEvidenceError,
)
from RGE.running.intermediate.GeneralEFT1C5Mixing import (
    Context, TERMS, calculate_master_beta,
)


def running_context(
    initial_context: Context,
    state: Mapping,
    basis: RunningQuarticBasis,
    *,
    expected_hhss=None,
) -> Context:
    """Preserve all Wilson conventions; replace ONLY the quartic tensor.

    `state` must contain instantaneous running RGBeta scalar coupling values.
    `basis` must cover all required quartic coupling responses in the fixed
    original Cartesian scalar basis, including Higgs and both heavy scalars.
    """
    if initial_context.model != basis.model:
        raise QuarticEvidenceError('Model mismatch between Wilson context and response basis')
    q = basis.evaluate(state)
    if q.shape != (initial_context.ns,) * 4:
        raise QuarticEvidenceError('Wilson and running tensor scalar dimensions differ')
    if expected_hhss is not None:
        basis.validate_hhss(state, expected_hhss)
    return replace(initial_context, quartic=q)


def running_seven_term_beta(
    initial_context: Context,
    master_coefficient: np.ndarray,
    state: Mapping,
    basis: RunningQuarticBasis,
    *,
    expected_hhss=None,
    extra_scalar_anomalous_dimension=None,
):
    """Evaluate the existing seven Wilson terms with a running quartic tensor.

    `extra_scalar_anomalous_dimension` optionally supplies an existing known
    (ns,ns) gamma matrix contribution, e.g. the caller's quark-Higgs correction.
    This function does not replace the caller's UV/EFT1 RGE or threshold matching.
    Return convention is 16*pi^2 dC/dln(mu), without multiplying by LOOP.
    """
    ctx = running_context(
        initial_context, state, basis, expected_hhss=expected_hhss)
    c = np.asarray(master_coefficient, dtype=complex)
    shape = (ctx.nf, ctx.nf, ctx.ns, ctx.ns)
    if c.shape != shape or not np.isfinite(c).all():
        raise QuarticEvidenceError('Invalid finite master Wilson tensor')
    terms = calculate_master_beta(ctx, c)
    if set(terms) != set(TERMS):
        raise QuarticEvidenceError('Existing seven-term Wilson inventory changed')
    if extra_scalar_anomalous_dimension is not None:
        gamma = np.asarray(extra_scalar_anomalous_dimension, complex)
        if gamma.shape != (ctx.ns, ctx.ns) or not np.isfinite(gamma).all():
            raise QuarticEvidenceError('Invalid additional scalar anomalous dimension')
        contribution = np.einsum('bd,ijad->ijab', gamma, c, optimize=True)
        terms['scalar_anomalous_dimension'] += (
            contribution + contribution.swapaxes(2, 3))
    beta = sum((terms[k] for k in TERMS), np.zeros(shape, complex))
    if not np.isfinite(beta).all():
        raise QuarticEvidenceError('Nonfinite seven-term Wilson beta')
    return beta, terms
