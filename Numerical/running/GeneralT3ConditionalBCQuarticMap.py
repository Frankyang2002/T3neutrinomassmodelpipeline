"""Conditional B(-1)/C(-1) Matchete-to-RGBeta HHSS tensor sign map.

The existing one-scale B/C evaluator assumes
    lambda_Matchete = -lambda_RGBeta
for the two direct general-scalar quartics below.  The scalar-pair EFT1
contraction must use the SAME conditional signs when it is compared against
that one-scale anchor.

This operation changes ONLY the two already-constructed HHSS tensor channels;
it does NOT change running scalar couplings, RGEs, masses, tree LLSS Wilson
sources, or the physical one-scale Matchete evaluator.  The map is a candidate
convention, NOT a certified Matchete/RGBeta basis identification.
"""
from __future__ import annotations

from dataclasses import replace
from typing import Any

import numpy as np

from Numerical.running.GeneralT3SequentialLeadingLog import TensorBundle, require

MIXED = frozenset(('B_-1', 'C_-1'))
DIRECT_QUARTICS = (
    'lambdaHHdagS2S2',
    'lambdaHHdagS1barS1bar',
)
CHANNEL_NAMES = tuple('general:' + key for key in DIRECT_QUARTICS)
CONVENTION = 'conditional lambda_Matchete = -lambda_RGBeta (two direct B/C HHSS quartics)'


def apply_conditional_bc_direct_quartic_map(model: str,
                                             bundle: TensorBundle) -> TensorBundle:
    """Return a NEW tensor bundle with both B/C direct quartic signs inverted.

    Fail closed on the wrong model, missing/duplicated channels, or a previously
    applied map.  Retain the other native/general/portal tensor contributions
    and the original mass eigenbasis exactly.
    """
    require(model in MIXED and bundle.model == model,
            'conditional direct quartic map is restricted to B_-1/C_-1')
    names = tuple(bundle.names)
    require(len(names) == len(set(names)),
            'duplicate HHSS channels before conditional sign map')
    require(all(names.count(name) == 1 for name in CHANNEL_NAMES),
            'required B/C direct quartic channels are missing or duplicated')
    provenance = dict(bundle.provenance)
    require('bc_direct_quartic_convention' not in provenance,
            'B/C direct quartic map already applied')

    old = np.asarray(bundle.channels, dtype=float)
    require(old.ndim == 5 and old.shape[0] == len(names) and
            old.shape[1:3] == (4, 4) and old.shape[3] == old.shape[4] and
            np.isfinite(old).all(), 'invalid HHSS channel tensor')
    original_total = np.asarray(bundle.quartic, dtype=float)
    require(original_total.shape == old.shape[1:] and
            np.isfinite(original_total).all() and
            np.allclose(original_total, old.sum(axis=0), rtol=1e-10, atol=1e-12),
            'initial quartic total is not the sum of its channels')

    updated = old.copy()
    for name in CHANNEL_NAMES:
        updated[names.index(name)] *= -1.0
    total = updated.sum(axis=0)
    require(np.allclose(total, total.transpose(1, 0, 2, 3),
                        rtol=1e-10, atol=1e-12) and
            np.allclose(total, total.transpose(0, 1, 3, 2),
                        rtol=1e-10, atol=1e-12),
            'conditional HHSS quartic violates scalar Bose symmetry')

    descriptions = []
    for item in provenance.get('quartic_channels', []):
        row = dict(item)
        if row.get('name') in DIRECT_QUARTICS:
            require(row.get('origin') != 'native Matchete source CG',
                    'refusing to sign-flip native source quartic')
            row['conditional_matchete_over_rgbeta_sign'] = -1
        descriptions.append(row)
    provenance.update({
        'quartic_channels': descriptions,
        'bc_direct_quartic_convention': CONVENTION,
        'bc_direct_quartics_signed': list(DIRECT_QUARTICS),
        'bc_direct_quartic_channel_signs': {name: -1 for name in CHANNEL_NAMES},
        'bc_direct_quartic_map_stage': 'EFT1 HHSS tensor construction only',
        'bc_direct_quartic_map_affects_rges': False,
        'bc_direct_quartic_map_physically_certified': False,
        'matchete_rgbeta_map_certified': False,
    })
    return replace(bundle, channels=updated, quartic=total, provenance=provenance)
