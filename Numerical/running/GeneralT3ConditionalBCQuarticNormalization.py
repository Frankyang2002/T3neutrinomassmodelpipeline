"""Conditional B(-1)/C(-1) channel normalization at the EFT1 HHSS tensor.

The *independent* common-mass-log audit supports the FIXED channel factors
    native:lambdaT3                -> sqrt(3)
    general:lambdaHHdagS2S2        -> 2
    general:lambdaHHdagS1barS1bar  -> 2
AFTER the existing sign-map operation flips the two direct RGBeta quartics.
The factors are numerical convention candidates, NOT a source-level proof of
Matchete/RGBeta or GammaCC/Weyl physical normalization.  No fitted factor is
used, and no UV, EFT1 scalar or SMEFT beta-function is changed.

UV -> integrate out ALL fermions -> EFT1(both scalars)
   -> integrate out ALL scalars -> SMEFT.
"""
from __future__ import annotations

from dataclasses import replace
import math

import numpy as np

from Numerical.running.GeneralT3SequentialLeadingLog import TensorBundle, require

MIXED = frozenset(('B_-1', 'C_-1'))
CHANNEL_FACTORS = {
    'native:lambdaT3': math.sqrt(3.0),
    'general:lambdaHHdagS2S2': 2.0,
    'general:lambdaHHdagS1barS1bar': 2.0,
}
DIRECT_CHANNELS = frozenset(CHANNEL_FACTORS) - {'native:lambdaT3'}
CONVENTION = 'conditional B/C common-log channel factors sqrt(3), 2, 2 after signed direct quartics'


def apply_conditional_bc_channel_normalization(model: str,
                                                bundle: TensorBundle) -> TensorBundle:
    """Return a new bundle with only the three fixed channel normalizations.

    Requires the existing B/C signed direct-quartic map to have already run.
    In particular, the direct channel factors are POSITIVE 2 here, not -2:
    the signs are owned by GeneralT3ConditionalBCQuarticMap.
    """
    from Numerical.running.GeneralT3ConditionalBCQuarticMap import (
        CHANNEL_NAMES, DIRECT_QUARTICS,
    )

    require(model in MIXED and bundle.model == model,
            'conditional channel normalization is restricted to B_-1/C_-1')
    names = tuple(bundle.names)
    require(len(names) == len(set(names)) and
            all(names.count(name) == 1 for name in CHANNEL_FACTORS),
            'missing or duplicated B/C normalization channel')
    provenance = dict(bundle.provenance)
    require('bc_channel_normalization_factors' not in provenance,
            'B/C channel normalization has already been applied')
    require(provenance.get('bc_direct_quartic_map_stage') ==
            'EFT1 HHSS tensor construction only' and
            provenance.get('bc_direct_quartic_channel_signs') ==
            {name: -1 for name in CHANNEL_NAMES} and
            provenance.get('bc_direct_quartics_signed') == list(DIRECT_QUARTICS),
            'conditional direct-quartic sign map must run first')
    old = np.asarray(bundle.channels, dtype=float)
    total_old = np.asarray(bundle.quartic, dtype=float)
    require(old.ndim == 5 and old.shape[0] == len(names) and
            old.shape[1:3] == (4, 4) and old.shape[3] == old.shape[4] and
            np.isfinite(old).all() and total_old.shape == old.shape[1:] and
            np.isfinite(total_old).all(), 'invalid B/C quartic tensor')
    require(np.allclose(total_old, old.sum(axis=0), rtol=2e-12, atol=1e-14),
            'quartic total and channel sum disagree before normalization')

    new = old.copy()
    for name, factor in CHANNEL_FACTORS.items():
        new[names.index(name)] *= factor
    total = new.sum(axis=0)
    require(np.allclose(total, total.transpose(1, 0, 2, 3),
                        rtol=1e-10, atol=1e-12) and
            np.allclose(total, total.transpose(0, 1, 3, 2),
                        rtol=1e-10, atol=1e-12),
            'normalized B/C quartic violates Bose symmetry')

    descriptions = []
    matched = set()
    for old_desc in provenance.get('quartic_channels', []):
        desc = dict(old_desc)
        source_name = desc.get('name')
        full_name = ('native:' if desc.get('origin') ==
                     'native Matchete source CG' else 'general:') + str(source_name)
        if full_name in CHANNEL_FACTORS:
            require(full_name not in matched,
                    f'duplicate normalization metadata for {full_name}')
            matched.add(full_name)
            if full_name in DIRECT_CHANNELS:
                require(desc.get('conditional_matchete_over_rgbeta_sign') == -1,
                        'direct quartic sign metadata missing')
            desc['conditional_log_normalization_factor'] = CHANNEL_FACTORS[full_name]
        descriptions.append(desc)
    require(matched == set(CHANNEL_FACTORS),
            'missing B/C normalization channel metadata')
    provenance.update({
        'quartic_channels': descriptions,
        'bc_channel_normalization_convention': CONVENTION,
        'bc_channel_normalization_factors': dict(CHANNEL_FACTORS),
        'bc_channel_normalization_stage': 'EFT1 HHSS tensor construction only',
        'bc_channel_normalization_affects_rges': False,
        'bc_channel_normalization_factors_fitted': False,
        'bc_channel_normalization_physical_certification': False,
        'matchete_rgbeta_map_certified': False,
        'spinor_normalization_certified': False,
    })
    return replace(bundle, channels=new, quartic=total, provenance=provenance)
