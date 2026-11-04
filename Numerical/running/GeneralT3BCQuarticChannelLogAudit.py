"""Read-only B(-1)/C(-1) Matchete-vs-EFT1 quartic-channel log audit.

For each of the three one-scale Matchete quartics, hold fitted masses, Yukawas
and scalar mixing fixed, zero the other two quartics in an endpoint COPY, and
shift every Log[MF^2/MS_{X,Y}^2] by log(2).  The difference / log(2) is the
independently extracted common-log coefficient for that quartic channel.

Compare with the separately saved, signed EFT1 per-quartic leading log divided
by 2 log(muF/muS).  Verify that BOTH independent three-channel sums recover
their total and that the original one-scale coefficient is unchanged.

No arbitrary normalization is applied.  Fitted factors are diagnostics only.
Neither the conditional Matchete/RGBeta sign map nor the GammaCC/Weyl factor is
physically certified.  No fitted files or RGE inputs are altered.
"""
from __future__ import annotations

from copy import deepcopy
import math
import re
from typing import Any, Mapping

import numpy as np

MODELS = ('B_-1', 'C_-1')
CHANNELS = (
    ('lambdaT3', 'lambdaT3', 'native:lambdaT3'),
    ('lambdaHHdagS2S2', 'lambdaT3ExtraQ020020I1',
     'general:lambdaHHdagS2S2'),
    ('lambdaHHdagS1barS1bar', 'lambdaT3ExtraQ020200I1',
     'general:lambdaHHdagS1barS1bar'),
)
FORMAT = 't3-bc-independent-three-quartic-log-audit-v1'
_DIRECT = frozenset((CHANNELS[1][0], CHANNELS[2][0]))
CHANNEL_BY_EFT = {eft: key for key, _symbol, eft in CHANNELS}


def require(ok: Any, message: str) -> None:
    if not bool(ok):
        raise ValueError('BC quartic-channel audit: ' + message)


def validate_source_inventory(source: str) -> dict[str, int]:
    """Fail closed if the one-scale source's scalar coupling set has changed."""
    require(isinstance(source, str) and source.strip(),
            'Matchete coefficient source is missing')
    found: dict[str, int] = {}
    for name in re.findall(r'Coupling\[(lambda[A-Za-z0-9]+),\s*\{\},\s*0\]', source):
        found[name] = found.get(name, 0) + 1
    required = {symbol for _key, symbol, _eft in CHANNELS}
    require(set(found) == required and all(found[s] > 0 for s in required),
            f'expected precisely three Matchete quartics {sorted(required)}, '
            f'observed {sorted(found)}')
    return found


def masked_endpoint(endpoint: Mapping[str, Any], active: str) -> dict:
    """Isolate one quartic in a deep copy; never touch running state in place."""
    require(active in {item[0] for item in CHANNELS},
            'unknown one-scale scalar quartic')
    require(isinstance(endpoint, Mapping) and
            'lambdaT3' in endpoint and isinstance(endpoint.get('reused_quartics'), Mapping),
            'saved EFT1 endpoint lacks B/C direct quartic inventory')
    require(all(name in endpoint['reused_quartics'] for name in _DIRECT),
            'missing a B/C direct quartic in EFT1 endpoint')
    out = deepcopy(dict(endpoint))
    zero = {'real': 0.0, 'imag': 0.0}
    if active != 'lambdaT3':
        out['lambdaT3'] = dict(zero)
    for name in _DIRECT:
        if active != name:
            out['reused_quartics'][name] = dict(zero)
    return out


def relative_error(a: np.ndarray, reference: np.ndarray) -> float:
    a = np.asarray(a, complex)
    ref = np.asarray(reference, complex)
    require(a.shape == ref.shape == (3, 3) and np.isfinite(a).all() and
            np.isfinite(ref).all(), 'nonfinite or wrong-shaped C5')
    return float(np.linalg.norm(a-ref) / max(float(np.linalg.norm(ref)), 1e-38))


def eft_quartic_channels(case: Mapping[str, Any]) -> tuple[dict[str, np.ndarray], np.ndarray, dict]:
    """Read signed, *saved* EFT1 decomposition; never recompute or fit it."""
    from Numerical.running.GeneralT3MatcheteLogCoefficientAudit import decode_matrix
    require(case.get('case') == '1111', 'need U=E=W=S=1 saved ablation case')
    p = case.get('quartic_provenance')
    require(isinstance(p, Mapping), 'missing saved quartic provenance')
    signs = p.get('bc_direct_quartic_channel_signs')
    require(isinstance(signs, Mapping) and all(
            signs.get('general:'+name) == -1 for name in _DIRECT),
            'saved case lacks the TWO required conditional B/C sign flips')
    require(p.get('bc_direct_quartic_map_physically_certified') is False,
            'unexpectedly certified B/C conditional sign map')
    entries = case.get('per_general_interaction')
    require(isinstance(entries, list), 'missing EFT1 per-channel decomposition')
    rows = [x for x in entries if x.get('category') == 'quartic_decomposition']
    counts = [x.get('name') for x in rows]
    require(len(counts) == len(set(counts)), 'duplicate EFT1 quartic channel')
    require(all(counts.count(channel) == 1 for channel in CHANNEL_BY_EFT),
            'missing or duplicate B/C quartic in EFT1 decomposition')
    log_mu = math.log(1e5/1e3)
    channels = {name: decode_matrix(next(x['c5_leading_log'] for x in rows
                       if x['name'] == name))/(2*log_mu)
                for name in CHANNEL_BY_EFT}
    total = decode_matrix(case['c5_induced_wilson_log'])/(2*log_mu)
    sum_of_all = sum((decode_matrix(x['c5_leading_log'])/(2*log_mu)
                      for x in rows), np.zeros((3, 3), complex))
    all_residual = relative_error(sum_of_all, total)
    relevant_sum = sum(channels.values(), np.zeros((3, 3), complex))
    three_residual = relative_error(relevant_sum, total)
    require(all_residual < 1e-8,
            f'EFT1 quartic decomposition does not close (error={all_residual:.3g})')
    require(three_residual < 1e-8,
            f'additional nonzero EFT1 HHSS channels detected (error={three_residual:.3g})')
    return channels, total, {
        'all_eft1_quartics_sum_relative_error': all_residual,
        'three_selected_eft1_quartics_sum_relative_error': three_residual,
        'all_eft1_channel_names': counts,
        'direct_quartics_have_conditional_negative_sign': True,
        'spinor_normalization_certified': False,
    }


def run_channel_comparison(model: str, *, case: dict, endpoint: dict,
                           heavy: tuple, coefficient_text: str,
                           original_compiled: Any, shifted_compiled: Any,
                           prior_matchete: np.ndarray,
                           hierarchy_factors: tuple[int, ...] = (1,)) -> dict:
    """Run independent one-scale channel toggles and saved-EFT1 comparisons."""
    from scripts import RunT3GeneralScalarBCExploratoryFit as bc
    from Numerical.running.GeneralT3MatcheteLogCoefficientAudit import (
        compare_coefficients, encode_matrix, decode_matrix,
    )
    require(model in MODELS, 'only B_-1 and C_-1 are supported')
    require(hierarchy_factors and hierarchy_factors[0] == 1 and
            len(set(hierarchy_factors)) == len(hierarchy_factors) and
            all(type(x) is int and x >= 1 for x in hierarchy_factors),
            'hierarchy factors must be unique positive integers beginning with 1')
    symbols = validate_source_inventory(coefficient_text)
    eft, eft_total, eft_audit = eft_quartic_channels(case)
    expected_one_scale = decode_matrix(case['c5_matching_anchored_reference'])
    prior = np.asarray(prior_matchete, complex)
    require(prior.shape == (3, 3) and np.isfinite(prior).all(),
            'prior independent Matchete derivative must be a finite 3x3 matrix')
    scans = []
    for hierarchy in hierarchy_factors:
        mass = np.asarray(heavy[0], float)*hierarchy
        current_heavy = (mass, heavy[1], heavy[2], heavy[3])
        full_original = bc.c5_at_endpoint(model,endpoint,current_heavy,original_compiled)
        full_shifted = bc.c5_at_endpoint(model,endpoint,current_heavy,shifted_compiled)
        full_derivative = (full_shifted-full_original)/math.log(2)
        one_scale_by_channel: dict[str, np.ndarray] = {}
        derivatives: dict[str, np.ndarray] = {}
        output_channels = []
        for key, symbol, eft_name in CHANNELS:
            masked = masked_endpoint(endpoint,key)
            original = bc.c5_at_endpoint(model,masked,current_heavy,original_compiled)
            shifted = bc.c5_at_endpoint(model,masked,current_heavy,shifted_compiled)
            matchete = (shifted-original)/math.log(2)
            eft1 = eft[eft_name]/hierarchy
            comparison = compare_coefficients(matchete, eft1)
            one_scale_by_channel[key] = original
            derivatives[key] = matchete
            output_channels.append({
                'quartic': key,
                'matchete_symbol': symbol,
                'eft1_channel': eft_name,
                'matchete_common_log_coefficient': encode_matrix(matchete),
                'eft1_common_log_coefficient': encode_matrix(eft1),
                'comparison': comparison,
            })
        sum_original = sum(one_scale_by_channel.values(), np.zeros((3,3), complex))
        sum_derivative = sum(derivatives.values(), np.zeros((3,3), complex))
        closure = {
            'matchete_original_three_channel_closure': relative_error(sum_original, full_original),
            'matchete_log_three_channel_closure': relative_error(sum_derivative, full_derivative),
            'matchete_vs_prior_independent_log': (relative_error(full_derivative,prior)
                                                 if hierarchy == 1 else None),
            'matchete_original_vs_signed_ablation_anchor': (
                relative_error(full_original,expected_one_scale)
                if hierarchy == 1 else None),
            'eft1_log_three_channel_closure': eft_audit['three_selected_eft1_quartics_sum_relative_error'],
        }
        require(closure['matchete_original_three_channel_closure'] < 1e-8 and
                closure['matchete_log_three_channel_closure'] < 1e-8,
                f'{model}: three Matchete quartics do not reproduce full one-scale C5')
        if hierarchy == 1:
            require(closure['matchete_vs_prior_independent_log'] < 1e-5 and
                    closure['matchete_original_vs_signed_ablation_anchor'] < 1e-5,
                    f'{model}: full one-scale result differs from saved signed audit')
        full_compare = compare_coefficients(full_derivative,eft_total/hierarchy)
        scans.append({
            'heavy_mass_scale_factor': hierarchy,
            'channels': output_channels,
            'total_comparison': full_compare,
            'checks': closure,
        })
    return {
        'format': FORMAT,
        'model': model,
        'status': 'Computed',
        'coupling_occurrences_in_one_scale_source': symbols,
        'hierarchy_scan': scans,
        'eft1_source_channel_diagnostics': eft_audit,
        'matchete_direct_quartic_map': 'lambda_Matchete=-lambda_RGBeta (conditional)',
        'normalization_factors_are_diagnostics_only': True,
        'independent_hard_matching_computed': False,
        'spinor_and_quartic_mapping_certified': False,
    }
