"""Pure mathematical and report-cache regressions for reverse mixing."""
from __future__ import annotations
import hashlib
import json
import sympy as sp
from RGE.group_factors.validation import C12ReverseMixingGroupFactors as reverse
from Reports import C12GroupFactorPDF as pdf

R = 8 * sp.sqrt(3) / 3


def test_exact_reverse_projection_aligns_full_tensor():
    t = {(1, 1, 5, 11): sp.Rational(1, 4),
         (1, 1, 6, 12): -sp.Rational(1, 4)}
    beta = {k: R * v for k, v in t.items()}
    all_keys = list(t) + [(1, 1, 5, 5)]
    result = reverse.project_reverse_mixing(t, beta, all_llss_components=all_keys)
    assert result['status'] == 'Verified'
    assert result['residual_count'] == 0
    assert sp.simplify(sp.sympify(result['factor']) - R) == 0


def test_outside_matched_direction_fails_verification():
    t = {(1, 1, 5, 11): 1, (1, 1, 6, 12): -1}
    beta = {(1, 1, 5, 11): 2, (1, 1, 6, 12): -2, (1, 1, 5, 5): 1}
    result = reverse.project_reverse_mixing(t, beta, all_llss_components=list(beta))
    assert result['status'] != 'Verified'
    assert result['factor'] is None
    assert result['residual_count'] == 1


def test_absent_or_incomplete_beta_cannot_claim_verified():
    t = {(1, 1, 5, 11): 1, (1, 1, 6, 12): -1}
    assert reverse.project_reverse_mixing(t, {}, all_llss_components=list(t))['status'] != 'Verified'
    partial = reverse.project_reverse_mixing(t, {(1, 1, 5, 11): 3}, all_llss_components=list(t))
    assert partial['residual_count'] == 1


def test_seed_normalisation_uses_matchete_mixed_yukawa_channel():
    u = sp.conjugate(sp.Symbol('y1'))
    v = sp.conjugate(sp.Symbol('y2'))
    mf = sp.Symbol('MF')
    matched = {(1, 1, 5, 11): u*v/(4*mf),
               (1, 1, 6, 12): -u*v/(4*mf)}
    got = reverse._normalised_c12_seed(matched)
    assert got == {(1, 1, 5, 11): sp.Rational(1, 4),
                   (1, 1, 6, 12): -sp.Rational(1, 4)}


def test_cache_accepted_only_when_verified_rep_and_hashes(tmp_path):
    d = tmp_path/'data'; d.mkdir()
    names = (reverse.SEED, reverse.QUARTIC, reverse.RGBETA)
    for name in names:
        (d/name).write_text('unchanged data', encoding='utf-8')
    rep = {'dS1': 1, 'dS2': 3, 'dF': 2, 'alpha': -2}
    result = {'schema': reverse.SCHEMA, 'status': 'Verified',
              'representation': rep, 'factor': '7*sqrt(3)/3',
              'coupling': 'conjugate(lambdaT3)', 'residual_count': 0,
              'source_files_sha256': {name: hashlib.sha256((d/name).read_bytes()).hexdigest()
                                       for name in names}}
    (d/reverse.OUTPUT).write_text(json.dumps(result), encoding='utf-8')
    item = {'model': 'T3-A-m2', 'status': 'Projected', 'representation': rep,
            'terms': {}, 'total_collinear': True, 'data_dir': str(d)}
    assert sp.simplify(pdf._reverse_mixing_factor(item) - 7*sp.sqrt(3)/3) == 0
    assert reverse.validate_cached_projection(d, {**rep, 'alpha': 0}) is None
    (d/reverse.QUARTIC).write_text('changed', encoding='utf-8')
    assert pdf._reverse_mixing_factor(item) is None
