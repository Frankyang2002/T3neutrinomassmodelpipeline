"""Pure synthetic parser and decomposition tests (no Wolfram runs)."""
from __future__ import annotations

import numpy as np
import pytest

from Numerical.diagnostics.BetaTermDominance import additive_terms, analyse_beta, category


def test_additive_terms_nested_and_exponents():
    assert additive_terms('g2^2*y1 - 2*Tr[yu . Trans[yu]]*y1 + lambdaT3*y2') == [
        'g2^2*y1', '- 2*Tr[yu . Trans[yu]]*y1', '+ lambdaT3*y2']
    assert additive_terms('1.3e-8*g2 + 2e+3*y1') == ['1.3e-8*g2', '+ 2e+3*y1']


def test_dependency_categories():
    assert category('g2^3') == 'gauge'
    assert category('lambdaT3*lambdaH') == 'scalar'
    assert category('g2^2*y1') == 'mixed'
    assert category('g2^2*y1', target='y1') == 'gauge'
    assert category('y1*Tr[yu . Trans[yu]]') == 'yukawa'


def test_sum_and_fractions_reproduce_canonical():
    env = {'g2': 2., 'y1': 3., 'lambdaT3': 5.}
    expression = 'g2^2 + lambdaT3 - y1'
    result = analyse_beta(expression, env, 6.)
    assert result['full_norm_16pi2_beta'] == pytest.approx(6.)
    assert result['categories']['gauge']['dominance'] == pytest.approx(4./12.)
    assert result['categories']['scalar']['dominance'] == pytest.approx(5./12.)
    assert result['categories']['yukawa']['dominance'] == pytest.approx(3./12.)


def test_detects_mismatched_total():
    with pytest.raises(ValueError, match='canonical'):
        analyse_beta('g2^2 + y1', {'g2': 2., 'y1': 3.}, 9.)
