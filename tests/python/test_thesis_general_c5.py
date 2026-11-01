"""Fast offline tests for full-general flavour channel assembly (no Wolfram)."""
import numpy as np
import pytest

from Numerical.matching.GeneralT3FlavourC5 import (
    flavour_c5_from_channels, prepare_heavy_flavour,
)


def row(i, channel, first, second, k=1.0):
    return {'ok': True, 'index': i, 'channel': channel,
            'factors': [first, second], 'real': k, 'imag': 0.0}


def test_all_flavour_channels_and_offdiagonal_structure():
    y1 = np.eye(3, dtype=complex)
    y2 = np.eye(3, dtype=complex)
    extra = np.array([[0, 2, 0], [4, 0, 0], [0, 0, 1]], complex)
    a = 'T3YUKy1Plain'; b = 'T3YUKy2Plain'; x = 'T3YUKy2ConjPlain'
    rows = [row(i, 'primary', a, b) for i in range(3)]
    rows += [row(i, 'additional', a, x, 0.5) for i in range(3)]
    c5 = flavour_c5_from_channels(rows, {'y1': y1, 'y2': y2,
                                          'y2Conj': extra}, ['primary', 'additional'])
    assert np.allclose(c5, c5.T)
    assert np.isclose(c5[0, 1], 0.5 * (extra[0,1] + extra[1,0]) / 2)
    assert not np.allclose(c5, np.eye(3))


def test_yukawa_bar_is_elementwise_conjugation_not_hermitian_transpose():
    y1 = np.eye(3, dtype=complex)
    y2 = 1j * np.eye(3, dtype=complex)
    rows = [row(i, 'phase', 'T3YUKy1Plain', 'T3YUKy2Bar') for i in range(3)]
    got = flavour_c5_from_channels(rows, {'y1':y1,'y2':y2}, ['phase'])
    assert np.allclose(got, -1j*np.eye(3))


def test_missing_heavy_channel_is_rejected():
    with pytest.raises(ValueError, match='grid incomplete'):
        flavour_c5_from_channels([row(0, 'a','T3YUKy1Plain','T3YUKy2Plain')],
                                 {'y1':np.eye(3),'y2':np.eye(3)},['a'])


def test_undeclared_extra_yukawa_cannot_be_silently_discarded():
    rows = [row(i, 'extra','T3YUKy1Plain','T3YUKy2ConjPlain') for i in range(3)]
    with pytest.raises(ValueError, match='missing Yukawa matrix'):
        flavour_c5_from_channels(rows, {'y1':np.eye(3),'y2':np.eye(3)},['extra'])


def test_duplicate_channel_is_rejected():
    rows = [row(0,'a','T3YUKy1Plain','T3YUKy2Plain')]*2
    with pytest.raises(ValueError, match='duplicated'):
        flavour_c5_from_channels(rows, {'y1':np.eye(3),'y2':np.eye(3)})
