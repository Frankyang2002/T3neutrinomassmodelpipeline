"""Offline structural tests; Matchete live evaluation occurs on Windows host."""
import numpy as np
import pytest

from scripts.RunT3Aminus2YukawaBenchmark import (
    complex_matrix, parse_fractions, params_from_config, spectral_masses,
)


def test_fractions():
    assert parse_fractions('0,0.05,0.1') == (0.0, 0.05, 0.1)
    for value in ('', '-0.1', '1.2', '0,0', 'potato'):
        with pytest.raises(ValueError):
            parse_fractions(value)


def test_reference_and_entry_matrices():
    real = [[1,2,3],[4,5,6],[7,8,9]]
    config = {'base_state': {'ordinary': {
        'y1': {'real': real, 'imag': [[0]*3 for _ in range(3)]},
        'y2': {'real': real, 'imag': [[0]*3 for _ in range(3)]}
    }}}
    params = params_from_config(config)
    assert params['y1_00'] == 1 and params['y2_22'] == 9
    entries = [[{'real': i*3+j, 'imag': 0} for j in range(3)] for i in range(3)]
    assert np.array_equal(complex_matrix(entries).real, np.arange(9).reshape(3,3))
    assert np.array_equal(complex_matrix(config['base_state']['ordinary']['y1']).real,
                          np.array(real))


def test_bilinear_splitting():
    data = spectral_masses({'mS1Sq': 1210000, 'mS2Sq': 1000000,
                            'muT3GeneralS2Sq': {'real': 10000, 'imag': 0}})
    assert data['mS2_minus_sq_gev2'] == 980000
    assert data['mS2_plus_sq_gev2'] == 1020000
    for mu in (-10000, 500001, 10000+100j):
        with pytest.raises(ValueError):
            spectral_masses({'mS1Sq': 1210000, 'mS2Sq': 1000000,
                             'muT3GeneralS2Sq': mu})
