"""Unit tests for isolated quadratic-mass prototype, not production RGE."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from Numerical.experimental.ScalarQuadraticPrototype import (
    mixed_t3_bc_doublet_spectrum,
    neutral_complex_scalar_splitting,
)


def test_no_doublet_mixing_recovers_input_masses():
    s = mixed_t3_bc_doublet_spectrum(700**2, 900**2, 0j)
    assert s.mass == pytest.approx((700, 900))
    assert np.allclose(s.unitary.conj().T @ s.mass_matrix @ s.unitary,
                       np.diag(s.mass_squared))


def test_complex_doublet_mixing_matches_analytic_formula():
    a, b, mu = 700.0**2, 900.0**2, complex(200.0**2, 150.0**2)
    s = mixed_t3_bc_doublet_spectrum(a, b, mu)
    analytic = (
        (a + b - np.sqrt((a - b)**2 + 4 * abs(mu)**2)) / 2,
        (a + b + np.sqrt((a - b)**2 + 4 * abs(mu)**2)) / 2,
    )
    assert s.mass_squared == pytest.approx(analytic)
    assert np.allclose(s.unitary.conj().T @ s.mass_matrix @ s.unitary,
                       np.diag(s.mass_squared), atol=1e-8)
    assert np.allclose(s.unitary.conj().T @ s.unitary, np.eye(2), atol=1e-12)


def test_doublet_tachyon_is_rejected():
    with pytest.raises(ValueError, match='positive'):
        mixed_t3_bc_doublet_spectrum(100, 100, 200)


def test_complex_real_representation_splitting_and_degeneracies():
    for dim in (1, 3):
        s = neutral_complex_scalar_splitting(750**2, complex(100**2, 7500),
                                             su2_dimension=dim)
        assert s.mass_squared == pytest.approx((750**2 - 12500, 750**2 + 12500))
        assert s.real_multiplet_degeneracy == (dim, dim)


def test_invalid_complex_scalar_representation():
    with pytest.raises(ValueError, match='singlet/triplet'):
        neutral_complex_scalar_splitting(1000, 50, su2_dimension=2)
