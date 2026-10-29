"""Fast regression checks for the candidate nonzero quadratic spectrum audit."""
from __future__ import annotations

import pytest

from scripts.RunT3QuadraticMixingValidation import run, EXPECTED_KINDS


def test_full_16_model_nonzero_quadratic_inventory():
    report = run()
    assert report['models_checked'] == 16
    assert report['extra_quadratic_model_slots_checked_nonzero'] == 8
    assert report['mixed_doublet_tensor_checks'] == 2
    assert report['global_matchete_rgbeta_normalisation_certified'] is False
    assert report['physical_scalar_threshold_matching_performed'] is False
    assert report['full_general_one_loop_c5_certified'] is False
    for model, entry in report['models'].items():
        assert entry['quadratic_kind'] == EXPECTED_KINDS.get(model, 'none')
        assert entry['analytic_hessian_agreement'] is True
        assert entry['phase_conjugation_sign_invariance'] is True
        assert entry['zero_boundary_limit_verified'] is True
        assert entry['minimum_mass_squared_gev2'] > 0


def test_phase_and_degenerate_diagonal_mass_inputs():
    report = run(m_s1_gev=1500., m_s2_gev=1500.,
                 mu2_real_gev2=13000., mu2_imag_gev2=-17000.)
    assert report['models_checked'] == 16
    for model in ('B_-1', 'C_-1'):
        assert report['models'][model]['sin_2theta'] == pytest.approx(1., abs=1e-12)


@pytest.mark.parametrize('invalid', [0., -1., float('nan'), float('inf')])
def test_invalid_scalar_mass_rejected(invalid):
    with pytest.raises(ValueError):
        run(m_s1_gev=invalid)


def test_zero_extra_coefficient_rejected_as_nonzero_validation():
    with pytest.raises(ValueError):
        run(mu2_real_gev2=0., mu2_imag_gev2=0.)
