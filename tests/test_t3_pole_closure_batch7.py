"""Isolated Batch7 regression tests; no Windows benchmark files required."""
from __future__ import annotations

import numpy as np
import pytest

from Numerical.parallel_t3.wilson.native_pole_closure_batch7 import (
    MODELS, COMPLETE_MATCHING_GATES, matrix_record, relative,
    symmetry_split, physical_total, source_by_native_label,
    source_cancellation_ledger, linear_coefficient_diagnostic,
    native_pole_evidence, scalar_integral_checks, certify_matching,
)


def sym(x):
    return (x + x.T)/2


def sample_total():
    z = np.array([[1.+1j, 2-1j, 0], [2-1j, 3, .5j], [0, .5j, 4]], complex)*1e-8
    return z


def sample_anti():
    return np.array([[0, 2j, 1], [-2j, 0, 3j], [-1, -3j, 0]],complex)*1e-8


def test_catalogue_all_sixteen():
    assert len(MODELS) == 16
    assert len(set(MODELS)) == 16
    assert 'B_-1' in MODELS


def test_complex_majorana_transpose_not_hermitian():
    s = sample_total()
    assert not np.allclose(s, s.conj().T, atol=1e-12)
    assert relative(physical_total(s, 'Majorana'),s) < 1e-15


def test_unsymmetric_source_is_serializable_unmodified():
    x = sample_total() + sample_anti()
    r = matrix_record(x, 'quartic source')
    y = np.asarray(r['real']) + 1j*np.asarray(r['imag'])
    assert np.array_equal(x, y)
    assert r['transpose_residual'] > .1
    assert r['physical_total_validation_applied'] is False


def test_physical_total_rejects_unsymmetric_matrix():
    with pytest.raises(ValueError, match='physical C5'):
        physical_total(sample_total()+sample_anti(), 'combined')
    with pytest.raises(ValueError):
        matrix_record(sample_total()+sample_anti(), 'physical', physical=True)


def test_matrix_split_preserves_full_complex_tensor():
    x = sample_total()+sample_anti()
    p = symmetry_split(x)
    assert relative(p['symmetric']+p['antisymmetric'], x) < 1e-14
    assert np.array_equal(p['antisymmetric'].T, -p['antisymmetric'])
    assert p['antisymmetric_fraction'] > 0


def test_source_antisymmetric_terms_cancel_exactly():
    s, a = sample_total(), sample_anti()
    raw = np.stack([.25*s+a, .75*s-a])
    ledger = source_cancellation_ledger(raw,['lambdaA','lambdaB'],s,label='EFT1')
    assert ledger['source_count'] == 2
    assert ledger['reassembly_residual'] < 1e-12
    assert ledger['sum_abs_antisymmetric_norm'] > 0
    assert abs(ledger['antisymmetric_cancellation_fraction']-1) < 1e-12
    assert ledger['per_source']['lambdaA']['transpose_residual']>0.1


def test_source_aggregate_symmetry_failure_blocks():
    s, a = sample_total(), sample_anti()
    with pytest.raises(ValueError):
        source_cancellation_ledger([s,a],['a','b'],s,label='bad')


def test_source_incomplete_sum_blocks():
    s = sample_total()
    with pytest.raises(ValueError, match='reassembly'):
        source_cancellation_ledger([s/2],['half'],s,label='bad')


def test_native_grouping_retains_distinct_labels_and_raw_phase():
    s, a = sample_total(), sample_anti()
    arr = np.zeros((3,3,3,3),complex)
    arr[0,0] = .25*s+a
    arr[1,0] = .75*s-a
    arr[2,2] = .1*s
    group = source_by_native_label(arr,['lambdaOne','lambdaOne','lambdaTwo'])
    assert list(group) == ['lambdaOne','lambdaTwo']
    assert relative(group['lambdaOne'],s) < 1e-14
    assert relative(group['lambdaTwo'],.1*s) < 1e-14


def test_native_group_rejects_missing_label_and_nonfinite():
    arr = np.zeros((1,3,3,3),complex)
    with pytest.raises(ValueError, match='label'):
        source_by_native_label(arr,['not_a_quartic'])
    arr[0,0,1,2] = np.nan
    with pytest.raises(ValueError, match='Nonfinite'):
        source_by_native_label(arr,['lambdaA'])


def test_scalar_pole_and_B0_are_independently_verified():
    s = sample_total()
    checks = scalar_integral_checks(explicit_ct=s, batch5_soft_pole=-s,
                                     b0_scale_ct=s)
    assert checks['EFT1_scalar_bubble_self_pole_cancelled']
    assert not checks['complete_UV_pole_cancelled']
    with pytest.raises(ValueError, match='disagrees'):
        scalar_integral_checks(explicit_ct=s, batch5_soft_pole=-s,
                               b0_scale_ct=1.2*s)


def test_complex_normalisation_fit_is_diagnostic_not_matching():
    s = sample_total()
    f = linear_coefficient_diagnostic(s,(2+3j)*s)
    assert f['fit_defined']
    assert f['postfit_residual_relative'] < 1e-14
    assert f['coefficient_real'] == pytest.approx(2)
    assert f['coefficient_imag'] == pytest.approx(3)
    assert not f['fit_used_in_physics']


def test_zero_reference_fit_is_undefined_not_silent_zero():
    r = linear_coefficient_diagnostic(np.zeros((3,3)),sample_total())
    assert not r['fit_defined']
    assert r['coefficient_real'] is None


def test_native_versus_ct_never_certifies_uv_counterterms():
    s = sample_total()
    d = native_pole_evidence(native=2*s, eft1_ct=s,
                             candidate_adm=s, candidate_scalar_adm=s)
    assert d['native_vs_scalar_CT_relative'] == pytest.approx(1)
    assert not d['native_is_proven_bare_UV_amplitude']
    assert not d['UV_pole_cancellation_verified']
    assert not d['counterterm_fitted_to_cancel_native']
    assert d['native_minus_EFT1_scalar_CT']['frobenius_GeV_inverse'] > 0


def test_certification_fails_closed():
    cert = certify_matching()
    assert len(cert['missing_gates']) == len(COMPLETE_MATCHING_GATES)
    assert not cert['complete_matching_claim_allowed']
    with pytest.raises(ValueError, match='Unknown'):
        certify_matching({'fictional_confirmation': True})


def test_certification_requires_all_explicit_verified_inputs():
    explicit = {k: True for k in COMPLETE_MATCHING_GATES}
    explicit['signed_GammaCC_Weyl_normalisation_verified'] = False
    result = certify_matching(explicit)
    assert not result['complete_matching_claim_allowed']
    assert result['missing_gates'] == ['signed_GammaCC_Weyl_normalisation_verified']


def test_finite_tensor_guard_and_shape_validation():
    with pytest.raises(ValueError):
        matrix_record(np.zeros((2,2)), 'wrong')
    with pytest.raises(ValueError):
        matrix_record(np.full((3,3), np.inf), 'inf')
    with pytest.raises(ValueError):
        relative(np.zeros((3,3)),np.zeros((2,2)))
