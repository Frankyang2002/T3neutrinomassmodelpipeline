import numpy as np
import pytest
from Numerical.diagnostics.T3BThresholdSubtractionAudit import calibrated_subtraction, alignment_failure_diagnostic


def block(a):
    a = np.asarray(a, dtype=complex)
    return {'real': a.real.tolist(), 'imag': a.imag.tolist()}


def fixture_payload():
    ll = np.diag([1., 2., 3.]).astype(complex) * 1e-8
    direct = 2 * ll
    hard = np.eye(3) * 4e-8
    fixed = hard + direct
    control = {'status': 'FixedOrderControlComplete', 'model':'T3-B alpha=-1',
               'mu_f_gev':1e5, 'mu_s_gev':1e3,'fixed_order_c5_at_ms':block(fixed)}
    bridge = {'status':'OffDiagonalFullFlavorDiagnosticComplete', 'model':'T3-B alpha=-1',
              'mu_f_gev':1e5, 'mu_s_gev':1e3, 'fixed_order_c5_at_ms':block(fixed),
              'c5_at_ms_frozen_mf_leading_log':block(ll),
              'c5_at_ms_forward_only_running_lambda':block(1.1 * ll),
              'c5_at_ms_coupled_offdiagonal_only':block(1.102 * ll)}
    return control, bridge, hard, direct


def test_subtraction_preserves_exact_fixed_order_at_equal_transfer():
    c, b, h, d = fixture_payload()
    b['c5_at_ms_forward_only_running_lambda'] = b['c5_at_ms_frozen_mf_leading_log']
    b['c5_at_ms_coupled_offdiagonal_only'] = b['c5_at_ms_frozen_mf_leading_log']
    r = calibrated_subtraction(c, b, h, d)
    assert r['status'] == 'OffDiagonalThresholdSubtractionDiagnosticComplete'
    assert r['relative_additional_offdiagonal_to_fixed_c5'] < 1e-14


def test_subtraction_detects_reverse_and_running_coupling():
    c, b, h, d = fixture_payload()
    r = calibrated_subtraction(c, b, h, d)
    assert abs(r['direct_vs_leading_log_complex_factor'][0] - 2) < 1e-10
    assert r['relative_additional_offdiagonal_to_fixed_c5'] > 0
    assert not r['neutrino_mass_prediction_available']


def test_fails_if_production_direct_misaligned():
    c, b, h, d = fixture_payload()
    d[0,1] = 1e-8
    d[1,0] = 1e-8
    c['fixed_order_c5_at_ms'] = block(h+d)
    b['fixed_order_c5_at_ms'] = block(h+d)
    with pytest.raises(ValueError, match='not proportional'):
        calibrated_subtraction(c,b,h,d)


def test_mismatch_diagnostic_never_produces_corrected_c5():
    c, b, h, d = fixture_payload()
    d[0,1] = 1e-8
    d[1,0] = 1e-8
    c['fixed_order_c5_at_ms'] = block(h+d)
    b['fixed_order_c5_at_ms'] = block(h+d)
    result = alignment_failure_diagnostic(c,b,h,d)
    assert result['status'] == 'AlignmentMismatch'
    assert result['direct_vs_leading_log_alignment_residual'] > 0
    assert result['subtraction_accepted'] is False
    assert 'comparison_proxy_c5_not_physical_prediction' not in result
