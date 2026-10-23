import numpy as np
import pytest

from Numerical.diagnostics.T3BKernelFromControl import evaluate


def sample(lam=0.1 + 0.2j):
    return {
        'status': 'FixedOrderControlComplete',
        'model': 'T3-B alpha=-1',
        'intermediate_running': {'samples': [
            {'mu_gev': 1e5, 'lambdaT3': [lam.real, lam.imag]},
            {'mu_gev': 1e4, 'lambdaT3': [lam.real, lam.imag]},
            {'mu_gev': 1e3, 'lambdaT3': [lam.real, lam.imag]},
        ]},
        'fixed_order_c5_at_ms': {'frobenius_norm': 4.1e-7},
    }


def test_constant_quartic_kernel():
    z = 0.1 + 0.2j
    result = evaluate(sample(z))
    f = 4 * np.sqrt(3) / 3 * z / (16 * np.pi**2)
    b = 8 * np.sqrt(3) / 3 * z.conjugate() / (16 * np.pi**2)
    dt = np.log(1e3 / 1e5)
    assert complex(*result['integrated_s_to_kappa']) == pytest.approx(f * dt)
    assert complex(*result['integrated_kappa_to_s']) == pytest.approx(b * dt)
    assert complex(*result['second_order_s_to_kappa_to_s']) == pytest.approx(b * f * dt**2 / 2)
    assert complex(*result['second_order_kappa_to_s_to_kappa']) == pytest.approx(f * b * dt**2 / 2)


def test_zero_quartic():
    r = evaluate(sample(0j))
    assert r['integrated_s_to_kappa'] == [0.0, 0.0]
    assert r['second_order_s_to_kappa_to_s'] == [0.0, 0.0]


def test_rejects_bad_scale_order():
    c = sample()
    c['intermediate_running']['samples'].reverse()
    with pytest.raises(ValueError, match='ordered'):
        evaluate(c)
