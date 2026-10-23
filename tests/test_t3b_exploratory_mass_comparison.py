import numpy as np
from Numerical.diagnostics.T3BExploratoryCoupledBenchmark import exploratory_mass_comparison


def test_identical_c5_gives_identical_masses():
    m = np.diag([1e-15, 2e-15, 4e-15]).astype(complex)
    result = exploratory_mass_comparison(m, m, np.diag([1e-5, 1e-3, 1e-2]))
    assert np.allclose(result['mass_differences_ev'], 0, atol=1e-11)
    assert result['validated_physical_prediction_available'] is False


def test_uniform_rescaling_rescales_all_masses():
    m = np.diag([1e-15, 2e-15, 4e-15]).astype(complex)
    result = exploratory_mass_comparison(m, 1.1*m, np.eye(3))
    assert np.allclose(result['relative_mass_differences'], 0.1, atol=1e-9)
    assert result['comparison_scale'] == 'MS (not low energy)'


def test_invalid_matrix_rejected():
    import pytest
    with pytest.raises(ValueError):
        exploratory_mass_comparison(np.zeros((2,2)), np.eye(3), np.eye(3))
