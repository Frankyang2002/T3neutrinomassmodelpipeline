"""Minimal local checks for the new A(-4)/D(-2) source handoff."""
import numpy as np
import pytest

from Numerical.parallel_t3.complete.two_model_joint_scalar_batch9 import (
    IntegrationBlocked, real_scalar, quartic_source_path, MODELS,
)


def test_models_and_unsplit_symbolic_source_paths(tmp_path):
    assert MODELS == ('A_-4', 'D_-2')
    assert str(quartic_source_path(tmp_path, 'A_-4')).endswith(
        'output/full16_symbolic/T3_A_alpha_m4/data/eft1_after_F_scalar_quartic_seed.json')
    assert str(quartic_source_path(tmp_path, 'D_-2')).endswith(
        'output/full16_symbolic/T3_D_alpha_m2/data/eft1_after_F_scalar_quartic_seed.json')
    with pytest.raises(IntegrationBlocked):
        quartic_source_path(tmp_path, 'B_-1')


def test_real_scalar_never_silently_discards_complex_mass():
    assert real_scalar(np.asarray(1_000_000+0j), 'msq', positive=True) == 1_000_000
    with pytest.raises(IntegrationBlocked):
        real_scalar(1+0.5j, 'complex coupl')
    with pytest.raises(IntegrationBlocked):
        real_scalar(-2, 'tachyonic mass', positive=True)
    with pytest.raises(IntegrationBlocked):
        real_scalar(np.nan, 'invalid')
