import numpy as np
from scripts.RefitT3FrozenFullScalarLow import trial_seeds,whiten,uv_with_yukawas


def test_seed_scale_quartic():
    rng=np.random.default_rng(1)
    seeds=trial_seeds(np.ones(18),{'delta_m3l_sq_ev2':16.},
                      {'delta_m3l_sq_ev2':1.},3.,3,rng)
    assert np.allclose(seeds[1],.5)
    assert len(seeds)==3


def test_whiten_covariance():
    class Target:
        observable_names=('a','b')
        covariance=np.diag([4.,9.])
        central_values=np.array([1.,2.])
    assert np.allclose(whiten({'a':3.,'b':5.},Target()),[1.,1.])


def test_default_source_root_has_benchmark_segment():
    from pathlib import Path
    from scripts.RefitT3FrozenFullScalarLow import INPUT
    import inspect
    from scripts import RefitT3FrozenFullScalarLow as mod
    source=inspect.getsource(mod.main)
    assert "default=INPUT / 'optimised_low'" in source


def test_run_one_uses_source_root_model_subfolder():
    from scripts.RefitT3FrozenFullScalarLow import load_case
    import inspect
    source=inspect.getsource(load_case)
    assert 'folder=source_root/tag' in source
