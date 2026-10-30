"""Unit-level guards for all-16 zero-extra scalar controls."""
from dataclasses import dataclass

import pytest

from scripts import RunT3AllZeroExtraScalarControl as ctrl


def test_model_and_mode_selection():
    assert ctrl.select_csv('A_+0,B_-1,E_-4', ctrl.policy.MODELS) == (
        'A_+0','B_-1','E_-4')
    assert ctrl.select_csv('held_zero,boundary_zero', ctrl.MODES) == ctrl.MODES
    for bad in ('', 'B_-1,', 'B_-1,B_-1', 'Z_+0'):
        with pytest.raises(ValueError):
            ctrl.select_csv(bad, ctrl.policy.MODELS)


def test_extra_beta_projection_only_changes_additional_coordinates():
    betas = {'lambdaT3':'normal topology beta', 'lambdaH':'normal Higgs beta',
             'newQ':'extra beta', 'reusedQ':'extra beta', 'mu12':'mass beta'}
    reduced = ctrl.zero_extra_betas(betas, ('newQ',), ('reusedQ',), 'mu12')
    assert reduced == {'lambdaT3':'normal topology beta',
        'lambdaH':'normal Higgs beta', 'newQ':'0', 'reusedQ':'0', 'mu12':'0'}
    assert betas['newQ'] == 'extra beta'
    with pytest.raises(ValueError):
        ctrl.zero_extra_betas(betas, ('lambdaT3',), (), None)
    with pytest.raises(ValueError):
        ctrl.zero_extra_betas(betas, ('missing',), (), None)


def test_zero_mass_for_all_modes():
    for mode in ctrl.MODES:
        settings = ctrl.parameters('optimised_low','A_+0',mode,16)
        assert settings['mu2_fraction'] == 0
        assert settings['frozen_uv'] is True
        assert settings['historical_yukawas_reoptimised_for_full_scalar'] is False


def test_zero_additional_quartics_preserves_original_parameters(monkeypatch):
    @dataclass(frozen=True)
    class Base:
        lambdaT3: complex = 5e-9+0j
        lambdaH1: float = .06
        lambdaHHdagS2S2: complex = -.01+0j
        mS1Sq: float = 1.25e6
        mS2Sq: float = 1.0e6
        def validated(self):
            return self

    @dataclass(frozen=True)
    class Scalar:
        baseline: Base
        extra: dict
        def validated(self, registry):
            return self

    monkeypatch.setattr(ctrl.historical,'inventory',
         lambda registry,model: ({},('newQ',),('lambdaHHdagS2S2',)))
    monkeypatch.setattr(ctrl.full,'GeneralT3ScalarRGState',Scalar)
    original = Scalar(Base(),{'newQ':-.013+0j})
    result, new, reused = ctrl.restrict_scalar(
        {'registry':{},'scalar':original}, 'B_-1')
    assert (new,reused) == (('newQ',),('lambdaHHdagS2S2',))
    assert result.baseline.lambdaT3 == original.baseline.lambdaT3
    assert result.baseline.lambdaH1 == original.baseline.lambdaH1
    assert result.baseline.mS1Sq == original.baseline.mS1Sq
    assert result.baseline.mS2Sq == original.baseline.mS2Sq
    assert result.baseline.lambdaHHdagS2S2 == 0j
    assert result.extra['newQ'] == 0j
    assert original.baseline.lambdaHHdagS2S2 == -.01+0j
    assert original.extra['newQ'] == -.013+0j


def test_extra_endpoint_for_new_and_reused():
    from types import SimpleNamespace
    state = SimpleNamespace(extra={'newQ':0.002j},
        baseline=SimpleNamespace(reusedQ=0.0+0j))
    out = ctrl.extra_endpoint_values(state, ('newQ',), ('reusedQ',))
    assert out == {'newQ':0.002, 'reusedQ':0.0}
