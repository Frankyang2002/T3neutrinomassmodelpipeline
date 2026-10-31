"""Tests for candidate B/C coupled UV -> EFT1 five-mass benchmark.

From repository root:
  python -X utf8 -m pytest -q tests/test_t3_bc_coupled_uv_eft1.py
"""
from pathlib import Path
from types import SimpleNamespace
import math

import numpy as np
import pytest

from scripts import RunT3BCCoupledUVToEFT1 as both
from scripts import RunT3BCFullScalarMassEFT1 as eft
from scripts import VerifyT3BCFullScalarMassClosure as closure


def test_beta_merging_keeps_every_uv_coupling_and_replaces_all_four_mass_betas():
    layout = ('gY', 'yu', 'y1', 'y2', 'MF', 'mS1Sq', 'mS2Sq',
              'lambdaH', 'lambdaT3GeneralQ110101I1')
    old = {x: x for x in layout}
    old[both.MIXED] = 'old_mixed'
    mass = dict(zip(closure.NAMES, ('new_higgs', 'new_s1', 'new_s2', 'new_mu')))
    merged = eft.merge_beta_reports(old, mass, layout)
    assert set(merged) == set(layout) | {'mHSq', both.MIXED}
    assert all(merged[k] == old[k] for k in ('gY', 'yu', 'y1', 'y2',
                                               'MF', 'lambdaH'))
    assert all(merged[k] == mass[k] for k in closure.NAMES)
    assert merged[both.MIXED] != 'old_mixed'


def test_majorana_projection_removes_only_roundoff():
    layout = SimpleNamespace(block=lambda name: SimpleNamespace(
        kind='complex_matrix', shape=(3, 3), start=0, stop=18))
    y = np.arange(23., dtype=float)
    y[:9] = np.diag([100., 110., 120.]).reshape(-1)
    y[9:18] = np.diag([1., 2., 3.]).reshape(-1)
    y[1] = 1.e-9
    copied = both._symmetrize_majorana(y, layout)
    assert copied[1] == pytest.approx(5.e-10)
    assert copied[3] == pytest.approx(5.e-10)
    np.testing.assert_array_equal(copied[18:], y[18:])
    assert y[1] == 1.e-9  # original solver vector not mutated
    broken = y.copy()
    broken[1] = 1.
    with pytest.raises(ValueError, match='antisymmetric'):
        both._symmetrize_majorana(broken, layout)


def test_uv_derivative_includes_heavy_couplings_and_loop_factor_once(monkeypatch):
    class Layout:
        size = 5
        names = ('MF', 'y1', 'mS1Sq', 'mS2Sq', 'gY')
        def block(self, name):
            # The symmetrisation is separately tested. Here use a scalar stub.
            raise AssertionError('not used in mocked symmetriser')

    class State:
        baseline = object()
        extra = {'new_quartic': 0.02 + 0.03j}

    monkeypatch.setattr(both, '_symmetrize_majorana', lambda y, layout: y)
    def unpack(y, layout, registry, *, uv, mu_gev):
        assert uv and mu_gev == 1.e7 and len(y) == 5
        return State()
    monkeypatch.setattr(both, 'unpack_general_scalar_state', unpack)
    monkeypatch.setattr(both, 'state_environment', lambda s: {
        'MF': np.eye(3)*1.e5, 'y1': np.eye(3)*0.02})

    def evaluate(expr, env):
        assert env['mHSq'] == -7800.
        assert env[both.MIXED] == 4200.+3100.j
        assert env['new_quartic'] == 0.02+0.03j
        assert 'MF' in env and 'y1' in env
        return {'mf':1., 'y':2., 's1':3., 's2':4., 'g':5.,
                'mixed':6.+7.j, 'h':8.}[expr]
    monkeypatch.setattr(both, 'evaluate_inputform_expression', evaluate)

    def pack(values, layout, *, divide_by_loop_factor):
        assert divide_by_loop_factor is True
        assert set(values) == set(layout.names)
        return np.array([values[name] for name in layout.names])/both.LOOP
    monkeypatch.setattr(both, 'beta_values_to_derivative', pack)
    betas = dict(zip(Layout.names, ('mf','y','s1','s2','g')))
    betas.update({both.MIXED:'mixed', 'mHSq':'h'})
    y = np.array([1., 2., 1.44e6, 1.96e6, 0.36, 4200., 3100., -7800.])
    dy = both.uv_derivative(y, mu_gev=1.e7, layout=Layout(),
                            registry={}, betas=betas)
    np.testing.assert_allclose(dy*both.LOOP, [1.,2.,3.,4.,5.,6.,7.,8.],
                               atol=1e-12, rtol=0)


def test_coupled_integrator_samples_both_endpoints():
    y0 = np.array([10., -7800.])
    rhs = lambda t, y: np.array([0., 2.])
    result = both._integrate(y0, rhs, 1e7, 1e5,
                            points=4, rtol=1e-8, atol=1e-11)
    assert result.success
    assert result.y.shape == (2, 4)
    np.testing.assert_allclose(result.y[0], np.repeat(10., 4))
    assert result.y[1,-1] == pytest.approx(-7800.+2.*math.log(1e5/1e7))


@pytest.mark.parametrize('model', both.MODELS)
def test_native_uv_mass_beta_affine_and_heavy_mass_dependence(model):
    tag = model.replace('-', 'm')
    mass_path = both.MASS_DIR / f'{tag}_uv_full_mass_v1.json'
    if not mass_path.is_file():
        pytest.skip('native full-mass UV beta JSON is not available locally')
    report = closure.check_payload(closure.read_payload(mass_path), model, 'UV')
    env = closure.benchmark_environment(model, 'UV')
    base = closure.evaluate_mass_vector(report, env)
    assert np.isfinite(base).all() and base.shape == (5,)
    # An MF->0 counterfactual removes the heavy-fermion inhomogeneous
    # beta contribution; it must change at least one BSM diagonal beta.
    without_mf = dict(env, MF=np.zeros((3, 3), dtype=complex))
    changed = closure.evaluate_mass_vector(report, without_mf)
    assert abs(changed[1]-base[1]) > 1e-6
    assert abs(changed[2]-base[2]) > 1e-6
    # Higgs portal feedback from the separately extracted *full* mass beta.
    shifted = closure.evaluate_mass_vector(
        report, dict(env, mHSq=env['mHSq']+2000.))
    np.testing.assert_allclose(shifted[1:3]-base[1:3],
        4*2000.*np.array([env['lambdaH1'],env['lambdaH2']]),
        rtol=1e-9, atol=1e-5)


@pytest.mark.parametrize('model', both.MODELS)
def test_existing_uv_beta_inventory_is_complete_at_fitted_start(model):
    stage29 = both.STAGE29_DIR / 'stage29_convention_uv_handoff.json'
    tag = model.replace('-', 'm')
    dimensionless = both.DIMENSIONLESS_DIR / f'{tag}_uv_quadratic_live_rgbeta.json'
    if not stage29.is_file() or not dimensionless.is_file():
        pytest.skip('Stage29 or native RGBeta UV report not available')
    registry = both.load_registry(both.legacy.REGISTRY)
    state, vector, layout, _ = both.load_uv_start(
        model, stage29_dir=both.STAGE29_DIR, registry=registry)
    report = both.checked_uv_dimensionless(
        model, registry=registry, layout=layout, path=dimensionless)
    assert 'MF' in report and 'y1' in report and 'y2' in report
    assert set(report) == set(layout.names) | {both.MIXED}
    assert vector.shape == (layout.size,)
    assert state.baseline.mu_gev == 1.e7
