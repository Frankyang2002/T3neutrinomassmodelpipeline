"""Small deterministic contracts for continuous all-scalar UV/EFT1 machinery.

These tests do not pretend to certify cross-package physical matching.
"""
from __future__ import annotations

import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from scripts import RunT3FullScalarContinuous as run


def _fixture(kind='mixed'):
    model = 'B_-1'
    symbol = run.MODEL_MASS[kind]
    registry = {'models': {model: {
        'total_additional_quartics': 1,
        'preexisting_rgbeta_symbols': [],
        'new_rgbeta_symbols': [{'name': 'newQ'}],
    }}}
    meta = {
        'GeneralScalarModel': model, 'GeneralScalarStage': 'UV',
        'GeneralScalarBasisStatus': 'CandidateNativeRGBetaUnmappedToMatchete',
        'GeneralScalarQuarticBetaBasis': 'RGBetaFullQuarticProjection',
        'GeneralScalarQuadraticMassMixing': 'IncludedCandidate',
        'GeneralScalarAdditionalQuartics': 1,
        'GeneralScalarExistingStructuralCount': 0,
        'QuadraticKind': kind, 'QuadraticMassSymbol': symbol,
        'QuadraticPrefactorConvention': 'V2 = mu2 * invariant + h.c.; no 1/2',
        'PhysicalMatcheteRGBetaMassBasisCertified': False,
        'ReportBetaConvention': run.EXPECTED_BETA_CONVENTION,
        'SharedScalar': False, 'dS1': 2, 'dS2': 2, 'dF': 1, 'alpha': -1,
        'GeneralScalarNewBetaKeys': ['newQ'],
    }
    betas = {'g2': '1', 'newQ': '0'}
    if symbol:
        betas[symbol] = '0'
    payload = {
        'format': 'general-t3-quadratic-rgbeta-experimental-v1',
        'status': 'SuccessCandidateMassBeta',
        'physical_mapping_certified': False,
        'mass_eigenstate_threshold_matching_performed': False,
        'metadata': meta, 'report_betas': betas,
        'beta_symbols': list(betas),
    }
    rep = SimpleNamespace(d_s1=2, d_s2=2, d_f=1, alpha=-1)
    layout = SimpleNamespace(names=('g2','newQ'), size=3)
    return payload, registry, rep, layout


def test_beta_inventory_complete_quadratic_uv():
    data, registry, rep, layout = _fixture()
    assert set(run.beta_inventory(data, model='B_-1', stage='UV',
          representation=rep, layout=layout, registry=registry,
          kind='mixed')) == {'g2','newQ','muT3GeneralS1S2'}


def test_beta_inventory_rejects_missing_quadratic_beta():
    data, registry, rep, layout = _fixture()
    data['report_betas'].pop('muT3GeneralS1S2')
    data['beta_symbols'].remove('muT3GeneralS1S2')
    with pytest.raises(ValueError, match='beta inventory'):
        run.beta_inventory(data, model='B_-1', stage='UV',
                           representation=rep, layout=layout, registry=registry,
                           kind='mixed')


def test_beta_inventory_rejects_wrong_stage_or_certification():
    data, registry, rep, layout = _fixture()
    data['metadata']['GeneralScalarStage'] = 'EFT1'
    with pytest.raises(ValueError, match='inconsistent GeneralScalarStage'):
        run.beta_inventory(data, model='B_-1', stage='UV',
                           representation=rep, layout=layout, registry=registry,
                           kind='mixed')
    data['metadata']['GeneralScalarStage'] = 'UV'
    data['physical_mapping_certified'] = True
    with pytest.raises(ValueError, match='certification'):
        run.beta_inventory(data, model='B_-1', stage='UV',
                           representation=rep, layout=layout, registry=registry,
                           kind='mixed')


def test_beta_inventory_none_quadratic_does_not_invent_mass_beta():
    data, registry, rep, layout = _fixture(kind='none')
    assert set(run.beta_inventory(data, model='B_-1', stage='UV',
          representation=rep, layout=layout, registry=registry,
          kind='none')) == {'g2', 'newQ'}


def test_continuous_integrator_includes_extra_coordinate(monkeypatch):
    """Verify no accidental omission of the last (complex mu2) state block."""
    def synthetic_rhs(t, y, **kwargs):
        assert y.shape == (4,)
        return np.array([0., 0., -y[2]/run.LOOP, -y[3]/run.LOOP])
    monkeypatch.setattr(run, 'coupled_rhs', synthetic_rhs)
    y0=np.array([1.,2.,3.,4.])
    sol=run.integrate_segment(y0, mu_start=1e7, mu_end=1e5,
        layout=SimpleNamespace(size=2), registry={}, stage='UV',
        kind='mixed', betas={}, points=8)
    t = math.log(1e5/1e7)
    assert sol.y.shape == (4,8)
    np.testing.assert_allclose(sol.y[:2,-1], y0[:2], rtol=0, atol=1e-11)
    np.testing.assert_allclose(sol.y[2:,-1], y0[2:]*math.exp(-t/run.LOOP), rtol=1e-8)


def test_resume_requires_identical_source_hashes(tmp_path):
    p = tmp_path/'saved.json'
    p.write_text('{"format":"'+run.FORMAT+'","status":"'+run.STATUS+'",'
                 '"input_sha256":{"beta":"123"},"settings":{"mu":1}}',encoding='utf-8')
    assert run.result_is_resumable(p, expected_inputs={'beta':'123'},
                                   expected_settings={'mu':1}) is True
    with pytest.raises(ValueError, match='changed'):
        run.result_is_resumable(p, expected_inputs={'beta':'new'},
                                expected_settings={'mu':1})


def test_run_id_and_model_path_are_strict():
    assert run.model_tag('B_-1') == 'B_m1'
    assert run.run_id('run_20261029') == 'run_20261029'
    with pytest.raises(ValueError):
        run.run_id('../unsafe')
    with pytest.raises(ValueError):
        run.model_tag('not_a_model')

@pytest.mark.parametrize('stage', ['UV', 'EFT1'])
def test_quadratic_beta_backreacts_on_mass_and_quartic(monkeypatch, stage):
    """A nonzero holomorphic mass must affect *other* beta coordinates as well."""
    class Base:
        mS1Sq = 10.
    def unpack(values, layout, registry, **kwargs):
        return SimpleNamespace(baseline=Base(), extra={'newQ':complex(values[1],values[2])})
    def env(base):
        return {'mS1Sq':base.mS1Sq}
    def evaluator(expr, env):
        mu=env['muT3GeneralS1S2']
        if expr=='mass': return 2*mu.real
        if expr=='quartic': return (1.+2.j)*mu
        if expr=='mix': return env['newQ']*env['mS1Sq']
        raise AssertionError(expr)
    def pack(values, layout, **kwargs):
        z=values['newQ']; factor=run.LOOP if kwargs['divide_by_loop_factor'] else 1.
        return np.array([float(values['mS1Sq']),complex(z).real,complex(z).imag])/factor
    monkeypatch.setattr(run, 'unpack_general_scalar_state', unpack)
    monkeypatch.setattr(run, 'state_environment', env)
    monkeypatch.setattr(run, 'intermediate_scalar_state_environment', env)
    monkeypatch.setattr(run, 'evaluate_inputform_expression', evaluator)
    monkeypatch.setattr(run, 'beta_values_to_derivative', pack)
    y=np.array([10.,0.2,0.3,0.1,0.4])
    ret=run.coupled_rhs(1.,y,layout=SimpleNamespace(size=3),registry={},
                        stage=stage,kind='mixed',betas={
                            'mS1Sq':'mass','newQ':'quartic',
                            'muT3GeneralS1S2':'mix'})
    assert ret.shape == (5,)
    np.testing.assert_allclose(ret[0], 2*.1/run.LOOP)
    np.testing.assert_allclose(ret[1:3], np.array([-.7, .6])/run.LOOP)
    np.testing.assert_allclose(ret[3:], np.array([2., 3.])/run.LOOP)


def test_atomic_checkpoint_retains_previous_json_on_failed_serialization(tmp_path):
    path=tmp_path/'checkpoint.json'
    run.write_json_atomic(path, {'status':'previous'})
    with pytest.raises(TypeError):
        run.write_json_atomic(path, {'non_jsonable': object()})
    assert run.read(path)=={'status':'previous'}
