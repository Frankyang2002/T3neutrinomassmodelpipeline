"""Stage 16 scope and threshold-ordering regression tests (no production solver mocked as proof)."""
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from Numerical.experimental.GeneralT3CoupledUV import construct_uv_state, ALL_COMPLEX
from Numerical.experimental.GeneralT3HierarchicalTrajectory import (
    _threshold_scales, _guard_running_hierarchy, _validate_payload,
    _transport_summary, _matrix_json, _complex_json, ThresholdScales,
)


def bench():
    return json.loads((ROOT/'docs/general_t3_threshold_benchmark.json').read_text())

def initial_proxy(config):
    return SimpleNamespace(mu_gev=float(config['scale_uv_gev']),
        MF=np.diag(np.asarray(config['MF'],dtype=complex)),
        mS1Sq=float(config['MS1'])**2, mS2Sq=float(config['MS2'])**2)


@pytest.mark.parametrize('model', ['B_m1','C_m1'])
def test_good_hierarchy_and_fields(model):
    config=bench()
    t=_threshold_scales(config,initial_proxy(config))
    assert t.mu_f > t.mu_s > t.mu_low
    assert min(config['MF'])>t.mu_f>max(config['MS1'],config['MS2'])
    if (ROOT/'Numerical/core/State.py').exists():
        state=construct_uv_state(model,config)
        assert state.representation.alpha==-1
        assert all(getattr(state,name)==0 for name in ALL_COMPLEX)


def test_previous_nonhierarchical_benchmark_rejected():
    data=json.loads((ROOT/'docs/general_c5_toy_benchmark.json').read_text())
    data['thresholds']={'mu_F':1900.,'mu_S':1100.,'mu_low':246.22}
    with pytest.raises(ValueError,match='fermion-first mass hierarchy'):
        _threshold_scales(data,initial_proxy(data))


@pytest.mark.parametrize('bad', [
    {'mu_F':1200,'mu_S':1100,'mu_low':246.22},
    {'mu_F':4500,'mu_S':1900,'mu_low':246.22},
    {'mu_F':4400,'mu_S':1100,'mu_low':2000},
    {'mu_F':float('nan'),'mu_S':1100,'mu_low':246.22},
])
def test_invalid_hierarchy_rejected(bad):
    data=bench();data['thresholds']=bad
    with pytest.raises(ValueError):
        _threshold_scales(data,initial_proxy(data))


def test_eft1_requires_all_six_beta_functions():
    payload={'status':'Success', 'metadata':{'ReportBetaConvention':'16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g',
        'dS1':2,'dS2':2,'dF':1,'alpha':-1, 'IntegratedField':'F', 'ActiveBSMFields':['S1','S2']},
             'report_betas':{name:'0' for name in ALL_COMPLEX}}
    _validate_payload(payload,'EFT1','B_m1')
    del payload['report_betas']['lambdaHHdagS2S2']
    with pytest.raises(ValueError,match='missing'):
        _validate_payload(payload,'EFT1','B_m1')


def test_wrong_beta_convention_rejected():
    payload={'status':'Success','metadata':{'ReportBetaConvention':'wrong'},
             'report_betas':{name:'0' for name in ALL_COMPLEX}}
    with pytest.raises(ValueError,match='convention'):
        _validate_payload(payload,'UV','B_m1')


def test_running_threshold_mass_order_guard():
    t=ThresholdScales(4400.,1100.,246.22)
    trajectory=SimpleNamespace(fermion_threshold=SimpleNamespace(masses_gev=np.array([5000.,6000.,7000.])),
                               scalar_threshold=SimpleNamespace(masses_gev=(1350.,1780.)))
    assert _guard_running_hierarchy(trajectory,t)['heavy_at_mu_f_gev'][0]==5000.
    trajectory.scalar_threshold.masses_gev=(1350.,4600.)
    with pytest.raises(ValueError,match='Running masses violate'):
        _guard_running_hierarchy(trajectory,t)


def test_no_wilson_insertion_even_when_baseline_available():
    a=bench(); t=ThresholdScales(4400.,1100.,246.22)
    def state(scale,val):
        return SimpleNamespace(mu_gev=scale,**{name:complex(val) for name in ALL_COMPLEX})
    uv=state(20000.,0)
    f=state(4400.,-1e-12)
    eft=state(4400.,-1e-12)
    s=state(1100.,-1.2e-12)
    fake=SimpleNamespace(uv_initial_state=uv,uv_threshold_state=f,
        intermediate_initial_state=eft,eft1_threshold_state=s,
        fermion_threshold=SimpleNamespace(masses_gev=[5000.,6000.,7000.]),
        scalar_threshold=SimpleNamespace(masses_gev=[1350.,1780.]),
        uv=SimpleNamespace(nfev=23),eft1=SimpleNamespace(nfev=31))
    k=np.eye(3,dtype=complex)*1e-15
    r=_transport_summary(fake,k,k,'B_m1',a,t)
    assert r['baseline_hierarchical_c5_included']
    assert r['extra_c5_at_mu_s_gev_inverse'] is None
    assert r['extra_neutrino_mass_low_ev'] is None
    assert r['extra_c5_status']=='BlockedMissingExtraEFT1WilsonMatchingAndAnomalousDimensions'
    assert r['quartics_rgbeta']['lambdaHHdagS2S2']['eft1_at_mu_s']['real']==-1.2e-12
    assert not r['strict_per_mass_thresholds']


def test_f_boundary_discontinuity_rejected():
    a=bench();t=ThresholdScales(4400.,1100.,246.22)
    def state(scale,value):
        return SimpleNamespace(mu_gev=scale,**{name:complex(value) for name in ALL_COMPLEX})
    fake=SimpleNamespace(uv_initial_state=state(20000.,0),
        uv_threshold_state=state(4400.,1e-12),
        intermediate_initial_state=state(4400.,2e-12),
        eft1_threshold_state=state(1100.,2e-12),
        fermion_threshold=SimpleNamespace(masses_gev=[5000.,6000.,7000.]),
        scalar_threshold=SimpleNamespace(masses_gev=[1350.,1780.]),
        uv=SimpleNamespace(nfev=23),eft1=SimpleNamespace(nfev=31))
    with pytest.raises(ValueError,match='unaccounted discontinuity'):
        _transport_summary(fake,None,None,'C_m1',a,t)


def test_nonfinite_matrix_rejected():
    with pytest.raises(ValueError):
        _matrix_json(np.eye(3)*float('nan'))
    with pytest.raises(ValueError):
        _complex_json(complex(float('inf'),0))


def test_production_scalar_eft_beta_files_if_available():
    for m in 'BC':
        f=ROOT/f'output/T3_{m}_alpha_m1/data/eft1_rgbeta_rge.json'
        if not f.exists():
            pytest.skip('Only user local repository contains EFT1 beta exports')
        _validate_payload(json.loads(f.read_text(encoding='utf-8-sig')),'EFT1',f'{m}_m1')


def test_wrong_model_and_wrong_active_fields_fail():
    payload={'status':'Success','metadata':{'ReportBetaConvention':'16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g',
            'dS1':2,'dS2':2,'dF':1,'alpha':-1,'IntegratedField':'F','ActiveBSMFields':['S1','S2']},
             'report_betas':{name:'0' for name in ALL_COMPLEX}}
    with pytest.raises(ValueError,match='metadata'):
        _validate_payload(payload,'EFT1','C_m1')
    payload['metadata']['ActiveBSMFields']=['S1','F']
    with pytest.raises(ValueError,match='only S1/S2'):
        _validate_payload(payload,'EFT1','B_m1')
