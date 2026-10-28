"""Batch20 fail-closed offline tests; no native Wolfram/real RGBeta invoked."""
import json
import math
from pathlib import Path
import sys
import types

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts import RunT3BPlus1GeneralEFT1 as b20


def valid_registry():
    r = {'models': {f'm{i}':{} for i in range(15)}}
    r['models'][b20.MODEL] = {'dimensions':[2,2,1],'alpha':1,
        'total_additional_quartics':1,'preexisting_rgbeta_symbols':[],
        'new_rgbeta_symbols':[{'name':b20.NEW}]}
    return r


def valid_payload():
    names = ['gY','g2','g3','yu','yd','ye','lambdaH',
             'mS1Sq','mS2Sq','lambdaS1','lambdaS2','lambdaH1','lambdaH2',
             'lambda12','lambdaT3','lambdaH1Adj','lambdaH2Adj','lambda12Adj']
    b = {x:'1' for x in names}
    b['lambdaH1'] = b20.NEW+'*Bar['+b20.NEW+']'
    b[b20.NEW] = '-2*'+b20.NEW
    return {'format':b20.FORMAT,'status':'Success','model':b20.MODEL,
        'metadata':{'SharedScalar':False,'dS1':2,'dS2':2,'dF':1,'alpha':1,
                    'IntegratedField':'F','ActiveBSMFields':['S1','S2'],
                    'ReportBetaConvention':b20.CONVENTION},
        'physical_coupling_mapping_certified':False,
        'wilson_matching_performed':False,
        'stage29_quartic_evolved_in_this_export':False,
        'base_beta_symbols':names,'new_quartic_symbols':[b20.NEW],
        'report_betas':b}


@pytest.fixture
def fake_packer(monkeypatch):
    layout = types.SimpleNamespace(size=18,names=tuple(valid_payload()['base_beta_symbols']))
    def pack(s):return np.ones(18), layout
    # Import in validated function requires package chain as in actual repo.
    numerical=types.ModuleType('Numerical');running=types.ModuleType('Numerical.running')
    vector=types.ModuleType('Numerical.running.IntermediateScalarStateVector')
    vector.pack_intermediate_scalar_state=pack
    for name, module in [('Numerical',numerical),('Numerical.running',running),
                         ('Numerical.running.IntermediateScalarStateVector',vector)]:
        monkeypatch.setitem(sys.modules,name,module)
    state=types.SimpleNamespace(representation=types.SimpleNamespace(
        d_s1=2,d_s2=2,d_f=1,alpha=1,shared_scalar=False))
    return state,layout


@pytest.mark.parametrize('value', ['a', 'x_2', 'b-3', 'RUN2026'])
def test_run_id_accept(value):assert b20.validate_run_id(value)==value

@pytest.mark.parametrize('value', ['', '../a', 'hello world', 'x/y', 'x'*65,'é', None])
def test_run_id_reject(value):
    with pytest.raises(ValueError):b20.validate_run_id(value)


def test_registry_correct():b20.validate_registry(valid_registry())

@pytest.mark.parametrize('change', [
    ('alpha',-1),('dimensions',[2,3,1]),('total_additional_quartics',2),
    ('preexisting_rgbeta_symbols',[{'existing_name':'bad'}]),
    ('new_rgbeta_symbols',[]),
])
def test_registry_reject(change):
    r=valid_registry();r['models'][b20.MODEL][change[0]]=change[1]
    with pytest.raises(ValueError):b20.validate_registry(r)


def test_registry_requires_16():
    r=valid_registry();r['models'].pop('m0')
    with pytest.raises(ValueError):b20.validate_registry(r)


def test_sidecar_complex():
    assert b20.sidecar_complex({'real':-0.02,'imag':0.005})==complex(-0.02,0.005)

@pytest.mark.parametrize('raw',[{}, {'real':0}, {'real':0,'imag':float('nan')},
                                 {'real':0,'imag':0,'wrong':0}, None])
def test_bad_sidecar(raw):
    with pytest.raises(ValueError):b20.sidecar_complex(raw)


def test_validate_export_passes(fake_packer):
    state,layout=fake_packer
    betas,actual=b20.validate_export(valid_payload(),state)
    assert actual is layout and len(betas)==19 and b20.NEW in betas

@pytest.mark.parametrize('field,value',[('model','C_-1'),('physical_coupling_mapping_certified',True),
    ('wilson_matching_performed',True),('stage29_quartic_evolved_in_this_export',True),
    ('new_quartic_symbols',[]),('base_beta_symbols',['gY'])])
def test_export_metadata_fails(fake_packer,field,value):
    p=valid_payload();p[field]=value
    with pytest.raises(ValueError):b20.validate_export(p,fake_packer[0])


def test_export_wrong_metadata(fake_packer):
    p=valid_payload();p['metadata']['alpha']=-1
    with pytest.raises(ValueError):b20.validate_export(p,fake_packer[0])


def test_export_missing_extra_rhs(fake_packer):
    p=valid_payload();p['report_betas'].pop(b20.NEW)
    with pytest.raises(ValueError):b20.validate_export(p,fake_packer[0])


def test_export_no_feedback(fake_packer):
    p=valid_payload();p['report_betas']['lambdaH1']='0'
    with pytest.raises(ValueError,match='coupled'):b20.validate_export(p,fake_packer[0])


def test_export_unresolved_beta(fake_packer):
    p=valid_payload();p['report_betas']['lambdaH']='BetaTerm[Quartic,1]'
    with pytest.raises(ValueError):b20.validate_export(p,fake_packer[0])


def test_export_rejects_wrong_state(fake_packer):
    state,_=fake_packer;state.representation.alpha=-1
    with pytest.raises(ValueError):b20.validate_export(valid_payload(),state)


def test_live_and_zero_ode(fake_packer,monkeypatch):
    layout=fake_packer[1]
    def toy(t,y,*,layout,betas):
        ret = np.zeros(y.size)
        ret[0]=-0.05*y[0]+0.2*y[-2]**2
        ret[-2]=-0.1*y[-2]
        ret[-1]=-0.1*y[-1]
        return ret
    monkeypatch.setattr(b20,'derivative',toy)
    v=np.ones(20);v[-2:]=[0.01,0]
    actual,rhs=b20.integrate(v,layout,valid_payload()['report_betas'],points=8)
    z=v.copy();z[-2:]=0
    control,_=b20.integrate(z,layout,valid_payload()['report_betas'],points=8)
    assert actual.success and actual.y.shape==(20,8)
    assert actual.y[-2,-1] != actual.y[-2,0]
    assert np.max(abs(control.y[-2:]))==0
    assert abs(actual.y[0,-1]-control.y[0,-1])>1e-6


def test_integrate_rejects_vector(fake_packer):
    with pytest.raises(ValueError):b20.integrate(np.ones(3),fake_packer[1],{},points=4)


def test_reference_matches(tmp_path):
    path=tmp_path/'original.json'
    obj={'model':b20.MODEL,'status':'SuccessSupportedRenormalisableIntermediateRunning',
         'input_sha256':{'output\\stage33_B_p1_fermion_threshold_handoff.json':'a'*64},
         'intermediate_final_state':{x:1 for x in
            ['mS1Sq','mS2Sq','lambdaS1','lambdaS2','lambdaH1','lambdaH2','lambda12']}}
    path.write_text(json.dumps(obj),encoding='utf-8')
    state=types.SimpleNamespace(**{x:1 for x in obj['intermediate_final_state']})
    assert b20.baseline_control_compare(path,state,'a'*64)['status'] == (
        'ZeroExtraAgreesWithSavedStage34')
    assert b20.baseline_control_compare(path,state,'b'*64)['status'] == (
        'ReferenceStaleOrIncompatible')
    state.lambdaS1=2
    with pytest.raises(ValueError):b20.baseline_control_compare(path,state,'a'*64)


def test_reference_missing(tmp_path):
    assert b20.baseline_control_compare(tmp_path/'missing',None,'x')['status'] == (
        'ReferenceNotAvailable')


def test_wolfram_source_is_native_eft1():
    path=Path(__file__).resolve().parents[2]/'scripts/ExportT3BPlus1EFT1GeneralOneLoop.wl'
    s=path.read_text(encoding='utf-8')
    assert 'T3RGBetaBuildEFT1[2,2,1,1]' in s
    assert 'Get[extraFile]' in s
    assert 'BetaTerm[Quartic,1]' in s
    assert 'T3RGBetaConventionalReportBetas[betas]' in s
    assert 'lambdaT3GeneralQ022000I1' in s


def test_stage34_reused_only_as_checked_input():
    source=(Path(__file__).resolve().parents[2]/
            'scripts/RunT3BPlus1GeneralEFT1.py').read_text(encoding='utf-8')
    assert 's34.verify_boundary(MODEL, Path(folder))' in source
    assert 's34.run_model(' not in source
    assert 'WOLFRAM = ROOT /' in source
