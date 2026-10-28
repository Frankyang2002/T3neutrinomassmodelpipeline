"""Offline tests for Batch 18. Real RGBeta replay must be run on Windows."""
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace
import sys
import types

import numpy as np
import pytest
from scripts import ValidateT3BCGeneralScalarSixFits as v


def target():
    return SimpleNamespace(observable_names=v.OBS,
        central_values=np.zeros(5),covariance=np.eye(5))


def zero_pred():return dict.fromkeys(v.OBS,0.)

def goodc5():return {'real':np.eye(3).tolist(),'imag':np.zeros((3,3)).tolist()}


def makehashes(tmp_path):
    d={}
    for i in range(4):
        p=tmp_path/f'input{i}.txt';p.write_text(f'input{i}')
        d[p.name]=hashlib.sha256(p.read_bytes()).hexdigest()
    return d


def base_cfg(model,lam):
    return {'representation':{'d_s1':2,'d_s2':2,
       'd_f':1 if model=='B_-1' else 3,'alpha':-1,'shared_scalar':False},
       'base_state':{'ordinary':{'lambdaT3':{'real':lam,'imag':0.},
          'mS1Sq':1e6,'mS2Sq':1.44e6}}}


def uv_data(tmp_path,model='B_-1',run_id='sample'):
    case='small_uv_running_5e-9'
    p,c=v.expected_paths(tmp_path,case,model,run_id,'uv');p.parent.mkdir(parents=True)
    params={'y1_00':0.01,'log10_lambdaT3':-4.}
    r={'format':'t3-bc-uv-hierarchy-constrained-sobol-refit-v1',
       'case':case,'model':model,'run_id':run_id,'status':'CandidateConstrainedAccepted',
       'general_scalar_sequential_matching_performed':False,
       'physical_matchete_rgbeta_map_certified':False,
       'historical_outputs_overwritten':False,
       'quartic_fraction':0.1,'mu2_fraction':0,
       'mass_screen':{'rhoF_max':0.5},
       'best_mass_safe':{'screen_pass':True,'chi2':0.,'prediction':zero_pred(),
           'rhoF':0.01,'scalar_light_muS_tev':1.2,'scalar_heavy_muS_tev':2.4,
           'parameters':params,'lambdaT3_uv':1e-4,
           'candidate_C5_muS_gev_inv':goodc5(),
           'candidate_C5_muLow_gev_inv':goodc5()},
       'source_sha256':makehashes(tmp_path)}
    p.write_text(json.dumps(r));c.write_text(json.dumps(base_cfg(model,1e-4)))
    return p,c


def frozen_data(tmp_path,model='B_-1',run_id='sample'):
    case='small_uv_frozen_5e-9'
    p,c=v.expected_paths(tmp_path,case,model,run_id,'frozen');p.parent.mkdir(parents=True)
    r={'format':'t3-bc-general-scalar-exploratory-yukawa-refit-v1',
       'case':case,'model':model,'run_id':run_id,'status':'ExploratoryAcceptedFit',
       'passes_approx_chi2_5p86':True,'UV_running_enabled':False,
       'full_sequential_general_scalar_matching_performed':False,
       'physical_matchete_rgbeta_map_certified':False,
       'historical_outputs_overwritten':False,'chi2_final':0.,'prediction':zero_pred(),
       'best_parameters':{'y1_00':.01},
       'final_diagnostics':{'trajectory':{'eft1_end':{
           'mS1Sq_gev2':1e6,'mS2Sq_gev2':1.44e6,
           'quadratic_mu2_gev2':{'real':0,'imag':0}}},
           'MF_muF_gev':[1e5,1.1e5,1.2e5],
           'C5_muS_gev_inv':goodc5(),'C5_low_gev_inv':goodc5()},
       'source_manifest':{'files_sha256':makehashes(tmp_path),
                        'fit_options':{'mu2_fraction':0,'quartic_fraction':.1}}}
    p.write_text(json.dumps(r));c.write_text(json.dumps(base_cfg(model,5e-9)))
    return p,c


def test_registry():
    assert len(v.CASES)==6
    assert len(set((a,b) for a,b,_,_ in v.CASES))==6
    assert sum(mode=='uv' for *_,mode in v.CASES)==2

@pytest.mark.parametrize('model,kind', [('B_-1','uv'),('C_-1','uv'),('B_-1','frozen'),('C_-1','frozen')])
def test_file_naming(model,kind,tmp_path):
    case='small_uv_running_5e-9' if kind=='uv' else 'small_uv_frozen_5e-9'
    path,config=v.expected_paths(tmp_path,case,model,'abc',kind)
    assert path.parent==config.parent
    assert path.suffix=='.json' and config.suffix=='.json'
    assert path.name==('fit.json' if kind=='uv' else ('B_m1' if model=='B_-1' else 'C_m1')+'_exploratory_fit.json')

@pytest.mark.parametrize('bad', [float('nan'),float('inf'),-1.])
def test_invalid_positive_scalar(bad):
    with pytest.raises(ValueError):v.frozen_masses({'mS1Sq_gev2':bad,
        'mS2Sq_gev2':1e6,'quadratic_mu2_gev2':{'real':0,'imag':0}})


def test_mixed_eigenmasses():
    a,b=v.frozen_masses({'mS1Sq_gev2':1e6,'mS2Sq_gev2':4e6,
                        'quadratic_mu2_gev2':{'real':3e5,'imag':4e5}})
    assert a>0 and a<1 and b>2
    with pytest.raises(ValueError):v.frozen_masses({'mS1Sq_gev2':1e6,'mS2Sq_gev2':1e6,
                    'quadratic_mu2_gev2':{'real':1e7,'imag':0}})


def test_nufit_recompute():
    prediction=dict(zip(v.OBS,[1.,2.,3.,4.,5.]))
    assert v.checked_chi2(prediction,55.,target())==55.
    with pytest.raises(ValueError,match='chi2'):v.checked_chi2(prediction,54.,target())


def test_missing_or_drifted_source(tmp_path):
    hashes=makehashes(tmp_path)
    assert v.verify_hashes(tmp_path,hashes)==4
    (tmp_path/'input0.txt').write_text('changed')
    with pytest.raises(ValueError,match='source changed'):v.verify_hashes(tmp_path,hashes)
    with pytest.raises(ValueError):v.verify_hashes(tmp_path,{})

@pytest.mark.parametrize('model',['B_-1','C_-1'])
def test_uv_saved_validation(tmp_path,monkeypatch,model):
    path,c=uv_data(tmp_path,model)
    monkeypatch.setattr(v,'config_yukawa',lambda cfg:{'y1_00':.01})
    out=v.validate_saved(tmp_path,'small_uv_running_5e-9',model,'sample','uv',target())
    assert out['chi2']==0 and out['rhoF']==.01
    assert out['sources_verified']==4
    r=json.loads(path.read_text());r['best_mass_safe']['rhoF']=.6
    path.write_text(json.dumps(r))
    with pytest.raises(ValueError,match='rhoF'):v.validate_saved(tmp_path,'small_uv_running_5e-9',model,'sample','uv',target())

@pytest.mark.parametrize('model',['B_-1','C_-1'])
def test_frozen_saved_validation(tmp_path,monkeypatch,model):
    frozen_data(tmp_path,model)
    monkeypatch.setattr(v,'config_yukawa',lambda cfg:{'y1_00':.01})
    out=v.validate_saved(tmp_path,'small_uv_frozen_5e-9',model,'sample','frozen',target())
    assert math.isclose(out['scalar_heavy_muS_tev'],1.2)
    assert out['sources_verified']==4


def test_uv_lambda_mutation_is_detected(tmp_path,monkeypatch):
    p,c=uv_data(tmp_path);monkeypatch.setattr(v,'config_yukawa',lambda cfg:{'y1_00':.01})
    cfg=json.loads(c.read_text());cfg['base_state']['ordinary']['lambdaT3']['real']=2e-4
    c.write_text(json.dumps(cfg))
    with pytest.raises(ValueError,match='lambda mismatch'):
        v.validate_saved(tmp_path,'small_uv_running_5e-9','B_-1','sample','uv',target())


def test_vector_coupling_magnitudes():
    blocks=[SimpleNamespace(name='lambdaH',kind='real_scalar',start=0,stop=1),
       SimpleNamespace(name='lambdaExtra',kind='complex_scalar',start=1,stop=3),
       SimpleNamespace(name='y1',kind='complex_matrix',start=3,stop=5,shape=(1,1))]
    out=v.couplings_of_vector(np.array([.2,3,4,2,0]),SimpleNamespace(blocks=blocks))
    assert out['max_abs_quartic']==5.
    assert out['max_yukawa_singular']==2.


def test_trace_record():
    row=dict(case='abc',model='B_-1',run_id='v1')
    s={'mu_gev':1000,'m_scalar_light_gev':1000,'m_scalar_heavy_gev':1500,
       'min_fermion_mass_gev':1e5,'max_scalar_over_min_fermion':0.015,'mu2_abs_gev2':0.}
    r=v.trace_row(row,'EFT1',s,{'max_abs_quartic':.1,'max_yukawa_singular':.2})
    assert tuple(r.keys()) and set(r)==set(v.TRACE_COLS)
    assert r['scalar_heavy_tev']==1.5


def test_summary_saved_only_distinct():
    rows=[dict(case='small_uv_running_5e-9',model='B_-1',chi2=1.,rhoF=.02,
               scalar_light_muS_tev=1.2,scalar_heavy_muS_tev=1.5,lambdaT3_uv=1e-5)]
    s=v.summary_md(rows,replay=False)
    assert 'no RG trajectory verification' in s
    assert 'not 48/48' in s


def test_cli_and_execution_writes_one_new_tree(tmp_path,monkeypatch):
    # Mock the import that validates run-id and the source-facing live methods.
    stub=types.ModuleType('scripts.RunT3GeneralScalarBCPilot')
    stub.parse_run_id=lambda s:s if s.isidentifier() else (_ for _ in ()).throw(ValueError())
    monkeypatch.setitem(sys.modules,stub.__name__,stub)
    monkeypatch.setattr(v,'load_target',lambda root:target())
    monkeypatch.setattr(v,'CASES',(('small_uv_running_5e-9','B_-1','runid','uv'),))
    def fake_validate(*args):
        return dict(case='small_uv_running_5e-9',model='B_-1',run_id='runid',
          kind='uv',chi2=.1,prediction=zero_pred(),rhoF=.02,
          scalar_light_muS_tev=1.,scalar_heavy_muS_tev=2.,lambdaT3_uv=.00001,
          fit_path=tmp_path/'fit.json',cfg_path=tmp_path/'cfg.json',_saved_fit={})
    monkeypatch.setattr(v,'validate_saved',fake_validate)
    called=[]
    def replay(*args,**kw):
        called.append(True)
        return {'replay_status':'CandidateReplayAgrees','max_abs_quartic':.2,
                'max_yukawa_singular':.01,'quartic_coordinate_screen_pass':True,
                'yukawa_singular_screen_pass':True},[dict.fromkeys(v.TRACE_COLS,1)]
    monkeypatch.setattr(v,'replay_point',replay)
    r=v.execute(tmp_path,'newfit',grid=4,saved_only=False)
    assert called and r['status']=='ReplayValidated'
    dst=tmp_path/v.OUT_BASE/'newfit'
    assert (dst/'report.json').exists() and (dst/'trajectory_samples.csv').exists()
    with pytest.raises(ValueError,match='already exists'):v.execute(tmp_path,'newfit',grid=4)
    s=v.execute(tmp_path,'savedfit',grid=4,saved_only=True)
    assert s['status']=='SavedOnlyValidated'
    assert not (tmp_path/v.OUT_BASE/'savedfit'/'trajectory_samples.csv').exists()


def test_execute_fail_closed_before_output(tmp_path,monkeypatch):
    stub=types.ModuleType('scripts.RunT3GeneralScalarBCPilot')
    stub.parse_run_id=lambda s:s
    monkeypatch.setitem(sys.modules,stub.__name__,stub)
    monkeypatch.setattr(v,'load_target',lambda root:target())
    monkeypatch.setattr(v,'validate_saved',lambda *a,**k:(_ for _ in ()).throw(ValueError('broken')))
    with pytest.raises(ValueError,match='broken'):v.execute(tmp_path,'badfit')
    assert not (tmp_path/v.OUT_BASE/'badfit').exists()


def test_cli_main_routes_options(monkeypatch):
    captured={}
    monkeypatch.setattr(v,'execute',lambda root,run_id,**opts:captured.update(
        root=root,run_id=run_id,**opts))
    monkeypatch.setattr(sys,'argv',['ValidateT3BCGeneralScalarSixFits',
        '--run-id','sample_v1','--grid','14','--saved-only'])
    v.main()
    assert captured['run_id']=='sample_v1'
    assert captured['grid']==14 and captured['saved_only'] is True


def test_complex_c5_symmetry_and_nonzero():
    assert v.check_c5_matrix(goodc5(),'C5')==pytest.approx(math.sqrt(3))
    bad=goodc5();bad['real'][0][1]=1
    with pytest.raises(ValueError,match='invalid complex symmetric'):
        v.check_c5_matrix(bad,'C5')
    with pytest.raises(ValueError,match='zero'):
        v.check_c5_matrix({'real':np.zeros((3,3)).tolist(),
                           'imag':np.zeros((3,3)).tolist()},'C5')


def test_negative_and_nonfinite_chi2_rejected(tmp_path,monkeypatch):
    p,c=uv_data(tmp_path);monkeypatch.setattr(v,'config_yukawa',lambda cfg:{'y1_00':.01})
    for x in [-1, float('inf')]:
        r=json.loads(p.read_text());r['best_mass_safe']['chi2']=x
        p.write_text(json.dumps(r))
        with pytest.raises(ValueError):v.validate_saved(
            tmp_path,'small_uv_running_5e-9','B_-1','sample','uv',target())


def test_frozen_ratio_uses_F_scale_not_S_scale(tmp_path,monkeypatch):
    p,c=frozen_data(tmp_path)
    monkeypatch.setattr(v,'config_yukawa',lambda cfg:{'y1_00':.01})
    cfg=json.loads(c.read_text());cfg['base_state']['ordinary']['mS2Sq']=4e6
    c.write_text(json.dumps(cfg))
    r=v.validate_saved(tmp_path,'small_uv_frozen_5e-9','B_-1','sample','frozen',target())
    assert r['rhoF']==pytest.approx(.02)
    assert r['scalar_heavy_muS_tev']==pytest.approx(1.2)
