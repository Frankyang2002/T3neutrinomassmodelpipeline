"""Offline checks of result classifications, provenance checks and figure creation."""
from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path

import pytest

SCRIPT=Path(__file__).resolve().parents[2]/'scripts/BuildT3ReducedScalarFinalResults.py'
spec=importlib.util.spec_from_file_location('finalreport',SCRIPT)
r=importlib.util.module_from_spec(spec)
spec.loader.exec_module(r)

MODELS=('A_-2','A_+0','A_-4','B_+1','B_-1','B_-3','C_+1','C_-1',
        'C_-3','D_+0','D_+2','D_-2','E_+0','E_+2','E_-2','E_-4')
SMALL_SHARED=(1e-14,.008,11518,1.47,37,12171,9273,9750,17818,.0002,232,.013,3922,2933,3915,15394)
LARGE_SHARED=(1e-14,.0075,11464,2.187,35,12124,9256,9739,17797,.0002,239,.012,4106,3106,4099,15437)


def source():
    rows=[]
    for case,shared in zip(r.FROZEN_CASES,(SMALL_SHARED,LARGE_SHARED)):
        for i,m in enumerate(MODELS):
            a=1e-13 if m=='A_-2' else .0001
            b=a if m=='A_-2' else shared[i]
            row={'model':m,'case':case,'status':'Accepted',
                 'shared_status':'NeutrinoEvaluated',
                 'chi2_optimised':a,'chi2_shared_a2':b,
                 'boundary_scale_gev':1e5,
                 'boundary_lambdaT3':(5e-9 if case==r.FROZEN_CASES[0] else .5)*(.5 if m.startswith('E_') else 1.),
                 'c5_source_kind':'hierarchical-hard-plus-direct'}
            for k in r.OBS:
                row['fit_'+k]=.1
                row['shared_'+k]=.1
            rows.append(row)
    return {'format':r.FROZEN_FORMAT,'cases':list(r.FROZEN_CASES),'rows':rows,
            'settings':{'a2_source':'new-fit','a2_as_start':True,
                        'extra_quartic_boundary_real_imag':[0.,0.],
                        'extra_quadratic_boundary':0.}}


def test_complete_frozen_counts_and_shared_acceptance():
    rows=r.validate_frozen_payload(source(),MODELS)
    assert len(rows)==32
    assert sum(float(z['chi2_shared_a2'])<=r.FIT_CUTOFF for z in rows.values())==10
    assert set(rows)=={(c,m) for c in r.FROZEN_CASES for m in MODELS}


@pytest.mark.parametrize('field,value',[
    ('a2_source','historical'),('a2_as_start',False),
    ('extra_quadratic_boundary',.01),
    ('extra_quartic_boundary_real_imag',[.01,0.])])
def test_reject_unjustified_boundary_or_wrong_a2(field,value):
    bad=source();bad['settings'][field]=value
    with pytest.raises(ValueError):r.validate_frozen_payload(bad,MODELS)


def test_reject_omitted_or_double_counted_models():
    bad=source();bad['rows'].pop()
    with pytest.raises(ValueError):r.validate_frozen_payload(bad,MODELS)
    bad=source();bad['rows'][0]['model']='E_-4'
    with pytest.raises(ValueError):r.validate_frozen_payload(bad,MODELS)


def test_reject_bad_chi2():
    bad=source();bad['rows'][1]['chi2_optimised']=9
    with pytest.raises(ValueError):r.validate_frozen_payload(bad,MODELS)
    bad=source();bad['rows'][1]['chi2_shared_a2']=None
    with pytest.raises(ValueError):r.validate_frozen_payload(bad,MODELS)


def test_reject_wrong_class_e_coupling_or_scale():
    bad=source();bad['rows'][12]['boundary_lambdaT3']=5e-9
    with pytest.raises(ValueError):r.validate_frozen_payload(bad,MODELS)
    bad=source();bad['rows'][1]['boundary_scale_gev']=1e7
    with pytest.raises(ValueError):r.validate_frozen_payload(bad,MODELS)


def test_extract_running_endpoint_from_historical_seed_only():
    d={'format':r.FROZEN_FORMAT,'model':'A_-2','case':r.RUNNING_CASE,
       'status':'NoMassSafeNeutrinoFit',
       'sobol_trials':[{'label':'historical_seed','chi2':None,
                        'status':'Rejected','scalar_masses_muS_tev':[40.0,41.0],
                        'rhoF':.38,'screen':False,'max_scalar_quartic_abs':.338},
                       {'label':'sobol_1','chi2':None,'status':'Rejected',
                        'scalar_masses_muS_tev':[60,61],'rhoF':.4}]}
    z=r.extract_running_endpoint(d,'A_-2')
    assert z['mS1_formal_muS_tev']==40.
    assert z['number_of_initial_trials']==2
    assert z['screen_pass']=='false'
    d['sobol_trials'][0]['chi2']=1.
    with pytest.raises(ValueError):r.extract_running_endpoint(d,'A_-2')


def test_per_model_sha_provenance_and_chi2(tmp_path):
    rows=r.validate_frozen_payload(source(),MODELS)
    for (case,model),z in rows.items():
        p=tmp_path/case/r.model_safe(model)/'fit.json'
        p.parent.mkdir(parents=True,exist_ok=True)
        p.write_text(json.dumps({'format':r.FROZEN_FORMAT,'model':model,'case':case,
          'status':'Accepted','extra_quadratic_boundary_exact_zero':True,
          'quartic_one_loop_running_included':True,
          'general_scalar_one_loop_matching_included':False,
          'best_mass_safe':{'chi2':z['chi2_optimised']},
          'shared_a2':{'chi2':z['chi2_shared_a2']}}),encoding='utf-8')
    paths=r.validate_frozen_fit_files(tmp_path,rows)
    assert len(paths)==32
    assert len(r.checksum(next(iter(paths.values()))))==64


def test_observable_csv_and_plot_outputs(tmp_path):
    rows=r.validate_frozen_payload(source(),MODELS)
    target={k:.1 for k in r.OBS};sig={k:.01 for k in r.OBS}
    items=r.make_observable_rows(rows,target,sig)
    assert len(items)==320
    assert {z['yukawa_choice'] for z in items}=={'individual','shared_A_-2'}
    r.write_csv(tmp_path/'pulls.csv',items)
    with (tmp_path/'pulls.csv').open(newline='') as f:assert len(list(csv.DictReader(f)))==320
    r.plot_frozen(rows,MODELS,tmp_path)
    running={m:{'mS1_formal_muS_tev':40.0+i,
               'mS2_formal_muS_tev':41.0+i,'rhoF_at_muF':.38} for i,m in enumerate(MODELS)}
    r.plot_running_endpoint(running,MODELS,tmp_path)
    assert (tmp_path/'frozen_chi2_comparison.png').stat().st_size>10000
    assert (tmp_path/'running_uv_formal_mass_endpoints.png').stat().st_size>10000


def test_uv_trajectory_plot_uses_real_sampled_rows(tmp_path):
    import numpy as np
    t=[]
    for model in MODELS:
        for mu in np.geomspace(1e7,1e5,80):
            t.append({'model':model,'mu_gev':float(mu),'mS1Sq_tev2':4.,'mS2Sq_tev2':5.,
                      'mS1_mass_tev_if_positive':2.,'mS2_mass_tev_if_positive':5**.5})
    r.plot_uv_trajectories(t,MODELS,tmp_path)
    assert (tmp_path/'running_uv_scalar_trajectories.png').stat().st_size>10000


def test_uv_beta_reused_b_c_symbols_not_dropped(tmp_path):
    """Regression: B/C Stage31 exports contain six already-in-state quartics."""
    from types import SimpleNamespace
    for model in ('B_-1','C_-1'):
        reused=tuple('old_'+str(i) for i in range(6))
        new=('new_quartic',)
        layout=SimpleNamespace(names=('mS1Sq','mS2Sq',*reused))
        entry={'model':model,'dimensions':(2,2,1 if model=='B_-1' else 3),'alpha':-1}
        beta_file=tmp_path/('stage31_'+r.model_safe(model)+'_full_uv_betas.json')
        beta_file.write_text(json.dumps({'model':model,'new_quartic_symbols':list(new),
                          'reused_quartic_symbols':list(reused)}),encoding='utf-8')
        seen=[]
        def validate(payload,received_entry,received_new,received_reused,base_names):
            seen.append((received_reused,tuple(base_names)))
            assert payload['model']==model
            assert set(payload['reused_quartic_symbols'])==set(received_reused)
            assert tuple(received_new)==new
            assert tuple(base_names)==layout.names
            return {'mS1Sq':'dummy'}
        def fail_old(*args):
            pytest.fail('Old bulk UV beta helper was called for reused quartics')
        bulk=SimpleNamespace(BETA_DIR=tmp_path,safe=r.model_safe,
                             read=r.read_json,uv_betas=fail_old)
        s31=SimpleNamespace(validate_export=validate)
        betas,path=r.select_uv_betas(model,entry,new,reused,layout,bulk=bulk,s31=s31)
        assert path==beta_file
        assert betas=={'mS1Sq':'dummy'}
        assert seen==[(reused,layout.names)]


def test_uv_beta_non_reused_preserves_existing_export_path(tmp_path):
    """All other models continue to use the exact previous Stage30/31 helper."""
    from types import SimpleNamespace
    calls=[]
    def existing(model,entry,new,layout):
        calls.append((model,new,layout.names))
        return {'g2':'beta'},tmp_path/'preexisting.json'
    bulk=SimpleNamespace(uv_betas=existing)
    layout=SimpleNamespace(names=('g2','mS1Sq','mS2Sq'))
    entry={'model':'A_+0'}
    betas,path=r.select_uv_betas('A_+0',entry,('extra',),(),layout,
                                  bulk=bulk,s31=None)
    assert calls==[('A_+0',('extra',),layout.names)]
    assert betas=={'g2':'beta'}
    assert path==tmp_path/'preexisting.json'


def test_uv_beta_reused_unexpected_model_fails_closed(tmp_path):
    from types import SimpleNamespace
    with pytest.raises(ValueError,match='unexpected reused UV quartics'):
        r.select_uv_betas('A_-2',{'model':'A_-2'},(),('unexpected',),
             SimpleNamespace(names=('unexpected',)),bulk=SimpleNamespace(),s31=SimpleNamespace())
