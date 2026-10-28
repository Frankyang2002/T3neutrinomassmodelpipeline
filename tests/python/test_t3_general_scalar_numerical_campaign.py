"""Offline tests: 48-case inventory, explicit support gate, and preservation.

The real RGBeta/Matchete integrations are executed only by the Windows user.
"""
from __future__ import annotations

from pathlib import Path
import json
import math
import sys
import types

import pytest

from scripts import RunT3GeneralScalarNumericalCampaign as p


def put(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2)+'\n',encoding='utf-8')


@pytest.fixture
def project(tmp_path):
    totals={'matchete_extra_quartics':50,'rgbeta_preexisting_structural':12,
            'rgbeta_new_vertices':38,'quadratic_couplings_blocked':8}
    put(tmp_path/'RGE/general/generated/general_rgbeta_quartic_registry.json',
        {'totals':totals,'models':{m:{} for m in p.MODEL_NAMES}})
    for case,(root,_,lam) in p.STUDIES.items():
        for model in p.MODEL_NAMES:
            stem=p.historical_stem(model)
            dims=p.MODEL_REP[model[0]]
            rep=dict(zip(('d_s1','d_s2','d_f'),dims))
            rep.update(alpha=int(model.split('_')[1]),shared_scalar=False)
            lt=lam/2 if model[0]=='E' else lam
            put(tmp_path/root/'fitted'/f'{stem}.json',
                {'representation':rep,'base_state':{'ordinary':{
                    'lambdaT3':{'real':lt,'imag':0.0}}}})
            put(tmp_path/root/'fits'/f'{stem}.json',
                {'model':stem,'status':'Success','chi2_fit':0.02})
    return tmp_path


def fake_fit(root,case,model,runid,chi2,*,q=0.1,mu=0.0,continuation=False):
    tag=model.replace('-','m')
    folder=root/p.OUT/case/runid
    base={'case':case,'model':model,'run_id':runid,
          'status':'ExploratoryAcceptedFit' if chi2<=p.CUTOFF else 'ExploratoryPoorFit',
          'chi2_final':chi2,'passes_approx_chi2_5p86':chi2<=p.CUTOFF,
          'physical_matchete_rgbeta_map_certified':False,
          'full_sequential_general_scalar_matching_performed':False}
    provenance={'prior_input_file_sha256' if continuation else 'files_sha256':
                {'RGE/general/generated/general_rgbeta_quartic_registry.json':
                 p.digest(root/'RGE/general/generated/general_rgbeta_quartic_registry.json')}}
    if continuation:
        base['source_manifest']={**provenance,'quartic_fraction':q,'mu2_fraction':mu}
    else:
        base['source_manifest']=provenance
        base['input_extra_quartic_fraction_of_abs_lambdaT3']=q
        base['input_mu2_fraction_of_min_scalar_m2']=mu
    put(folder/f'{tag}_exploratory_fit.json',base)
    put(folder/f'{tag}_exploratory_fitted_config.json',
        {'exploratory_general_scalar_fit':{'model':model,'case':case,'run_id':runid}})


def test_model_inventory_count():
    assert len(p.MODEL_NAMES)==16
    assert set(p.SUPPORTED)=={'B_-1','C_-1'}
    assert len(p.STUDIES)==3
    assert p.historical_stem('B_-1')=='T3_dS1_2_dS2_2_dF_1_alpha_m1'
    assert p.historical_stem('C_-1')=='T3_dS1_2_dS2_2_dF_3_alpha_m1'
    assert p.historical_stem('A_+0')=='T3_dS1_1_dS2_3_dF_2_alpha_p0'
    assert p.historical_stem('E_-4')=='T3_dS1_3_dS2_3_dF_2_alpha_m4'


def test_campaign_id_path_safety():
    for bad in ('','../x','a/b','a b','a'*36,'C:\\bad','☃'):
        with pytest.raises(ValueError):
            p.run_id(bad)
    assert p.run_id('gen26_01')=='gen26_01'


def test_48_case_inventory(project):
    inv=p.inventory(project)
    assert len(inv['rows'])==48
    assert sum(x['status']=='ReadyForExploratoryFit' for x in inv['rows'])==6
    assert sum(x['status']=='BlockedGeneralMatching' for x in inv['rows'])==42
    assert inv['registry']['totals']['rgbeta_new_vertices']==38
    assert not inv['physical_coupling_map_certified']


def test_separate_scales_and_class_e_factor(project):
    inv=p.inventory(project)
    cases={r['case'] for r in inv['rows']}
    assert cases==set(p.STUDIES)
    e=[r for r in inv['rows'] if r['model']=='E_+0']
    assert e[0]['historical_lambdaT3']==2.5e-9
    assert e[1]['historical_lambdaT3']==2.5e-9
    assert e[2]['historical_lambdaT3']==0.25
    assert [x['historical_uv_running'] for x in e]==[True,False,False]


def test_preserve_poor_historical_fit(project):
    case='large_uv_frozen_0p5';model='B_-1';stem=p.historical_stem(model)
    path=project/p.STUDIES[case][0]/'fits'/f'{stem}.json'
    put(path,{'model':stem,'status':'PoorFit','chi2_fit':32.294})
    row=p.historical_row(project,case,model)
    assert row['historical_chi2']==32.294
    assert row['historical_fit_screen_pass'] is False
    assert row['status']=='ReadyForExploratoryFit'


def test_detect_existing_accepted_batch11(project):
    case='small_uv_frozen_5e-9';model='B_-1'
    fake_fit(project,case,model,'b_small_fit_v1',3e-12)
    inv=p.inventory(project)
    row=next(r for r in inv['rows'] if r['case']==case and r['model']==model)
    assert row['status']=='ExistingExploratoryAccepted'
    assert row['existing_fit']['chi2']==3e-12
    assert row['existing_fit']['run_id']=='b_small_fit_v1'


def test_detect_existing_accepted_continuation(project):
    case='large_uv_frozen_0p5';model='C_-1'
    fake_fit(project,case,model,'c_large_fit_v2',4.58e-11,continuation=True)
    inv=p.inventory(project)
    row=next(r for r in inv['rows'] if r['case']==case and r['model']==model)
    assert row['status']=='ExistingExploratoryAccepted'
    assert row['existing_fit']['run_id']=='c_large_fit_v2'


def test_wrong_fraction_or_unaccepted_not_reused(project):
    case='large_uv_frozen_0p5';model='B_-1'
    fake_fit(project,case,model,'wrong_q',0.1,q=0.2)
    fake_fit(project,case,model,'bad_chi2',100.0)
    assert p.available_fits(project,case,model)==[]
    fake_fit(project,case,model,'right_q',1.0)
    assert len(p.available_fits(project,case,model))==1


def test_new_ids_unique_and_safe():
    names=[name for model in p.SUPPORTED for name in p.new_ids('general_v1',model)]
    assert len(names)==len(set(names))
    assert all(len(name)<=64 for name in names)
    assert all(p.run_id(name[:35]) for name in names)


def test_preflight_refuses_existing_campaign(project):
    inv=p.inventory(project)
    d=project/p.OUT/'campaigns'/'batch'
    d.mkdir(parents=True)
    with pytest.raises(ValueError,match='already exists'):
        p.validate_outputs_free(project,'batch',inv,fresh=False)


def test_preflight_refuses_preexisting_new_fit(project):
    inv=p.inventory(project)
    fid,_=p.new_ids('batch','B_-1')
    (project/p.OUT/'small_uv_running_5e-9'/fid).mkdir(parents=True)
    with pytest.raises(ValueError,match='already exists'):
        p.validate_outputs_free(project,'batch',inv,fresh=True)


def test_invalid_registry_fails_closed(project):
    put(project/'RGE/general/generated/general_rgbeta_quartic_registry.json',
        {'models':{m:{} for m in p.MODEL_NAMES},'totals':{'rgbeta_new_vertices':37}})
    with pytest.raises(ValueError,match='inventory unexpectedly changed'):
        p.inventory(project)


def test_missing_historical_input_not_mislabelled(project):
    case='small_uv_running_5e-9';model='A_+0'
    path=project/p.STUDIES[case][0]/'fits'/f'{p.historical_stem(model)}.json'
    path.unlink()
    inv=p.inventory(project)
    row=next(r for r in inv['rows'] if r['case']==case and r['model']==model)
    assert row['status']=='MissingHistoricalInput'
    with pytest.raises(ValueError,match='missing historical input'):
        p.validate_outputs_free(project,'batch',inv,fresh=False)


def test_report_new_only_and_42_blocks(project):
    inv=p.inventory(project)
    folder=project/p.OUT/'campaigns'/'report1';folder.mkdir(parents=True)
    original=project/p.STUDIES['small_uv_running_5e-9'][0]/'fits'/f'{p.historical_stem("B_-1")}.json'
    before=p.digest(original)
    p.write_report(project,'report1',inv)
    assert (folder/'campaign.json').exists()
    s=(folder/'summary.md').read_text(encoding='utf-8')
    assert 'Other 42 cells remain explicitly blocked' in s
    assert s.count('BlockedGeneralMatching')>=42
    assert p.digest(original)==before


def test_run_campaign_refuses_invalid_options_before_writes(project):
    with pytest.raises(ValueError):
        p.run_campaign(project,'batch',quartic_fraction=-0.2)
    with pytest.raises(ValueError):
        p.run_campaign(project,'batch',sobol_points=7)
    assert not (project/p.OUT/'campaigns'/'batch').exists()


def test_report_flags_uncertified(project):
    inv=p.inventory(project)
    assert inv['sequential_general_scalar_matching_performed'] is False
    assert inv['physical_coupling_map_certified'] is False
    assert inv['supported_cell_count']==6
    assert inv['general_c5_not_implemented_cell_count']==42


def test_accepted_existing_can_block_new_output_allocation(project):
    case='large_uv_frozen_0p5';model='C_-1'
    fake_fit(project,case,model,'accepted',1e-12)
    inv=p.inventory(project)
    fit,_=p.new_ids('batch',model)
    (project/p.OUT/case/fit).mkdir(parents=True)
    # Ordinary resume-like campaign reuses accepted result, does not touch fit path.
    # But fresh mode must fail rather than overwrite it.
    with pytest.raises(ValueError,match='already exists'):
        p.validate_outputs_free(project,'batch',inv,fresh=True)


def test_reject_wrong_representation(project):
    case='small_uv_frozen_5e-9';model='B_-1'
    path=project/p.STUDIES[case][0]/'fitted'/f'{p.historical_stem(model)}.json'
    put(path,{'representation':{'d_s1':1},'base_state':{'ordinary':{'lambdaT3':{'real':5e-9,'imag':0}}}})
    with pytest.raises(ValueError,match='representation mismatch'):
        p.inventory(project)


def test_chosen_existing_is_smallest_chi2(project):
    case='small_uv_running_5e-9';model='B_-1'
    fake_fit(project,case,model,'v1',3.0)
    fake_fit(project,case,model,'v2',0.1,continuation=True)
    rows=p.available_fits(project,case,model)
    assert [r['run_id'] for r in rows]==['v2','v1']


def install_fake_fit_modules(monkeypatch,root,*,initial_accepted=True):
    """Exercise campaign state machine without RG/Matchete dependencies.

    Mock both sys.modules AND package attributes: `from scripts import module`
    can otherwise retrieve the already-imported genuine fitter from scripts.
    The production same-checkout guard is deliberately retained.
    """
    import scripts as scripts_package
    name1='scripts.RunT3GeneralScalarBCExploratoryFit'
    name2='scripts.ContinueT3GeneralScalarBCExploratoryFit'
    r=types.ModuleType(name1);r.ROOT=root
    c=types.ModuleType(name2)
    calls=[]

    def simulate(case,model,identifier,*,accepted):
        chi=1e-8 if accepted else 100.
        fake_fit(root,case,model,identifier,chi)
        return {'case':case,'model':model,'run_id':identifier,
                'chi2_final':chi,'passes_approx_chi2_5p86':accepted}

    def new_run(case,model,identifier,**kwargs):
        calls.append(('fit',case,model,identifier))
        return simulate(case,model,identifier,accepted=initial_accepted)

    def continue_run(case,model,previous_run_id,identifier,**kwargs):
        calls.append(('continue',case,model,identifier))
        return simulate(case,model,identifier,accepted=True)

    r.new_run=new_run
    c.continue_run=continue_run
    monkeypatch.setitem(sys.modules,name1,r)
    monkeypatch.setitem(sys.modules,name2,c)
    monkeypatch.setattr(scripts_package,'RunT3GeneralScalarBCExploratoryFit',r,
                        raising=False)
    monkeypatch.setattr(scripts_package,'ContinueT3GeneralScalarBCExploratoryFit',c,
                        raising=False)
    return calls


def test_campaign_run_reuses_all_six_successes(project,monkeypatch):
    for case in p.STUDIES:
        for model in p.SUPPORTED:
            fake_fit(project,case,model,f'previous_{model[0]}',1e-8)
    calls=install_fake_fit_modules(monkeypatch,project)
    inv=p.run_campaign(project,'reuse_all')
    assert calls==[]
    assert sum(r['status']=='ExistingExploratoryAccepted' for r in inv['rows'])==6
    assert sum(r['status']=='BlockedGeneralMatching' for r in inv['rows'])==42
    assert (project/p.OUT/'campaigns'/'reuse_all'/'campaign.json').exists()
    assert (project/p.OUT/'campaigns'/'reuse_all'/'summary.md').exists()


def test_campaign_continues_failed_first_fit(project,monkeypatch):
    for case in p.STUDIES:
        for model in p.SUPPORTED:
            if (case,model)==('small_uv_running_5e-9','B_-1'):
                continue
            fake_fit(project,case,model,f'previous_{model[0]}',1e-8)
    calls=install_fake_fit_modules(monkeypatch,project,initial_accepted=False)
    inv=p.run_campaign(project,'one_repair')
    assert [x[0] for x in calls]==['fit','continue']
    row=next(r for r in inv['rows'] if r['case']=='small_uv_running_5e-9'
             and r['model']=='B_-1')
    assert row['status']=='NewExploratoryAccepted'
    assert row['new_fit']['run_id']=='one_repair_B_m1_cont'
    assert row['new_fit']['chi2']==1e-8
    assert sum(r['status']=='BlockedGeneralMatching' for r in inv['rows'])==42
    assert (project/p.OUT/'campaigns'/'one_repair'/'campaign.json').exists()


def test_stale_existing_fit_not_reused(project):
    case='small_uv_running_5e-9';model='B_-1'
    fake_fit(project,case,model,'former_fit',1e-10)
    assert len(p.available_fits(project,case,model))==1
    config=project/'RGE/general/generated/general_rgbeta_quartic_registry.json'
    raw=json.loads(config.read_text(encoding='utf-8'))
    raw['extra_metadata']='changed'
    put(config,raw)
    assert p.available_fits(project,case,model)==[]
