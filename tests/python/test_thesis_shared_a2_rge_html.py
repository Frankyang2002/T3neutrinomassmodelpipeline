"""Offline checks for true solver-sampled shared-A2 UV/EFT1 dashboard."""
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from scripts import RunT3SharedA2General16RGE as rge
from scripts import RunT3SharedA2General16Reports as shared


def state(uv=True, mu=-42.5):
    sm = SimpleNamespace(gY=.35,g2=.65,g3=1.1,lambdaH=.13)
    extras={'muT3GeneralS1Sq':complex(mu,0.),
            'y2Conj':np.eye(3)*.02}
    return SimpleNamespace(sm=sm,lambdaT3=.1+0j,mS1Sq=1e6,
       mS2Sq=1.2e6, general_couplings=extras,
       y1=np.eye(3)*1e-4, y2=np.eye(3)*2e-4)


class FakeResult:
    def __init__(self,uv=True):
        self.mu_gev=np.array([1e7,1e6,1e5]) if uv else np.array([1e5,1e4,1e3])
        self.solver_success=True
        self.nfev=23
        self.final_state=state(uv)
    def state_at_index(self,i):
        return state()


def test_gamma_grid_twenty_one_center():
    g=shared.gammas(.8,1.2,21)
    assert len(g)==21 and g[0]==pytest.approx(.8) and g[-1]==pytest.approx(1.2)
    assert g[10]==1.0 and all(a<b for a,b in zip(g,g[1:]))


def test_sampled_real_rge_not_endpoint_interpolation():
    stage=rge.sampled_stage(FakeResult(),uv=True)
    assert stage['mu_gev']==[1e7,1e6,1e5]
    assert stage['series']['y2_norm']==pytest.approx([np.sqrt(3)*2e-4]*3)
    assert stage['series']['muS1_real_gev2']==[-42.5]*3
    assert stage['nfev']==23
    eft=rge.sampled_stage(FakeResult(False),uv=False)
    assert 'y1_norm' not in eft['series']
    assert eft['series']['g2']==[.65,.65,.65]


def test_endpoint_check_detects_changed_ode(tmp_path):
    curve={'uv':rge.sampled_stage(FakeResult(),uv=True),
           'eft1':rge.sampled_stage(FakeResult(False),uv=False)}
    c={'numerical_details': {
       'uv_at_fermion_threshold':{'gY':.35,'g2':.65,'g3':1.1,'lambdaH':.13,'mS1Sq':1e6,'mS2Sq':1.2e6},
       'eft1_at_scalar_threshold':{'gY':.35,'g2':.65,'g3':1.1,'lambdaH':.13,'mS1Sq':1e6,'mS2Sq':1.2e6}}}
    rge.validate_endpoints(curve,c)
    c['numerical_details']['eft1_at_scalar_threshold']['mS1Sq']=2e6
    with pytest.raises(ValueError,match='endpoint'):
        rge.validate_endpoints(curve,c)


def test_generate_writes_real_data_and_labels_BC_unavailable(tmp_path,monkeypatch):
    out=tmp_path
    statuses=[{'model':'A_-2','status':'Success'}, {'model':'B_-1','status':'Success'}]
    mdir=out/'A_m2';mdir.mkdir()
    shared.campaign.save(mdir/'centre.json',{'numerical_details':{'dummy':True}})
    def compute(_model,_params,_settings,_ref,_conj,_centre,n):
        return {'format':rge.FORMAT,'model':'A_-2','status':'Success',
                'n_save_points_per_stage':n,
                'uv':rge.sampled_stage(FakeResult(),uv=True),
                'eft1':rge.sampled_stage(FakeResult(False),uv=False)}
    monkeypatch.setattr(rge,'compute_one',compute)
    assert rge.generate(out,statuses,{}, {},{},'a2-reference',n_points=41)
    assert (out/'interactive_RGE.html').is_file()
    page=(out/'interactive_RGE.html').read_text()
    assert 'A_-2' in page and 'B_-1' in page
    assert 'Actual' not in page or 'saved' in page
    status=json.loads((out/'rge_status.json').read_text())
    assert status['rge_curves_available']==1
    assert status['unavailable'][0]['model']=='B_-1'
    assert not (out/'B_m1'/'rge_trajectories.json').exists()


def test_report_index_links_both_dashboards(tmp_path):
    shared.report_html([{'model':'A_-2','status':'Success'}],tmp_path,
        donor_run='fit',gamma_grid=[.8,1.,1.2],oval=True,rge=True)
    text=(tmp_path/'index.html').read_text()
    assert 'sensitivity_oval.html' in text
    assert 'interactive_RGE.html' in text


def test_no_fake_weinberg_rge_in_source():
    assert 'c5_norm' not in rge.SERIES
    assert rge.FORMAT.endswith('trajectories-v1')
