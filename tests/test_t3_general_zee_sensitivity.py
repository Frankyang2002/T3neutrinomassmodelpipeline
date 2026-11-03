"""Pure contract tests: no expensive RG/matching or Wolfram required."""
from copy import deepcopy
import numpy as np
import pytest

from scripts.RunT3GeneralZeeSensitivity import (
    deform, leading_tensor, verify_nominal_invariance, gamma_grid,
    params_from, compensated_quartic_fraction, baseline_consistency,
)


def sample():
    return {'base_state': {'ordinary': {
        'lambdaT3': {'real': 0.1, 'imag': 0.},
        'y1': {'real': [[0.01,0.02,0.03], [0.05,0.04,0.06], [0.08,0.07,0.09]],
               'imag': np.zeros((3,3)).tolist()},
        'y2': {'real': [[0.03,0.08,0.05], [0.02,0.04,0.01], [0.07,0.06,0.09]],
               'imag': np.zeros((3,3)).tolist()},
        'lambdaT3GeneralQ020101I1': {'real': 0.01, 'imag': 0.},
        'muT3GeneralS2Sq': {'real': 100, 'imag': 0.},
    }}}


@pytest.mark.parametrize('direction',('T1','T2','T3'))
@pytest.mark.parametrize('gamma',(0.99,1.,1.01,1.2))
def test_nominal_tensor_and_independent_scalar_invariant(direction,gamma):
    base=sample(); copy=deepcopy(base)
    changed=deform(base,direction,gamma)
    verify_nominal_invariance(base,changed)
    np.testing.assert_allclose(leading_tensor(base),leading_tensor(changed),rtol=1e-13,atol=1e-18)
    assert base==copy  # no mutation
    old=base['base_state']['ordinary']
    new=changed['base_state']['ordinary']
    for name in ('lambdaT3GeneralQ020101I1','muT3GeneralS2Sq'):
        assert old[name] == new[name]


def test_rescaling_directions_distinct():
    b=sample()
    a=deform(b,'T1',4.)['base_state']['ordinary']
    assert a['lambdaT3']['real']==pytest.approx(.4)
    assert a['y1']['real'][0][0]==pytest.approx(.005)
    assert a['y2']['real'][0][0]==pytest.approx(.015)
    a=deform(b,'T2',4.)['base_state']['ordinary']
    assert a['y1']['real'][0][0]==pytest.approx(.01)
    assert a['y2']['real'][0][0]==pytest.approx(.0075)
    a=deform(b,'T3',4.)['base_state']['ordinary']
    assert a['lambdaT3']['real']==pytest.approx(.1)
    assert a['y1']['real'][0][0]==pytest.approx(.04)
    assert a['y2']['real'][0][0]==pytest.approx(.0075)


def test_gamma_grid_and_validation():
    vals=gamma_grid(.99,1.01,11)
    assert len(vals)==11 and vals[5]==1.
    assert all(x<y for x,y in zip(vals,vals[1:]))
    with pytest.raises(ValueError): gamma_grid(.9,1.1,10)
    with pytest.raises(ValueError): gamma_grid(1.,1.1,11)


def test_yukawa_params():
    p=params_from(sample())
    assert len(p)==18 and p['y1_01']==pytest.approx(.02)
    assert p['y2_20']==pytest.approx(.07)


def test_bc_quartic_compensating_rule():
    for gam in (.8,1.,1.2):
        q=compensated_quartic_fraction(.1,.5,.5*gam)
        assert q * abs(.5*gam)==pytest.approx(.05)
    with pytest.raises(ValueError): compensated_quartic_fraction(.5,.5,.1)


def test_baseline_check_rejects_wrong_observable():
    from scripts.RunT3GeneralZeeSensitivity import OBS
    baseline={'prediction':{key:1. for key in OBS}}
    rows=[{'status':'Success',**{key:1. for key in OBS}} for _ in range(3)]
    baseline_consistency(rows,baseline)
    rows[0][OBS[0]]=2.
    with pytest.raises(RuntimeError): baseline_consistency(rows,baseline)

# The continuation campaign is NOT the same schema as the original general fit.
# These tests prevent accidental fallback to the old, non-optimised input.

def test_merged_selection_uses_best_sources(monkeypatch, tmp_path):
    import scripts.RunT3GeneralZeeSensitivity as z
    import json
    rows=[]
    continued={'A_-2','B_+1','C_-1','D_+0','D_+2'}
    for name in z.MODEL_ORDER:
        source=('GeneralFitContinuation' if name in continued else 'OriginalGeneralCampaign')
        run=z.MERGED_SOURCES[source]
        rows.append({'model':name,'source':source,'status':'ExploratoryAcceptedFit',
                     'chi2_final':0.0001,
                     'result_path':str(z.FIT_ROOT.relative_to(z.ROOT) / run / z.tag(name) /
                                       'result.json').replace('/','\\')})
    manifest=tmp_path/'merged_all16.json'
    manifest.write_text(json.dumps({'format':z.MERGED_FORMAT,'rows':rows,
                                    'accepted_count':16,'blocked_count':0}),encoding='utf-8')
    monkeypatch.setattr(z,'MERGED_FILE',manifest)
    selected=z.fit_selections(['A_-2','B_-1','C_-1','D_+0','E_+2'],'merged')
    assert selected['A_-2']['run_id']=='continue_general16_v1'
    assert selected['C_-1']['run_id']=='continue_general16_v1'
    assert selected['D_+0']['run_id']=='continue_general16_v1'
    assert selected['B_-1']['run_id']=='thesis16_general_c5_full_v1'
    assert selected['E_+2']['run_id']=='thesis16_general_c5_full_v1'
    assert selected['C_-1']['manifest']['chi2_final']==pytest.approx(.0001)
    # Explicit old campaign remains possible only when deliberately requested.
    legacy=z.fit_selections(['A_-2'],'thesis16_general_c5_full_v1')
    assert legacy['A_-2']['source']=='ExplicitCampaign'


def test_merged_manifest_rejects_wrong_path_and_poor_fit(monkeypatch,tmp_path):
    import scripts.RunT3GeneralZeeSensitivity as z
    import json
    rows=[]
    for name in z.MODEL_ORDER:
        run='thesis16_general_c5_full_v1'
        rows.append({'model':name,'source':'OriginalGeneralCampaign',
                     'status':'ExploratoryAcceptedFit','chi2_final':1.,
                     'result_path':str(z.FIT_ROOT.relative_to(z.ROOT)/run/z.tag(name)/'result.json')})
    manifest=tmp_path/'merged_all16.json'
    monkeypatch.setattr(z,'MERGED_FILE',manifest)
    payload={'format':z.MERGED_FORMAT,'rows':rows,'accepted_count':16,'blocked_count':0}
    rows[0]['result_path']='output/../../../../wrong.json'
    manifest.write_text(json.dumps(payload),encoding='utf-8')
    with pytest.raises(ValueError,match='expected path'):
        z.fit_selections(['A_-4'],'merged')
    rows[0]['result_path']=str(z.FIT_ROOT.relative_to(z.ROOT)/'thesis16_general_c5_full_v1'/'A_m4'/'result.json')
    rows[0]['chi2_final']=34.5
    manifest.write_text(json.dumps(payload),encoding='utf-8')
    with pytest.raises(ValueError,match='fit screen'):
        z.fit_selections(['A_-4'],'merged')


def test_continued_bc_fractions_and_yukawa_integrity():
    import scripts.RunT3GeneralZeeSensitivity as z
    fit={'matching_and_rge_provenance': {'method':'mixed-scalar-conditional-one-scale',
            'quartic_fraction':0.1,'mu2_fraction':0.01}}
    assert z._bc_input_fractions(fit)==pytest.approx((.1,.01))
    assert z._bc_input_fractions({'extra_quartic_fraction':.2,
                                 'extra_mixed_scalar_bilinear_fraction':.03})==pytest.approx((.2,.03))
    data=sample()
    fit['best_parameters']=z.params_from(data)
    z._check_fit_yukawas(data,fit)
    fit['best_parameters']['y1_00']*=2
    with pytest.raises(ValueError,match='Fitted Yukawa mismatch'):
        z._check_fit_yukawas(data,fit)


def test_continued_nonmixed_sources_without_legacy_fit_flags(monkeypatch,tmp_path):
    import json
    from types import SimpleNamespace
    import scripts.RunT3GeneralZeeSensitivity as z
    monkeypatch.setattr(z,'FIT_ROOT',tmp_path/'fits')
    files={k:tmp_path/(k+'.json') for k in
           ('config','fit','overrides','uv_beta','eft1_beta','matching','matching_status')}
    # The underlying general Matchete source is text, not JSON.
    files['matching'].write_text('general coefficient',encoding='utf-8')
    cfg=sample()
    cfg['general_fit_continuation']={'source_model':'A_-2',
        'source_campaign':'thesis16_general_c5_full_v1',
        'source_scalar_boundary_modified':False}
    fitted=z.params_from(cfg)
    fit={'model':'A_-2','chi2_final':1e-11,'best_parameters':fitted,
        'physical_general_c5_certified':False,
        'full_general_sequential_matching_performed':False,
        'matching_and_rge_provenance':{'method':'full-general-matchete-one-scale-channels',
            'modified_scalar_boundary':False,
            'general_symbolic_provenance':{'source':'canonical_full_general_symbolic',
                'model':'A_-2','run_id':'full16_general_v1'}}}
    def entries(mat):
        return [[{'real':float(x),'imag':0.} for x in row] for row in mat]
    overrides={'lambdaT3':cfg['base_state']['ordinary']['lambdaT3'],
        'y1':entries(cfg['base_state']['ordinary']['y1']['real']),
        'y2':entries(cfg['base_state']['ordinary']['y2']['real'])}
    for name,obj in [('config',cfg),('fit',fit),('overrides',overrides),
                     ('uv_beta',{'report_betas':{k:'beta' for k in overrides}}),
                     ('eft1_beta',{}),('matching_status',{'status':'Matched'})]:
        files[name].write_text(json.dumps(obj),encoding='utf-8')
    monkeypatch.setattr(z,'source_paths',lambda *args:files.copy())
    sel={'run_id':'continue_general16_v1','source':'GeneralFitContinuation',
         'manifest':{'chi2_final':1e-11}}
    config,record,paths=z.sources('A_-2',SimpleNamespace(canonical_run_id='full16_general_v1'),sel)
    assert record['chi2_final']==1e-11
    overrides['y1'][0][0]['real']=999
    files['overrides'].write_text(json.dumps(overrides),encoding='utf-8')
    with pytest.raises(ValueError,match='UV override/Yukawa'):
        z.sources('A_-2',SimpleNamespace(canonical_run_id='full16_general_v1'),sel)


def test_both_general_matching_status_schemas():
    import scripts.RunT3GeneralZeeSensitivity as z
    z._verify_conditional_matching({'full_general_sequential_matching_performed':False})
    z._verify_conditional_matching({'full_sequential_general_scalar_matching_performed':False})
    with pytest.raises(ValueError,match='certification'):
        z._verify_conditional_matching({'full_general_sequential_matching_performed':True})
    with pytest.raises(ValueError,match='certification'):
        z._verify_conditional_matching({})
