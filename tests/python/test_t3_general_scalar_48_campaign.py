"""Offline unit checks for the fast 16x3 campaign, no native Wolfram required."""
from __future__ import annotations

import ast
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/'scripts/RunT3GeneralScalar48Campaign.py'
WOLFRAM=ROOT/'scripts/ExportT3GeneralScalarEFT1Bulk.wl'


def pure():
    """Load the pure functions without importing the expensive project stack."""
    tree=ast.parse(SOURCE.read_text(encoding='utf-8'))
    names={'safe','stem','inventory','row_of','imaginary_zero'}
    nodes=[x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name in names]
    assert {x.name for x in nodes}==names
    namespace={}
    exec(compile(ast.Module(body=nodes,type_ignores=[]),str(SOURCE),'exec'),namespace)
    return namespace


@pytest.mark.parametrize('model,expected',[
    ('A_+0','A_p0'),('A_-2','A_m2'),('A_-4','A_m4'),
    ('B_+1','B_p1'),('B_-1','B_m1'),('B_-3','B_m3'),
    ('C_+1','C_p1'),('C_-1','C_m1'),('C_-3','C_m3'),
    ('D_+0','D_p0'),('D_+2','D_p2'),('D_-2','D_m2'),
    ('E_+0','E_p0'),('E_+2','E_p2'),('E_-2','E_m2'),('E_-4','E_m4')])
def test_safe_model_names(model,expected):
    assert pure()['safe'](model)==expected


@pytest.mark.parametrize('dims,alpha,expected',[
    ([1,3,2],0,'T3_dS1_1_dS2_3_dF_2_alpha_p0.json'),
    ([2,2,1],1,'T3_dS1_2_dS2_2_dF_1_alpha_p1.json'),
    ([2,2,3],-3,'T3_dS1_2_dS2_2_dF_3_alpha_m3.json'),
    ([3,3,2],-4,'T3_dS1_3_dS2_3_dF_2_alpha_m4.json')])
def test_historical_filename(dims,alpha,expected):
    assert pure()['stem']({'dimensions':dims,'alpha':alpha})==expected


def test_inventory_zero_additional():
    obj={'models':{'E_-4':{'new_rgbeta_symbols':[],
          'preexisting_rgbeta_symbols':[],'total_additional_quartics':0}}}
    _,new,reused=pure()['inventory'](obj,'E_-4')
    assert new==() and reused==()


def test_inventory_valid_complex():
    obj={'models':{'B_+1':{'new_rgbeta_symbols':[{'name':'lambdaT3GeneralQ022000I1'}],
          'preexisting_rgbeta_symbols':[],'total_additional_quartics':1}}}
    _,new,reused=pure()['inventory'](obj,'B_+1')
    assert new==('lambdaT3GeneralQ022000I1',) and not reused


def test_inventory_inconsistent_fails():
    obj={'models':{'B_+1':{'new_rgbeta_symbols':[],
          'preexisting_rgbeta_symbols':[],'total_additional_quartics':1}}}
    with pytest.raises(ValueError):
        # The extracted inventory uses the source's fail() error helper.
        fn=pure()['inventory'];fn.__globals__['fail']=lambda msg: (_ for _ in ()).throw(ValueError(msg))
        fn(obj,'B_+1')


def test_row_of_failed_source_not_accepted():
    r=pure()['row_of']({'model':'A_+0','case':'small_uv_running_5e-9',
        'status':'FailedOrUnavailable','reason':'missing live RGBeta projection'})
    assert r['status']=='FailedOrUnavailable' and r['chi2'] is None
    assert r['reason']=='missing live RGBeta projection'


def test_row_of_has_all_five_observables():
    values=dict(zip(('delta_m21_sq_ev2','delta_m3l_sq_ev2','sin2_theta12',
             'sin2_theta13','sin2_theta23'),(7.53e-5,.00251,.3088,.02248,.47)))
    rec={'model':'B_+1','case':'large_uv_frozen_0p5','status':'CandidateConstrainedAccepted',
         'best_mass_safe':{'chi2':.001,'rhoF':.01,'scalar_masses_muS_tev':[1.,1.1],
                           'prediction':values}}
    row=pure()['row_of'](rec)
    assert row['chi2']==.001
    for k,v in values.items():assert row[k]==v


def test_no_silent_wolfram_fallback_and_per_model_failure_recording():
    text=SOURCE.read_text(encoding='utf-8')
    assert 'FAILED EXPORT' in text and 'FailedOrUnavailable' in text
    assert 'campaign.csv' in text and 'campaign.json' in text
    assert 'NoMassSafeNeutrinoFit' in text
    assert 'BetaTerm[Quartic,1]' in WOLFRAM.read_text(encoding='utf-8')

# Regression: historical E-class benchmarks are half of the non-E lambdaT3
# boundaries, while the folder/study labels remain the same.
@pytest.mark.parametrize('dimensions,nominal,expected', [
    ([1,3,2],5e-9,5e-9),([2,2,1],0.5,0.5),
    ([3,3,2],5e-9,2.5e-9),([3,3,2],0.5,0.25),
])
def test_original_e_class_lambdaT3_half_convention(dimensions,nominal,expected):
    import math
    tree=ast.parse(SOURCE.read_text(encoding='utf-8'))
    node=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='scenario_lambda_t3')
    scope={'fail':lambda msg: (_ for _ in ()).throw(ValueError(msg))}
    exec(compile(ast.Module(body=[node],type_ignores=[]),str(SOURCE),'exec'),scope)
    assert math.isclose(scope['scenario_lambda_t3']({'dimensions':dimensions},nominal),expected,rel_tol=0,abs_tol=1e-20)


def bc_functions(root):
    import json, math, hashlib
    tree=ast.parse(SOURCE.read_text(encoding='utf-8'))
    names={'validated_bc_fit','existing_fit'}
    nodes=[x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name in names]
    namespace={
        'ROOT':root,'REUSED':('B_-1','C_-1'),'math':math,'Path':Path,
        'read': lambda p: json.loads(Path(p).read_text(encoding='utf-8')),
        'sha':lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest(),
        'fail':lambda msg: (_ for _ in ()).throw(ValueError(msg))}
    exec(compile(ast.Module(body=nodes,type_ignores=[]),str(SOURCE),'exec'),namespace)
    return namespace


def save_validation_fixture(root,case='small_uv_frozen_5e-9',modified=False):
    import hashlib, json
    dest=root/'output/studies/general_scalar_extended'
    folder=dest/case/'b_small_fit_v2'
    folder.mkdir(parents=True)
    fit=folder/'B_m1_exploratory_fit.json'; fit.write_text('{"status":"ExploratoryAcceptedFit"}\n')
    cfg=folder/'B_m1_exploratory_fitted_config.json';cfg.write_text('{"config":"real"}\n')
    digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    lambda_uv=8.997e-6 if modified else 5e-9
    row={'case':case,'model':'B_-1','status':'SavedFitValidated',
         'replay_status':'CandidateReplayAgrees','chi2':0.001,'replayed_chi2':0.001,
         'quartic_coordinate_screen_pass':True,'yukawa_singular_screen_pass':True,
         'positivity_sampled':True,'scalar_light_muS_tev':1.1,
         'scalar_heavy_muS_tev':1.3,'rhoF':0.011,
         'max_yukawa_singular':1.2,'max_abs_quartic':0.33,
         'lambdaT3_uv':lambda_uv,
         'prediction':dict(zip(('delta_m21_sq_ev2','delta_m3l_sq_ev2','sin2_theta12',
                      'sin2_theta13','sin2_theta23'),(7.53e-5,0.00251,0.3088,0.02248,0.47))),
         # Simulate Windows path separators in validation report.
         'fit_path':str(fit.relative_to(root)).replace('/','\\'),
         'cfg_path':str(cfg.relative_to(root)).replace('/','\\'),
         'fit_sha256':digest(fit),'config_sha256':digest(cfg)}
    proof={'format':'t3-bc-six-fit-conditional-validation-v1','status':'ReplayValidated',
        'candidate_count':6,'perturbative_coordinate_screens_all_pass':True,
        'rg_replay_performed':True,'matching_physically_certified':False,
        'full_sequential_matching_performed':False,'results':[row]}
    place=dest/'validation/bc_six_validation_v1/report.json';place.parent.mkdir(parents=True)
    place.write_text(json.dumps(proof,indent=2)+'\n')
    return fit,cfg


def test_recovers_frozen_bc_fit_from_validated_nonstandard_filename(tmp_path):
    save_validation_fixture(tmp_path)
    result=bc_functions(tmp_path)['existing_fit']('small_uv_frozen_5e-9','B_-1')
    assert result['status']=='ExistingConditionalAccepted'
    assert result['best_mass_safe']['chi2']==0.001
    assert result['best_mass_safe']['prediction']['sin2_theta13']==0.02248
    assert result['original_boundary_relaxed'] is False
    assert result['source_fit'].endswith('B_m1_exploratory_fit.json')


def test_fails_closed_if_validation_source_changed(tmp_path):
    fit,_=save_validation_fixture(tmp_path)
    fit.write_text(fit.read_text()+'tampering!')
    with pytest.raises(ValueError,match='differs from validated SHA-256'):
        bc_functions(tmp_path)['existing_fit']('small_uv_frozen_5e-9','B_-1')


def test_marks_modified_bc_uv_fit_separately(tmp_path):
    save_validation_fixture(tmp_path,case='small_uv_running_5e-9',modified=True)
    result=bc_functions(tmp_path)['existing_fit']('small_uv_running_5e-9','B_-1')
    assert result['original_boundary_relaxed'] is True
    assert result['historical_lambdaT3_fixed']!=5e-9


def test_e_case_metadata_keeps_original_folder_label():
    source=SOURCE.read_text(encoding='utf-8')
    assert "meta.get('case') != ('0p5' if nominal_lam==.5 else '5e-9')" in source
    assert 'scenario_lambda_t3(entry,nominal_lam)' in source


def test_c_plus1_historical_seed_search_not_artificially_rejected():
    """C(+1) frozen small has a verified historical component |y|=2.746132..."""
    import math, numpy as np
    source=SOURCE.read_text(encoding='utf-8')
    tree=ast.parse(source)
    node=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='run_case')
    call=next(n for n in ast.walk(node) if isinstance(n,ast.Call) and
        isinstance(n.func,ast.Attribute) and n.func.attr=='search_box')
    bound=next(k.value for k in call.keywords if k.arg=='abs_yukawa_max')
    cap=eval(compile(ast.Expression(body=bound),str(SOURCE),'eval'),{'math':math})
    assert math.isclose(cap,math.sqrt(4*math.pi),rel_tol=1e-14)
    assert 2.746132227654911 < cap
    assert cap < 3.6  # still a perturbativity-derived search cap, not arbitrary 10


def test_current_errors_are_separate_from_append_only_history(tmp_path):
    import ast
    tree=ast.parse(SOURCE.read_text(encoding='utf-8'))
    run=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='run')
    inner=next(x for x in ast.walk(run) if isinstance(x,ast.FunctionDef) and x.name=='record_error')
    history=tmp_path/'errors.log'
    history.write_text('FIT E_+0: stale error\n',encoding='utf-8')
    recent=tmp_path/'errors_latest.log'
    recent.write_text('',encoding='utf-8')
    ns={'base':tmp_path,'fresh_errors':recent}
    exec(compile(ast.Module(body=[inner],type_ignores=[]),str(SOURCE),'exec'),ns)
    ns['record_error']('FIT','C_+1 small_uv_frozen_5e-9',ValueError('new error'))
    assert 'stale error' in history.read_text(encoding='utf-8')
    assert 'new error' in history.read_text(encoding='utf-8')
    assert 'stale error' not in recent.read_text(encoding='utf-8')
    assert 'new error' in recent.read_text(encoding='utf-8')
