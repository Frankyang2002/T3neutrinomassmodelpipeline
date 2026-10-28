"""Batch12 continuation offline regression tests; numerical physics runs on Windows."""
from __future__ import annotations

from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
import sys
from types import ModuleType, SimpleNamespace

import numpy as np
import pytest

from scripts import ContinueT3GeneralScalarBCExploratoryFit as c


def parameters(scale):
    return {k:float((i + 1) * scale / 18) for i,k in enumerate(c.NAMES)}


def write_case(root:Path, case:str, model:str, *,scale:float=1e-4):
    tag=c.old.b9.safe_tag(model)
    folder=root/'output/studies/general_scalar_extended'/case/'prior_v1'
    folder.mkdir(parents=True)
    v=parameters(scale)
    ordinary={}
    for k in ('y1','y2'):
        ordinary[k]={'real':[[v[f'{k}_{i}{j}'] for j in range(3)] for i in range(3)],
                     'imag':np.zeros((3,3)).tolist()}
    config={'base_state':{'ordinary':ordinary},'exploratory_general_scalar_fit':{
       'case':case,'model':model,'run_id':'prior_v1'}}
    hashes={}
    for i in range(5):
        relative=f'dummy/input{i}.txt'
        full=root/relative
        full.parent.mkdir(parents=True,exist_ok=True)
        full.write_text(f'input {i}',encoding='utf8')
        hashes[relative]=sha256(full.read_bytes()).hexdigest()
    record={'format':c.old.FORMAT,'case':case,'model':model,'run_id':'prior_v1',
       'status':'ExploratoryPoorFit','chi2_final':116.,
       'historical_fit_chi2':0.,'best_parameters':v,
       'input_extra_quartic_fraction_of_abs_lambdaT3':0.1,
       'input_mu2_fraction_of_min_scalar_m2':0.,
       'physical_matchete_rgbeta_map_certified':False,
       'full_sequential_general_scalar_matching_performed':False,
       'UV_running_enabled':c.old.b9.CASES[case][1],
       'source_manifest':{'files_sha256':hashes,
                          'fit_options':{'quartic_fraction':0.1,'mu2_fraction':0.}},
       'historical_fitted_config':'history/model.json'}
    (folder/f'{tag}_exploratory_fit.json').write_text(json.dumps(record),encoding='utf8')
    (folder/f'{tag}_exploratory_fitted_config.json').write_text(json.dumps(config),encoding='utf8')
    history=root/'history/model.json';history.parent.mkdir(parents=True)
    history.write_text('unmodified',encoding='utf8')
    return record,config


@pytest.mark.parametrize(('case','model'),[
 ('small_uv_frozen_5e-9','B_-1'),('large_uv_frozen_0p5','C_-1'),
 ('small_uv_running_5e-9','B_-1')])
def test_load_prior_and_real_yukawa_match(tmp_path,case,model):
    p,_=write_case(tmp_path,case,model)
    got,config,fp,cp=c.load_previous(case,model,'prior_v1',root=tmp_path)
    assert got['best_parameters']==p['best_parameters']
    assert config['base_state']['ordinary']['y2']['real'][2][2]==p['best_parameters']['y2_22']
    assert fp.is_file() and cp.is_file()
    c.check_previous_source_hashes(got,root=tmp_path)


def test_reject_prior_probe(tmp_path):
    p,_=write_case(tmp_path,'small_uv_frozen_5e-9','B_-1')
    path,_=c.previous_paths('small_uv_frozen_5e-9','B_-1','prior_v1',root=tmp_path)
    d=json.loads(path.read_text());d['status']='ExploratoryProbe'
    path.write_text(json.dumps(d))
    with pytest.raises(ValueError,match='not a completed'):
        c.load_previous('small_uv_frozen_5e-9','B_-1','prior_v1',root=tmp_path)


def test_reject_previous_yukawa_mismatch(tmp_path):
    write_case(tmp_path,'small_uv_frozen_5e-9','B_-1')
    _,cp=c.previous_paths('small_uv_frozen_5e-9','B_-1','prior_v1',root=tmp_path)
    d=json.loads(cp.read_text());d['base_state']['ordinary']['y1']['real'][0][0]+=1e-5
    cp.write_text(json.dumps(d))
    with pytest.raises(ValueError,match='does not exactly match'):
        c.load_previous('small_uv_frozen_5e-9','B_-1','prior_v1',root=tmp_path)


def test_reject_stale_hashes_and_traversal(tmp_path):
    p,_=write_case(tmp_path,'small_uv_frozen_5e-9','B_-1')
    (tmp_path/'dummy/input0.txt').write_text('mutated')
    with pytest.raises(ValueError,match='source changed'):
        c.check_previous_source_hashes(p,root=tmp_path)
    p['source_manifest']['files_sha256']={'../../etc/passwd':'a'*64}
    with pytest.raises(ValueError,match='missing original physics-source hash|escapes'):
        c.check_previous_source_hashes(p,root=tmp_path)


def test_previous_path_no_escape(tmp_path):
    for value in ['../other','a/b','bad id','']:
        with pytest.raises(ValueError):
            c.previous_paths('small_uv_frozen_5e-9','B_-1',value,root=tmp_path)


def fake_benchmark_module(monkeypatch):
    pkg=ModuleType('Numerical');fit=ModuleType('Numerical.fitting')
    mod=ModuleType('Numerical.fitting.SobolBenchmarkSearch')
    class Parameter:
        def __init__(self,name,binding,low,high):
            self.name,self.binding,self.low,self.high=name,binding,low,high
    mod.BenchmarkParameter=Parameter
    mod.generate_sobol_parameter_points=lambda definitions,n_points,seed: []
    for name,value in [('Numerical',pkg),('Numerical.fitting',fit),
                       ('Numerical.fitting.SobolBenchmarkSearch',mod)]:
        monkeypatch.setitem(sys.modules,name,value)
    return mod


@pytest.mark.parametrize('scale',[2.0,3.1e-4])
def test_adaptive_bounds_follow_yukawa_magnitude(monkeypatch,scale):
    fake_benchmark_module(monkeypatch)
    seed=parameters(scale)
    out=c.adaptive_definitions(seed,bound=3,relative_window=.5,
                              floor_ratio=.02,absolute_floor=1e-8)
    assert len(out)==18 and len({x.binding for x in out})==18
    for x in out:
        assert x.low<seed[x.name]<x.high
        assert -3<=x.low and x.high<=3
    assert out[0].high-out[0].low<3*scale
    if scale<0.01:
        assert max(x.high-x.low for x in out)<0.001


def test_invalid_adaptive_parameters(monkeypatch):
    fake_benchmark_module(monkeypatch)
    seed=parameters(0.1)
    with pytest.raises(ValueError):c.adaptive_definitions(seed,bound=3,relative_window=0,floor_ratio=.02,absolute_floor=1e-8)
    with pytest.raises(ValueError):c.adaptive_definitions(seed,bound=.001,relative_window=.5,floor_ratio=.02,absolute_floor=1e-8)


def test_candidate_collection_seed_source(monkeypatch):
    fake_benchmark_module(monkeypatch)
    class Evaluator:
        def __call__(self,point):
            return SimpleNamespace(chi2=100.,ordering='NO',
              observable_names=c.FIELDS,prediction=np.ones(5),
              whitened_residual=np.ones(5))
    seed=parameters(0.1)
    out,failed=c.evaluate_candidates(Evaluator(),seed,(),sobol_points=0,sobol_seed=5)
    assert len(out)==1 and failed==[] and out[0]['source']=='saved_previous_best_recomputed'


def test_result_markdown_uncertified():
    value={'case':'small_uv_frozen_5e-9','model':'B_-1','run_id':'v2',
         'previous_run_id':'prior_v1','chi2_previous':120.,'chi2_start_recomputed':120.,
         'chi2_final':100.,'passes_approx_chi2_5p86':False,
         'prediction':dict.fromkeys(c.FIELDS,0.5),
         'candidate_evaluation_count':20,'local_sobol_failure_count':0}
    s=c.result_markdown(value)
    assert '18 real Yukawa' in s and 'UNCERTIFIED' not in s.upper() or 'NOT a certified' in s
    assert 'Prior fit chi2: 120' in s and 'Continued best chi2: 100' in s


def test_new_only_destination(tmp_path):
    dest=c.old.fit_destination('small_uv_frozen_5e-9','v2',root=tmp_path)
    dest.mkdir(parents=True)
    with pytest.raises(ValueError,match='already exists'):
        c.old.fit_destination('small_uv_frozen_5e-9','v2',root=tmp_path)


def test_available_uploaded_results_schema():
    """Optional local verification; never needed in user's repository pytest."""
    root=Path('/mnt/data/general_scalar_extended.zip')
    if not root.exists():pytest.skip('local analysis ZIP unavailable')
    import zipfile
    with zipfile.ZipFile(root) as z:
        for case,run,tag in [('small_uv_frozen_5e-9','b_small_fit_v1','B_m1'),
                             ('large_uv_frozen_0p5','c_large_fit_v1','C_m1')]:
            p=f'general_scalar_extended/{case}/{run}/'
            d=json.loads(z.read(p+f'{tag}_exploratory_fit.json'))
            assert d['format']==c.old.FORMAT
            assert len(c.values_from_record(d))==18


def test_full_continuation_wiring_with_fake_physics(monkeypatch,tmp_path):
    case,model='small_uv_frozen_5e-9','B_-1'
    prior,cfg=write_case(tmp_path,case,model,scale=.3)
    fake_benchmark_module(monkeypatch)
    def fake_refine(**kwargs):
        assert kwargs['active_counts']==(18,)
        assert kwargs['max_nfev']==(25,)
        assert kwargs['seed_count']==1
        assert len(kwargs['definitions'])==18
        assert kwargs['candidates'][0]['source']=='saved_previous_best_recomputed'
        p=deepcopy(kwargs['candidates'][0]['parameters'])
        p['y1_00']+=.001
        return {'parameters':p,'chi2':5.,'source':'all18_stage'}, {'enabled':True,'stages':[{'stage':'all18'}]}
    monkeypatch.setattr(sys.modules['Numerical.fitting.SobolBenchmarkSearch'],
                        'refine_sobol_candidates',fake_refine,raising=False)
    registry=ModuleType('Numerical.experimental.GeneralT3ScalarRGState')
    registry.load_registry=lambda path: {}
    exp=ModuleType('Numerical.experimental')
    monkeypatch.setitem(sys.modules,'Numerical.experimental',exp)
    monkeypatch.setitem(sys.modules,'Numerical.experimental.GeneralT3ScalarRGState',registry)
    observables=ModuleType('Numerical.fitting.OscillationFit')
    observables.OscillationFitTarget=SimpleNamespace(from_json=lambda path:'target')
    monkeypatch.setitem(sys.modules,'Numerical.fitting.OscillationFit',observables)
    base=ModuleType('Numerical.orchestration.RunSmallLambdaRefit')
    base.TARGET=tmp_path/'dummy/input4.txt'
    base.setmat=lambda d,z:d.update(real=np.asarray(z.real).tolist(),imag=np.asarray(z.imag).tolist())
    orch=ModuleType('Numerical.orchestration');orch.RunSmallLambdaRefit=base
    monkeypatch.setitem(sys.modules,'Numerical.orchestration',orch)
    monkeypatch.setitem(sys.modules,'Numerical.orchestration.RunSmallLambdaRefit',base)
    coeff=ModuleType('scripts.RunGeneralT3NumericalC5')
    coeff.DEFAULT_SOURCE=tmp_path/'fake_coeff'
    monkeypatch.setitem(sys.modules,'scripts.RunGeneralT3NumericalC5',coeff)
    monkeypatch.setattr(c.old.b9,'study_source',lambda *a,**k:{
        'cfg_path':tmp_path/'history/model.json',
        'uv_on':False,
        'cfg':cfg})
    monkeypatch.setattr(c.old.b10,'check_convention_evidence',lambda *a,**k: True)
    monkeypatch.setattr(c.old.b7,'read_json',lambda path:{})
    monkeypatch.setattr(c.old,'compile_source',lambda path:'compiled')
    monkeypatch.setattr(c.old,'beta_bundle',lambda *a,**k:{})
    class FakeEvaluator:
        def __init__(self,*a,**kw):
            assert kw['quartic_fraction']==.1 and kw['mu2_fraction']==0.
            self.calls=0
        def evaluate(self,point,*,detail=False):
            self.calls+=1
            chi=5. if point['y1_00']!=prior['best_parameters']['y1_00'] else 116.
            fit=SimpleNamespace(chi2=chi,ordering='NO',
                observable_names=c.FIELDS,prediction=np.ones(5)*.3,
                whitened_residual=np.full(5,np.sqrt(chi/5)))
            if detail:return fit,{'trajectory':'fake'}
            return fit
        def __call__(self,point):return self.evaluate(point)
    monkeypatch.setattr(c.old,'CandidateEvaluator',FakeEvaluator)
    r=c.continue_run(case,model,'prior_v1','v2',root=tmp_path,max_nfev=25)
    assert r['chi2_final']==5. and r['passes_approx_chi2_5p86']
    assert r['previous_run_id']=='prior_v1' and r['optimiser']['enabled']
    dest=tmp_path/'output/studies/general_scalar_extended'/case/'v2'
    assert (dest/'summary.md').is_file()
    newcfg=json.loads((dest/'B_m1_exploratory_fitted_config.json').read_text())
    assert newcfg['exploratory_general_scalar_fit']['continued_from']=='prior_v1'
    assert prior['best_parameters']['y1_00']!=r['best_parameters']['y1_00']
    # Historical reference and original run remain unchanged.
    assert (tmp_path/'history/model.json').read_text()=='unmodified'
    assert json.loads(c.previous_paths(case,model,'prior_v1',root=tmp_path)[0].read_text())==prior


def test_continuation_rejects_existing_destination_before_work(tmp_path):
    case,model='large_uv_frozen_0p5','C_-1'
    write_case(tmp_path,case,model)
    out=c.old.fit_destination(case,'v2',root=tmp_path)
    out.mkdir(parents=True)
    with pytest.raises(ValueError,match='already exists'):
        c.continue_run(case,model,'prior_v1','v2',root=tmp_path)
