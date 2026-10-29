"""Stage20 general-quartic source and registry tests; no Wolfram required."""
from pathlib import Path
import importlib.util
import json
import re
import sys
from copy import deepcopy

import pytest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from scripts import GenerateGeneralT3RGBetaBasis as gen
from scripts import RunGeneralT3GeneralQuarticRGBeta as runner

@pytest.fixture(scope='module')
def manifest():
    return json.loads(gen.MANIFEST.read_text(encoding='utf-8'))

@pytest.fixture(scope='module')
def registry(tmp_path_factory,manifest):
    return gen.generate(manifest,tmp_path_factory.mktemp('stage20'))


def test_16_model_manifest_and_50_quartics(manifest):
    assert manifest['model_count']==16
    assert len(manifest['models'])==16
    assert sum(v['extra_quartic_complex_couplings'] for v in manifest['models'].values())==50
    assert sum(v['quadratic_complex_couplings_not_enabled'] for v in manifest['models'].values())==8


def test_all_16_generators_close_rank_exactly(manifest):
    for name,model in manifest['models'].items():
        result=gen.make_model(name,model)
        assert result['total_additional_quartics']==model['extra_quartic_complex_couplings']
        assert len(result['new_rgbeta_symbols'])+len(result['preexisting_rgbeta_symbols'])==result['total_additional_quartics']
        assert result['quadratics_not_implemented']==model['quadratic_complex_couplings_not_enabled']


def test_aggregate_closure_and_no_double_count(registry):
    assert registry['totals']==dict(matchete_extra_quartics=50,
          rgbeta_preexisting_structural=12,rgbeta_new_vertices=38,quadratic_couplings_blocked=8)
    assert len(registry['models'])==16
    for name in ['B_-1','C_-1']:
        point=registry['models'][name]
        assert len(point['preexisting_rgbeta_symbols'])==6
        assert len(point['new_rgbeta_symbols'])==1
        assert point['new_rgbeta_symbols'][0]['counts']==[1,1,0,1,0,1]
    for name,point in registry['models'].items():
        if name not in ['B_-1','C_-1']:
            assert not point['preexisting_rgbeta_symbols']


def test_every_new_vertex_has_exactly_four_fields_and_native_invariant(registry):
    for body in registry['models'].values():
        assert len({o['name'] for o in body['new_rgbeta_symbols']})==len(body['new_rgbeta_symbols'])
        for v in body['new_rgbeta_symbols']:
            assert len(v['fields'])==4
            assert sum(v['counts'])==4
            assert v['invariant']
            assert v['mapping_status']=='NotNormalised'
            assert 'lambdaT3Extra' not in v['name']  # do not falsely identify with Matchete
            assert 'T3RGBetaS2TripleInvariant' in v['invariant'] or any(
                name in v['invariant'] for name in ('eps[','del[','delS2[')) or v['invariant']=='1'


def test_generated_wolfram_files_have_no_explicit_production_edits(tmp_path,manifest):
    registry=gen.generate(manifest,tmp_path)
    assert len(list(tmp_path.glob('*_extra_quartics.wl')))==16
    for model,info in registry['models'].items():
        code=(tmp_path/f'{model}_extra_quartics.wl').read_text(encoding='utf-8')
        assert code.count('AddQuartic[')==len(info['new_rgbeta_symbols'])
        assert code.count('SelfConjugate -> False')==len(info['new_rgbeta_symbols'])
        assert 'ResetModel[' not in code
        assert 'Abort[]' in code


def test_doublet_bose_sym_eliminates_antisymmetric_pair():
    assert gen.independent_weights(['Hdag','Hdag','S1','S1'],2)
    assert len(gen.independent_weights(['Hdag','Hdag','S1','S1'],2))==1
    assert len(gen.independent_weights(['H','Hdag','S1dag','S2dag'],2))==2
    assert len(gen.independent_weights(['S1dag']*4,2))==0


def test_triplet_bose_symmetry_yields_one_or_two_singlets():
    assert len(gen.independent_weights(['S2dag']*4,3))==1
    assert len(gen.independent_weights(['S1','S1dag','S2dag','S2dag'],3))==2


def test_quartics_do_not_change_original_catalogue(manifest,tmp_path):
    previous=deepcopy(manifest)
    gen.generate(manifest,tmp_path)
    assert previous==manifest


def test_registry_refuses_missing_model(manifest,tmp_path):
    bad=deepcopy(manifest)
    bad['models'].pop('B_-1')
    with pytest.raises(ValueError,match='16-model'):
        gen.generate(bad,tmp_path)


def test_existing_source_check_fails_closed(registry,tmp_path):
    missing=runner.verify_existing_rgbeta(registry,tmp_path/'missing.wl')
    assert missing['status']=='BlockedMissingProductionRGBetaSource'
    (tmp_path/'partial.wl').write_text('AddQuartic[lambdaHHdagS2S2, ...]')
    with pytest.raises(ValueError,match='not present'):
        runner.verify_existing_rgbeta(registry,tmp_path/'partial.wl')


def test_export_checker_rejects_missing_beta(tmp_path):
    p=tmp_path/'export.json'
    p.write_text(json.dumps({'status':'Success','metadata':{'GeneralScalarModel':'A_-2','GeneralScalarStage':'UV',
                'GeneralScalarNewBetaKeys':['lambdaX']},'report_betas':{'g2':'0'}}))
    with pytest.raises(ValueError,match='missing'):
        runner.check_export(p,'A_-2','UV',['lambdaX'])


def test_export_checker_accepts_complete_candidate_beta(tmp_path):
    p=tmp_path/'export.json'
    p.write_text(json.dumps({'status':'Success','metadata':{'GeneralScalarModel':'A_-2','GeneralScalarStage':'UV',
                'GeneralScalarNewBetaKeys':['lambdaX']},'report_betas':{'g2':'0','lambdaX':'lambdaX^2'}}))
    assert runner.check_export(p,'A_-2','UV',['lambdaX'])['new_complex_quartics']==1


def test_generated_b_c_missing_mixed_basis_is_not_same_contraction(registry):
    b=registry['models']['B_-1']['new_rgbeta_symbols'][0]['invariant']
    c=registry['models']['C_-1']['new_rgbeta_symbols'][0]['invariant']
    assert b==c
    assert 'del[SU2L @ fund, #4, #1]' in b
    assert 'eps[SU2L @ fund, #2, #3]' in b
    assert 'eps[SU2L @ fund, #1, #4]' not in b


def test_registry_marks_all_quadratic_mixing_unimplemented(registry):
    assert registry['totals']['quadratic_couplings_blocked']==8
    assert any(v['quadratics_not_implemented'] for v in registry['models'].values())


def test_ordinary_general_numerical_module_optional_production_import():
    try:
        spec = importlib.util.find_spec('Numerical.core.State')
    except ModuleNotFoundError:
        spec = None
    if spec is None:
        pytest.skip('Real production Numerical/core/State.py is local to user repository')
    from Numerical.experimental.GeneralT3ScalarRGState import _model_name
    from Numerical.core.State import T3Representation
    assert _model_name(T3Representation(d_s1=2,d_s2=2,d_f=1,alpha=-1))=='B_-1'
    assert _model_name(T3Representation(d_s1=3,d_s2=3,d_f=2,alpha=-2))=='E_-2'


def test_generic_B_state_roundtrip_and_synthetic_complete_rge(registry):
    try:
        spec = importlib.util.find_spec('Numerical.core.State')
    except ModuleNotFoundError:
        spec = None
    if spec is None:
        pytest.skip('Production numerical modules are only present in the local full repository.')
    from Numerical.experimental.GeneralT3CoupledUV import construct_uv_state
    from Numerical.experimental.GeneralT3ScalarRGState import (
        GeneralT3ScalarRGState,pack_general_scalar_state,unpack_general_scalar_state,
        evaluate_general_scalar_derivative,run_general_scalar_segment,
        project_general_after_fermion,
    )
    from Numerical.core.StateVector import StateVectorLayout
    config=json.loads((ROOT/'docs/general_t3_threshold_benchmark.json').read_text())
    baseline=construct_uv_state('B_m1',config)
    name=registry['models']['B_-1']['new_rgbeta_symbols'][0]['name']
    initial=GeneralT3ScalarRGState(baseline,{name:3e-3+2e-3j})
    v,layout=pack_general_scalar_state(initial,registry)
    assert isinstance(layout,StateVectorLayout)
    assert layout.names[-1]==name
    restored=unpack_general_scalar_state(v,layout,registry,uv=True,mu_gev=baseline.mu_gev)
    assert restored.extra[name]==initial.extra[name]
    eft=project_general_after_fermion(initial,registry)
    assert eft.extra==initial.extra
    assert eft.baseline.mu_gev==baseline.mu_gev
    expressions={block.name:(f'0*{block.name}' if block.kind=='complex_matrix' else '0')
                 for block in layout.blocks}
    meta={'dS1':2,'dS2':2,'dF':1,'alpha':-1,
          'GeneralScalarModel':'B_-1','GeneralScalarStage':'UV',
          'GeneralScalarBasisStatus':'CandidateNativeRGBetaUnmappedToMatchete',
          'GeneralScalarNewBetaKeys':[name],
          'ReportBetaConvention':'16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g'}
    payload={'status':'Success','metadata':meta,'report_betas':expressions}
    dy,dlay=evaluate_general_scalar_derivative(payload,initial,registry)
    assert dlay==layout
    assert not any(dy)
    evolved=run_general_scalar_segment(initial,payload,registry,9500)
    assert evolved['solver_success'] is True
    assert evolved['state'].extra[name]==initial.extra[name]
    assert evolved['physical_neutrino_mass_prediction'] is False


def test_model_selection_accepts_powershell_comma_separated_names():
    assert runner.normalise_model_selection(['B_-1,C_-1'])==['B_-1','C_-1']
    assert runner.normalise_model_selection(['B_-1', 'C_-1'])==['B_-1','C_-1']
    assert runner.normalise_model_selection([' B_-1, C_-1 ', 'B_-1'])==['B_-1','C_-1']


def test_model_selection_empty_means_all_and_malformed_fails():
    assert runner.normalise_model_selection(None) is None
    assert runner.normalise_model_selection([]) is None
    with pytest.raises(ValueError, match='Empty model name'):
        runner.normalise_model_selection(['B_-1,'])


def test_model_selection_unknown_names_are_still_rejected():
    selection=runner.normalise_model_selection(['B_-1,not_a_model'])
    assert sorted(set(selection)-{'B_-1','C_-1'})==['not_a_model']


def test_stage20_run_accepts_two_models_in_one_argument(monkeypatch, tmp_path):
    from scripts import GenerateGeneralT3RGBetaBasis as gen
    monkeypatch.setattr(runner, 'ROOT', tmp_path)
    monkeypatch.setattr(runner, 'verify_existing_rgbeta',
        lambda registry, source: {'status': 'SourceStructuralReuseVerified'})
    result=runner.run_stage20(
        models=['B_-1,C_-1'], output_dir=tmp_path/'generated'
    )
    assert result['selected_models']==['B_-1','C_-1']
    assert result['status']=='CandidateBasisGenerated'
    assert result['reports']==[]
    assert (tmp_path/'output/general_t3_stage20/stage20_summary.json').exists()


# Regressions for Windows WolframScript -file argument handling.
def test_wolframscript_negative_alpha_uses_production_token():
    assert runner.wolfram_integer_token(-1) == 'm1'
    assert runner.wolfram_integer_token(-3) == 'm3'
    assert runner.wolfram_integer_token(0) == '0'
    assert runner.wolfram_integer_token(2) == '2'


def test_generated_wolfram_command_never_passes_negative_flag(tmp_path):
    args = runner.general_rgbeta_command(
        'wolframscript',tmp_path/'RGE/running/rgbeta/RunGeneralT3RGBeta.wl',
        'UV', 2, 2, 1, -1, 'B_-1', tmp_path/'b.json',
    )
    assert args[0] == 'wolframscript'
    assert args[1] == '-file'
    assert args[3:9] == ['UV', '2', '2', '1', 'm1', 'B_-1']
    assert '-1' not in args
    with pytest.raises(ValueError, match='Unsupported general RGBeta stage'):
        runner.general_rgbeta_command('wolframscript',tmp_path/'test.wl','bad',2,2,1,-1,'B_-1',tmp_path/'b.json')


def test_wolfram_failure_includes_command_and_log_tail(monkeypatch,tmp_path):
    import types
    monkeypatch.setattr(runner,'ROOT',tmp_path)
    monkeypatch.setattr(runner,'verify_existing_rgbeta',
        lambda registry, source: {'status':'SourceStructuralReuseVerified'})
    monkeypatch.setattr(runner.shutil,'which',lambda executable: 'wolframscript')
    calls=[]
    def fake_run(command,**kwargs):
        calls.append(command)
        return types.SimpleNamespace(returncode=2,stdout='Usage: wolframscript -file ...',stderr='unknown option')
    monkeypatch.setattr(runner.subprocess,'run',fake_run)
    with pytest.raises(RuntimeError,match='unknown option'):
        runner.run_stage20(export_betas=True,models=['B_-1'],stage='UV',
                           output_dir=tmp_path/'generated')
    assert calls[0][7]=='m1'
    p=tmp_path/'output/general_t3_stage20/B_-1_uv_rgbeta.log'
    log=p.read_text(encoding='utf-8')
    assert 'Command:' in log
    assert 'm1' in log
    assert 'Usage: wolframscript' in log
    assert 'unknown option' in log


def test_no_generated_rgbeta_invariance_option_in_all_16_models(registry, tmp_path):
    gen.generate(json.loads(gen.MANIFEST.read_text(encoding='utf-8')),tmp_path)
    for model in registry['models']:
        text=(tmp_path/f'{model}_extra_quartics.wl').read_text(encoding='utf-8')
        assert 'CheckInvariance ->' not in text
        assert text.count('AddQuartic[')==len(registry['models'][model]['new_rgbeta_symbols'])


def test_independent_gauge_charge_parity_and_rank_audit(registry):
    result=gen.audit_candidate_gauge_structures(registry)
    assert result['candidate_quartics_checked']==38
    assert result['status']=='AnalyticChargeParityFieldTypingAndInvariantRankVerified'
    assert result['native_rgbeta_checkinvariance_used'] is False
    assert result['native_rgbeta_su2_ward_identity_verified'] is False


def test_charge_check_catches_invalid_candidate(registry):
    broken=deepcopy(registry)
    broken['models']['B_-1']['new_rgbeta_symbols'][0]['counts']=[2,0,0,1,0,1]
    with pytest.raises(ValueError,match='hypercharge'):
        gen.audit_candidate_gauge_structures(broken)


def test_su2_b_c_mixed_delta_epsilon_orientation(registry):
    # Diagonal SU(2) transformation: fundamental weights (+1,-1), dual (-1,+1).
    # delta(H, S2dag) and epsilon(Hdag, S1dag) are invariant.
    # epsilon(H, S2dag) is *not*; its entries carry nonzero weight.
    from itertools import product
    eps=lambda a,b: { (0,1):1,(1,0):-1 }.get((a,b),0)
    for model in ('B_-1','C_-1'):
        inv=registry['models'][model]['new_rgbeta_symbols'][0]['invariant']
        assert inv=='(del[SU2L @ fund, #4, #1] * eps[SU2L @ fund, #2, #3])'
        for a,b,c,d in product(range(2),repeat=4):
            tensor=(a==d)*eps(b,c)
            weight=(1-2*a)-(1-2*b)-(1-2*c)-(1-2*d)
            assert not tensor or weight==0
    assert (1-2*0)-(1-2*0)-(1-2*1)-(1-2*1)!=0  # old mixed epsilon pattern


def test_runner_check_gauge_is_honest_about_native_wolfram(monkeypatch,tmp_path):
    monkeypatch.setattr(runner,'ROOT',tmp_path)
    monkeypatch.setattr(runner,'verify_existing_rgbeta',lambda registry,source:{'status':'SourceStructuralReuseVerified'})
    monkeypatch.setattr(runner.shutil,'which',lambda exe:'wolframscript')
    def fake_run(cmd,**kwargs):
        assert 'T3_GENERAL_RGBETA_VERIFY_GAUGE' not in kwargs['env']
        dest=Path(cmd[-1]); dest.parent.mkdir(parents=True,exist_ok=True)
        model,stage=cmd[-2],cmd[3]
        manifest=json.loads(gen.MANIFEST.read_text())
        names=[q['name'] for q in gen.make_model(model,manifest['models'][model])['new_rgbeta_symbols']]
        dest.write_text(json.dumps({'status':'Success','metadata':{
            'GeneralScalarModel':model,'GeneralScalarStage':stage,'GeneralScalarNewBetaKeys':names},
            'report_betas':{name:'0' for name in names}}))
        import types
        return types.SimpleNamespace(returncode=0,stdout='OK',stderr='')
    monkeypatch.setattr(runner.subprocess,'run',fake_run)
    result=runner.run_stage20(export_betas=True,check_gauge=True,models=['B_-1'],stage='UV',output_dir=tmp_path/'extra')
    assert result['status']=='RGBetaCandidateBetasExportedForSelectedModels'
    assert result['candidate_gauge_structure_audit']['candidate_quartics_checked']==38
    assert result['native_rgbeta_checkinvariance_used'] is False
    assert result['gauge_invariance_checked'] is False


def test_b_c_corrected_su2_contraction_under_non_diagonal_rotation(registry):
    import numpy as np
    eps=np.array([[0,1],[-1,0]],dtype=complex)
    delta=np.eye(2,dtype=complex)
    right=np.einsum('da,bc->abcd',delta,eps)
    old=np.einsum('ad,bc->abcd',eps,eps)
    theta=0.37
    sx=np.array([[0,1],[1,0]],dtype=complex)
    u=np.cos(theta)*delta+1j*np.sin(theta)*sx
    ub=np.conjugate(u)
    def rotate(tensor):
        return np.einsum('abcd,ae,bf,cg,dh->efgh',tensor,u,ub,ub,ub)
    assert np.max(np.abs(rotate(right)-right))<1e-12
    assert np.max(np.abs(rotate(old)-old))>1e-2
    for model in ('B_-1','C_-1'):
        assert registry['models'][model]['new_rgbeta_symbols'][0]['invariant']==(
            '(del[SU2L @ fund, #4, #1] * eps[SU2L @ fund, #2, #3])'
        )


# Regression for multi-invariant quartic mixing. B/C(-1) share the field
# content of lambdaHHdagS1barS2barCross with a new, independent quartic.
def test_stage20_uses_full_rgbeta_quartic_projection():
    wl=(ROOT/'RGE/running/rgbeta/RunGeneralT3RGBeta.wl').read_text(encoding='utf-8')
    assert 'BetaTerm[Quartic, 1]' in wl
    assert 'BetaTerm[coupling,1]' not in wl
    assert 'BetaTerm[coupling, 1]' not in wl
    assert 'CheckPotentialCouplingMixing::mixes' in wl
    assert 'Join[baseBetas, KeyTake[quarticNames, baseQuarticKeys]]' in wl
    assert 'extraBetas = KeyTake[quarticNames, newNames]' in wl
    assert 'GeneralScalarQuarticBetaBasis' in wl


def test_stage20_quartic_projection_guards_missing_or_failed_betas():
    wl=(ROOT/'RGE/running/rgbeta/RunGeneralT3RGBeta.wl').read_text(encoding='utf-8')
    assert 'missingNames = Complement[newNames, Keys[quarticNames]]' in wl
    assert 'If[!AssociationQ[projectedQuartics]' in wl
    assert 'AnyTrue[Values[quarticNames]' in wl
    assert 'Length[extraBetas] =!= Length[newNames]' in wl


def test_stage20_existing_bc_extra_mixes_same_fields_in_rgbeta_registry(registry):
    for model in ('B_-1', 'C_-1'):
        rows=registry['models'][model]
        new=rows['new_rgbeta_symbols'][0]
        existing=[x for x in rows['preexisting_rgbeta_symbols']
                  if x['counts']==new['counts']]
        assert [x['existing_name'] for x in existing]==['lambdaHHdagS1barS2barCross']
