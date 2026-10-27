"""Offline contract tests; these do not claim live Matchete EFT1 success."""
from __future__ import annotations
import json
from pathlib import Path
from types import SimpleNamespace
import pytest

from scripts import RunGeneralT3FermionFirstEFT1Tree as task


def certified(model='B_-1', stem='B_m1_final_eft1_fermion_first_tree'):
    return dict(
        format=task.FORMAT, status=task.STATUS, model=model,
        matching_approach='FermionFirstTreeOnlyScalarsRemainLight',
        matching_EFT_order=5, matching_loop_order=0,
        mass_basis_ordinary_check_passed=True,
        fermion_heavy_mass_basis_detailed_check_passed=True,
        source_quartic_count=17, source_yukawa_count=2,
        scalar_light_fields=['NewScalarX','NewScalarY'],
        integrated_heavy_fields=['NewFermion'], tree_matching_completed=True,
        fermion_field_absent=True, x_scalar_present=True,
        y_scalar_present=True, llss_contains_x=True,llss_contains_y=True,
        llss_term_count=4, eft1_rge_performed=False,
        scalar_threshold_matching_performed=False,
        one_loop_hard_fermion_matching_performed=False,
        physical_C5_computed=False, sequential_EFT1_completed=False,
        physical_matching_performed=False,
        matchete_rgbeta_quartic_normalization_certified=False,
        matched_eft_file=stem+'_eft1_tree.inputform.txt',
        transition_representation='EOMEFT',
        transition_eom_eft1_file=stem+'_eft1_transition_eom.inputform.txt',
        llss_inputform_file=stem+'_llss_tree.inputform.txt',
    )


def write_artifacts(directory, stem):
    (directory/(stem+'_eft1_tree.inputform.txt')).write_text('Field[NewScalarX] + Field[NewScalarY]')
    (directory/(stem+'_llss_tree.inputform.txt')).write_text('GammaCC Field[NewScalarX] Field[NewScalarY]')
    (directory/(stem+'_eft1_transition_eom.inputform.txt')).write_text('GammaCC Field[NewScalarX] Field[NewScalarY]')


@pytest.mark.parametrize('model',['B_-1','C_-1'])
def test_stem(model):
    assert task._stem(model,'final')==model.replace('-','m')+'_final_eft1_fermion_first_tree'

@pytest.mark.parametrize('model,endpoint',[('B_+1','final'),('B_-1','banana')])
def test_reject_unsupported(model,endpoint):
    with pytest.raises(ValueError,match='Unsupported'):
        task._stem(model,endpoint)


def test_certificate_and_exports(tmp_path):
    stem='B_m1_final_eft1_fermion_first_tree'
    write_artifacts(tmp_path,stem)
    a,b=task.validate_tree_report(certified(), 'B_-1',tmp_path,stem)
    assert a.name.endswith('_eft1_tree.inputform.txt')
    assert b.name.endswith('_llss_tree.inputform.txt')


@pytest.mark.parametrize('field,wrong',[
    ('status','OneLoopDirectSymbolicWeinbergProbePassedNumericalC5Pending'),
    ('format','wrong'),('model','C_-1'),('matching_EFT_order',6),
    ('matching_loop_order',1),('mass_basis_ordinary_check_passed',False),
    ('fermion_heavy_mass_basis_detailed_check_passed',False),
    ('source_quartic_count',16),('scalar_light_fields',['NewScalarX']),
    ('integrated_heavy_fields',['NewScalarX','NewScalarY','NewFermion']),
    ('tree_matching_completed',False),('fermion_field_absent',False),
    ('x_scalar_present',False),('y_scalar_present',False),
    ('llss_contains_x',False),('llss_contains_y',False),
    ('llss_term_count',0),('eft1_rge_performed',True),
    ('scalar_threshold_matching_performed',True),
    ('one_loop_hard_fermion_matching_performed',True),
    ('physical_C5_computed',True),('sequential_EFT1_completed',True),
    ('physical_matching_performed',True),
    ('matchete_rgbeta_quartic_normalization_certified',True),
])
def test_reject_invalid_result(tmp_path,field,wrong):
    stem='B_m1_final_eft1_fermion_first_tree'
    write_artifacts(tmp_path,stem)
    report=certified()
    report[field]=wrong
    with pytest.raises(ValueError):
        task.validate_tree_report(report,'B_-1',tmp_path,stem)


def test_reject_invalid_artifact(tmp_path):
    stem='B_m1_final_eft1_fermion_first_tree'
    report=certified()
    with pytest.raises(ValueError,match='missing or empty'):
        task.validate_tree_report(report,'B_-1',tmp_path,stem)
    write_artifacts(tmp_path,stem)
    (tmp_path/(stem+'_llss_tree.inputform.txt')).write_text('No scalars')
    with pytest.raises(ValueError,match='lacks certified'):
        task.validate_tree_report(report,'B_-1',tmp_path,stem)
    (tmp_path/(stem+'_llss_tree.inputform.txt')).write_text('GammaCC Field[NewScalarX] Field[NewScalarY] Field[NewFermion]')
    with pytest.raises(ValueError,match='integrated fermion'):
        task.validate_tree_report(report,'B_-1',tmp_path,stem)


def test_wolfram_source_f_first_and_native_invariants():
    s=(Path(__file__).resolve().parents[2]/'Lagrangian/RunGeneralT3FermionFirstEFT1Tree.wl').read_text()
    assert 'Mass->{Light,MSX}' in s and 'Mass->{Light,MSY}' in s
    assert 'Mass->{Heavy,MSX}' not in s and 'Mass->{Heavy,MSY}' not in s
    assert 'RunT3Matching[lag,5,0]' in s
    assert 'matchingResult["EOMEFT"]' in s
    assert 'transition_eom_eft1_file' in s
    assert 'RunT3Matching[lag,5,1]' not in s
    assert 'T3MECheckpoint["fermion-first-LLSS-sector-extracted"]' in s
    for f in ('T3MEBuildNativeTopology','T3MEBuildNativeH2Portals',
              'T3MEBuildNativeS2ScalarQuartics','T3MEBuildNativeYukawa1',
              'T3MEBuildNativeYukawa2','T3MEBuildNativeExtraQuartics'):
        assert f in s
    assert '"sequential_EFT1_completed"->False' in s


def test_missing_valid_prior_certificate_fail_closed(monkeypatch,tmp_path):
    monkeypatch.setattr(task,'prepare_model',lambda *a,**k:{'payload':{'unitary_U':{}},'parameters':tmp_path/'parameters.json'})
    monkeypatch.setattr(task,'_prior_certificate',lambda *a,**k: (_ for _ in ()).throw(ValueError('prior check absent')))
    with pytest.raises(ValueError,match='prior check absent'):
        task.run_model('B_-1',outdir=tmp_path/'result')


def test_python_runner_mock_subprocess(tmp_path,monkeypatch):
    root=tmp_path/'project'
    (root/'Lagrangian').mkdir(parents=True)
    (root/'output').mkdir()
    wl=root/'Lagrangian/RunGeneralT3FermionFirstEFT1Tree.wl';wl.write_text('WL contents')
    matching=root/'Lagrangian/RunMatching.wl';matching.write_text('Match[]')
    parms=root/'output/parameters.json';parms.write_text('{}')
    prior=root/'output/prior.json';prior.write_text('CERTIFIED')
    monkeypatch.setattr(task,'ROOT',root)
    monkeypatch.setattr(task,'WOLFRAM',wl)
    monkeypatch.setattr(task,'MATCHING_SOURCE',matching)
    monkeypatch.setattr(task,'prepare_model',lambda *a,**k:{'payload':{'unitary_U':{}},'parameters':parms})
    monkeypatch.setattr(task,'_prior_certificate',lambda *a,**k:prior)
    monkeypatch.setattr(task.shutil,'which',lambda _: 'wolframscript')
    def pretend(command,**kwargs):
        assert 'T3_EXPERIMENTAL_GENERAL_QUARTICS' not in kwargs['env']
        target=Path(command[-1]);target.write_text(json.dumps(certified(stem=target.stem)))
        write_artifacts(target.parent,target.stem)
        return SimpleNamespace(returncode=0,stdout='tree EFT1 passed',stderr='')
    monkeypatch.setattr(task.subprocess,'run',pretend)
    result=task.run_model('B_-1',outdir=root/'output/eft1')
    assert result['llss_term_count']==4
    assert result['sequential_EFT1_completed'] is False
    assert result['physical_C5_computed'] is False
    assert (root/'output/eft1/B_m1_final_eft1_fermion_first_tree.log').exists()


def test_nonzero_process_error_refused(monkeypatch,tmp_path):
    root=tmp_path/'project';(root/'Lagrangian').mkdir(parents=True);(root/'output').mkdir()
    wl=root/'Lagrangian/RunGeneralT3FermionFirstEFT1Tree.wl';wl.write_text('wl')
    mt=root/'Lagrangian/RunMatching.wl';mt.write_text('matching')
    parms=root/'output/input.json';parms.write_text('{}')
    prior=root/'output/prior.json';prior.write_text('check')
    for field,value in [('ROOT',root),('WOLFRAM',wl),('MATCHING_SOURCE',mt)]:monkeypatch.setattr(task,field,value)
    monkeypatch.setattr(task,'prepare_model',lambda *a,**k:{'payload':{'unitary_U':{}},'parameters':parms})
    monkeypatch.setattr(task,'_prior_certificate',lambda *a,**k:prior)
    monkeypatch.setattr(task.shutil,'which',lambda _:'wolframscript')
    def nope(command,**kwargs):
        target=Path(command[-1]);target.write_text(json.dumps(certified(stem=target.stem)))
        return SimpleNamespace(returncode=5,stdout='',stderr='Failed to match')
    monkeypatch.setattr(task.subprocess,'run',nope)
    with pytest.raises(ValueError,match='EFT1 incomplete'):
        task.run_model('B_-1',outdir=root/'output/eft1')


def test_eom_transition_required_and_clean(tmp_path):
    stem='B_m1_final_eft1_fermion_first_tree'
    write_artifacts(tmp_path,stem)
    p=tmp_path/(stem+'_eft1_transition_eom.inputform.txt')
    p.unlink()
    with pytest.raises(ValueError,match='missing canonical'):
        task.validate_tree_report(certified(), 'B_-1',tmp_path,stem)
    p.write_text('GammaCC Field[NewScalarX] Field[NewScalarY] Field[NewFermion]')
    with pytest.raises(ValueError,match='integrated F'):
        task.validate_tree_report(certified(), 'B_-1',tmp_path,stem)
    p.write_text('GammaCC Field[NewScalarX] Field[NewScalarY]')
    bad=certified()
    bad['transition_representation']='MatchedEFT'
    with pytest.raises(ValueError):
        task.validate_tree_report(bad,'B_-1',tmp_path,stem)
