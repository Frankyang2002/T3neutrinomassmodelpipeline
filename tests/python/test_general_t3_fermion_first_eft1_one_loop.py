"""Offline contracts only; full Matchete Match requires the user's local kernel."""
from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from scripts import RunGeneralT3FermionFirstEFT1OneLoop as task


def certificate(model='B_-1', stem='B_m1_final_eft1_fermion_first_oneloop', n=2):
    return {
        'format': task.FORMAT,
        'status': task.STATUS, 'model': model,
        'matching_approach': 'FermionFirstTreeAndOneLoopScalarsRemainLight',
        'matching_EFT_order': 5, 'matching_loop_order': 1,
        'mass_basis_ordinary_check_passed': True,
        'fermion_heavy_mass_basis_detailed_check_passed': True,
        'source_quartic_count': 17, 'source_yukawa_count': 2,
        'scalar_light_fields': ['NewScalarX', 'NewScalarY'],
        'integrated_heavy_fields': ['NewFermion'],
        'tree_matching_completed': True, 'fermion_field_absent': True,
        'x_scalar_present': True, 'y_scalar_present': True,
        'llss_contains_x': True, 'llss_contains_y': True,
        'llss_term_count': 20, 'one_loop_matching_completed': True,
        'loop_difference_hbar_present': True,
        'loop_difference_fermion_absent': True,
        'loop_difference_term_count': 7,
        'raw_hard_llhh_candidate_term_count': n,
        'raw_hard_llhh_no_fermion_field': True,
        'raw_hard_llhh_scheme_renormalized': False,
        'raw_hard_llhh_pole_rge_certified': False,
        'eft1_rge_performed': False,
        'scalar_threshold_matching_performed': False,
        'one_loop_hard_fermion_matching_performed': False,
        'physical_C5_computed': False,
        'sequential_EFT1_completed': False,
        'physical_matching_performed': False,
        'matchete_rgbeta_quartic_normalization_certified': False,
        'matched_eft_file': stem + '_eft1_tree.inputform.txt',
        'llss_inputform_file': stem + '_llss_tree.inputform.txt',
        'one_loop_matched_eft_file': stem + '_eft1_tree_plus_loop.inputform.txt',
        'one_loop_difference_file': stem + '_eft1_one_loop_correction.inputform.txt',
        'raw_hard_llhh_file': stem + '_raw_hard_llhh_sector.inputform.txt',
    }


def write_artifacts(folder: Path, stem: str, n: int = 2):
    folder.mkdir(parents=True, exist_ok=True)
    (folder/(stem+'_eft1_tree.inputform.txt')).write_text('Field[NewScalarX] + Field[NewScalarY]')
    (folder/(stem+'_llss_tree.inputform.txt')).write_text('GammaCC Field[NewScalarX] Field[NewScalarY]')
    (folder/(stem+'_eft1_tree_plus_loop.inputform.txt')).write_text('Field[NewScalarX] + hbar F')
    (folder/(stem+'_eft1_one_loop_correction.inputform.txt')).write_text('hbar GammaCC Field[H, Scalar] Field[H, Scalar]')
    (folder/(stem+'_raw_hard_llhh_sector.inputform.txt')).write_text('hbar GammaCC Field[H, Scalar] Field[H, Scalar]' if n else '0')


@pytest.mark.parametrize('model', ['B_-1','C_-1'])
def test_allowed_stems(model):
    assert task._stem(model, 'final') == model.replace('-','m') + '_final_eft1_fermion_first_oneloop'


@pytest.mark.parametrize('model,which', [('A_+0','final'),('B_-1','wrong')])
def test_other_models_refused(model,which):
    with pytest.raises(ValueError,match='Unsupported'):
        task._stem(model,which)


def test_valid_certificate_and_artifacts(tmp_path):
    stem=task._stem('B_-1','final')
    write_artifacts(tmp_path,stem)
    result=task.validate_one_loop_report(certificate(), 'B_-1',tmp_path,stem)
    assert len(result)==5
    assert result['raw_hard_llhh_file'].is_file()


def test_zero_hard_sector_is_explicit_and_accepted(tmp_path):
    stem=task._stem('B_-1','final')
    write_artifacts(tmp_path,stem,n=0)
    assert len(task.validate_one_loop_report(certificate(n=0),'B_-1',tmp_path,stem))==5


@pytest.mark.parametrize('key,bad', [
    ('format','wrong'), ('status','Success'),('model','C_-1'),
    ('matching_EFT_order',6),('matching_loop_order',0),
    ('fermion_field_absent',False),('x_scalar_present',False),
    ('y_scalar_present',False),('llss_term_count',0),
    ('one_loop_matching_completed',False),('loop_difference_hbar_present',False),
    ('loop_difference_fermion_absent',False),('loop_difference_term_count',0),
    ('raw_hard_llhh_candidate_term_count',-1),
    ('raw_hard_llhh_scheme_renormalized',True),
    ('raw_hard_llhh_pole_rge_certified',True),
    ('physical_C5_computed',True),('eft1_rge_performed',True),
    ('scalar_threshold_matching_performed',True),
    ('sequential_EFT1_completed',True),('physical_matching_performed',True),
    ('matchete_rgbeta_quartic_normalization_certified',True),
])
def test_report_fail_closed(tmp_path,key,bad):
    stem=task._stem('B_-1','final')
    write_artifacts(tmp_path,stem)
    report=certificate()
    report[key]=bad
    with pytest.raises(ValueError):
        task.validate_one_loop_report(report,'B_-1',tmp_path,stem)


@pytest.mark.parametrize('filename,newtext',[
    ('_eft1_one_loop_correction.inputform.txt','GammaCC'),
    ('_eft1_one_loop_correction.inputform.txt','hbar Field[NewFermion]'),
    ('_llss_tree.inputform.txt','No GammaCC'),
    ('_raw_hard_llhh_sector.inputform.txt','hbar Field[NewScalarY]'),
    ('_raw_hard_llhh_sector.inputform.txt','hbar Field[H]'),
])
def test_bad_artifact_refused(tmp_path,filename,newtext):
    stem=task._stem('B_-1','final')
    write_artifacts(tmp_path,stem)
    (tmp_path/(stem+filename)).write_text(newtext)
    with pytest.raises(ValueError):
        task.validate_one_loop_report(certificate(),'B_-1',tmp_path,stem)


def test_script_has_matching_power_counting_and_no_physical_claim():
    script=(Path(__file__).resolve().parents[2]/'Lagrangian/RunGeneralT3FermionFirstEFT1OneLoop.wl').read_text()
    assert 'Mass->{Light,MSX}' in script and 'Mass->{Light,MSY}' in script
    assert 'RunT3Matching[lag,5,0]' in script
    assert 'RunT3Matching[lag,5,1]' in script
    assert 'loopDifference=Expand[oneLoopMatchedEFT-matchedEFT]' in script
    assert 'raw_hard_llhh_scheme_renormalized"->False' in script
    assert '"physical_C5_computed"->False' in script
    assert 'RunSequentialT3Matching[lag' not in script
    assert 'T3MECheckpoint["fermion-first-one-loop-Match-started"]' in script
    assert 'T3MECheckpoint["fermion-first-one-loop-sector-extracted"]' in script


def test_tree_dependency_fail_closed(monkeypatch,tmp_path):
    monkeypatch.setattr(task,'prepare_model',lambda *a,**k: {'payload':{'unitary_U':{}},'parameters':tmp_path/'parameters.json'})
    monkeypatch.setattr(task,'_prior_certificate',lambda *a,**k:tmp_path/'prior.json')
    with pytest.raises(ValueError,match='run the successful tree EFT1 probe'):
        task.run_model('B_-1',tree_dir=tmp_path/'missing',outdir=tmp_path/'out')


def test_fake_successful_wolfram_execution(monkeypatch,tmp_path):
    root=tmp_path/'repo';(root/'Lagrangian').mkdir(parents=True);(root/'output').mkdir()
    wl=root/'Lagrangian/RunGeneralT3FermionFirstEFT1OneLoop.wl';wl.write_text('Matchete matching script')
    match=root/'Lagrangian/RunMatching.wl';match.write_text('MatchingStages[]')
    parms=root/'output/input.json';parms.write_text('{}')
    cert=root/'output/mass.json';cert.write_text('{}')
    tree=root/'output/tree.json';tree.write_text('{}')
    for key,value in [('ROOT',root),('WOLFRAM',wl),('MATCHING_SOURCE',match)]:
        monkeypatch.setattr(task,key,value)
    monkeypatch.setattr(task,'prepare_model',lambda *a,**k:{'payload':{'unitary_U':{}},'parameters':parms})
    monkeypatch.setattr(task,'_prior_certificate',lambda *a,**k:cert)
    monkeypatch.setattr(task,'_certify_prior_tree',lambda *a,**k:tree)
    monkeypatch.setattr(task.shutil,'which',lambda _:'wolframscript')
    def fake_run(command,**kwargs):
        assert 'T3_EXPERIMENTAL_GENERAL_QUARTICS' not in kwargs['env']
        target=Path(command[-1]);stem=target.stem
        target.write_text(json.dumps(certificate(stem=stem)))
        write_artifacts(target.parent,stem)
        return SimpleNamespace(returncode=0,stdout='one-loop completed',stderr='')
    monkeypatch.setattr(task.subprocess,'run',fake_run)
    result=task.run_model('B_-1',outdir=root/'output/eft1_oneloop')
    assert result['loop_difference_terms']==7
    assert result['raw_hard_llhh_candidate_terms']==2
    assert not result['physical_C5_computed']


def test_nonzero_exit_refused(monkeypatch,tmp_path):
    root=tmp_path/'repo';(root/'Lagrangian').mkdir(parents=True);(root/'output').mkdir()
    wl=root/'Lagrangian/RunGeneralT3FermionFirstEFT1OneLoop.wl';wl.write_text('Matchete script')
    match=root/'Lagrangian/RunMatching.wl';match.write_text('MatchingStages[]')
    parms=root/'output/input.json';parms.write_text('{}')
    cert=root/'output/mass.json';cert.write_text('{}')
    tree=root/'output/tree.json';tree.write_text('{}')
    for k,v in [('ROOT',root),('WOLFRAM',wl),('MATCHING_SOURCE',match)]:
        monkeypatch.setattr(task,k,v)
    monkeypatch.setattr(task,'prepare_model',lambda *a,**k:{'payload':{'unitary_U':{}},'parameters':parms})
    monkeypatch.setattr(task,'_prior_certificate',lambda *a,**k:cert)
    monkeypatch.setattr(task,'_certify_prior_tree',lambda *a,**k:tree)
    monkeypatch.setattr(task.shutil,'which',lambda _:'wolframscript')
    def failed(command,**kwargs):
        target=Path(command[-1]);stem=target.stem
        target.write_text(json.dumps(certificate(stem=stem)))
        write_artifacts(target.parent,stem)
        return SimpleNamespace(returncode=5,stdout='',stderr='Matchete failed')
    monkeypatch.setattr(task.subprocess,'run',failed)
    with pytest.raises(ValueError,match='one-loop Matchete incomplete'):
        task.run_model('B_-1',outdir=root/'output/eft1')
