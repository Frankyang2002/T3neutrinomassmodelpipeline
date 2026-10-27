"""Offline boundary tests; genuine Matchete scalar Match must run locally."""
from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from scripts import RunGeneralT3EFT1ScalarThresholdProbe as task


def write_artifacts(folder: Path, stem: str, *, count=2):
    folder.mkdir(parents=True, exist_ok=True)
    (folder/(stem+'_sm_tree.inputform.txt')).write_text('Field[H, Scalar]')
    (folder/(stem+'_hermiticity_diagnostics.json')).write_text(json.dumps({
        'final_hermiticity': True,
        'final_ordinary_CheckLagrangian': 'True',
        'equivalence_rewrite_accepted': 'None',
        'physical_sequential_C5_computed': False,
    }))
    (folder/(stem+'_sm_tree_plus_loop.inputform.txt')).write_text('hbar GammaCC Field[H, Scalar]')
    (folder/(stem+'_scalar_one_loop_correction.inputform.txt')).write_text('hbar GammaCC')
    (folder/(stem+'_raw_scalar_llhh.inputform.txt')).write_text(
        'hbar*GammaCC*Field[H, Scalar]*Field[l, Fermion]' if count else '0'
    )


def certificate(model='B_-1', stem='B_m1_final_eft1_scalar_threshold_probe',
                input_name='source_eft1_tree.inputform.txt', count=2):
    return {
        'format': task.FORMAT, 'status': task.STATUS, 'model': model,
        'input_eft1_tree': input_name,
        'input_eft1_tree_character_count': 40,
        'eigenfield_input_yukawas_and_quartics_reconstructed': True,
        'fermion_removed_before_scalar_matching': True,
        'heavy_scalar_fields': ['NewScalarX', 'NewScalarY'],
        'scalar_heavy_mass_basis_detailed_check_passed': True,
        'hermiticity_check_passed': True,
        'canonical_EOMEFT_hermiticity_rewrite': 'None',
        'hermiticity_diagnostic_file': stem+'_hermiticity_diagnostics.json',
        'tree_matching_completed': True, 'one_loop_matching_completed': True,
        'scalar_one_loop_correction_hbar_present': True,
        'raw_scalar_llhh_term_count': count,
        'scalar_loop_difference_term_count': 8,
        'sm_tree_file': stem+'_sm_tree.inputform.txt',
        'sm_one_loop_file': stem+'_sm_tree_plus_loop.inputform.txt',
        'scalar_correction_file': stem+'_scalar_one_loop_correction.inputform.txt',
        'raw_scalar_llhh_file': stem+'_raw_scalar_llhh.inputform.txt',
        'no_scalar_fields_in_matched_eft': True,
        'eft1_running_performed': False,
        'hard_pole_counterterms_combined': False,
        'hard_pole_rge_certified': False,
        'physical_sequential_C5_computed': False,
        'physical_matching_performed': False,
    }


def make_tree(tmp_path):
    path = tmp_path/'source_eft1_tree.inputform.txt'
    path.write_text('0123456789' * 4)
    return path


@pytest.mark.parametrize('model', ['B_-1', 'C_-1'])
def test_model_stem(model):
    assert task._stem(model, 'final').startswith(model.replace('-', 'm'))


@pytest.mark.parametrize('model,which', [('A_+0','final'),('B_-1','unknown')])
def test_wrong_model(model, which):
    with pytest.raises(ValueError,match='Unsupported'):
        task._stem(model,which)


def test_report_accepted(tmp_path):
    stem=task._stem('B_-1','final')
    tree=make_tree(tmp_path)
    write_artifacts(tmp_path,stem)
    files=task.validate_report(certificate(), 'B_-1', tmp_path,stem,tree)
    assert len(files)==4
    assert all(x.exists() for x in files.values())


def test_zero_raw_sector_allowed(tmp_path):
    stem=task._stem('B_-1','final')
    tree=make_tree(tmp_path)
    write_artifacts(tmp_path,stem,count=0)
    assert len(task.validate_report(certificate(count=0),'B_-1',tmp_path,stem,tree))==4


@pytest.mark.parametrize('name,bad', [
    ('status','Done'),('format','Wrong'),('model','C_-1'),
    ('fermion_removed_before_scalar_matching',False),
    ('heavy_scalar_fields',['NewScalarX']),
    ('scalar_heavy_mass_basis_detailed_check_passed',False),
    ('hermiticity_check_passed',False),
    ('canonical_EOMEFT_hermiticity_rewrite','ManualHc'),
    ('one_loop_matching_completed',False),
    ('raw_scalar_llhh_term_count',-1),
    ('raw_scalar_llhh_term_count',True),
    ('scalar_loop_difference_term_count',0),
    ('scalar_loop_difference_term_count',False),
    ('scalar_one_loop_correction_hbar_present',False),
    ('no_scalar_fields_in_matched_eft',False),
    ('eft1_running_performed',True),
    ('hard_pole_counterterms_combined',True),
    ('hard_pole_rge_certified',True),
    ('physical_sequential_C5_computed',True),
    ('physical_matching_performed',True),
    ('input_eft1_tree_character_count',41),
])
def test_report_fail_closed(tmp_path,name,bad):
    stem=task._stem('B_-1','final')
    tree=make_tree(tmp_path)
    write_artifacts(tmp_path,stem)
    report=certificate();report[name]=bad
    with pytest.raises(ValueError):
        task.validate_report(report,'B_-1',tmp_path,stem,tree)


@pytest.mark.parametrize('artifact,text', [
    ('_sm_tree_plus_loop.inputform.txt','hbar Field[NewScalarX]'),
    ('_sm_tree_plus_loop.inputform.txt','hbar Field[NewFermion]'),
    ('_scalar_one_loop_correction.inputform.txt','tree only'),
    ('_raw_scalar_llhh.inputform.txt','Field[NewScalarY] GammaCC'),
    ('_raw_scalar_llhh.inputform.txt','hbar GammaCC'),
])
def test_report_artifacts_fail_closed(tmp_path,artifact,text):
    stem=task._stem('B_-1','final')
    tree=make_tree(tmp_path)
    write_artifacts(tmp_path,stem)
    (tmp_path/(stem+artifact)).write_text(text)
    with pytest.raises(ValueError):
        task.validate_report(certificate(),'B_-1',tmp_path,stem,tree)


def test_wrong_tree_path_or_name(tmp_path):
    stem=task._stem('B_-1','final')
    tree=make_tree(tmp_path)
    write_artifacts(tmp_path,stem)
    c=certificate();c['input_eft1_tree']='different_tree.inputform.txt'
    with pytest.raises(ValueError,match='input_eft1_tree'):
        task.validate_report(c,'B_-1',tmp_path,stem,tree)


def test_wl_source_is_second_threshold_and_fail_closed():
    path=Path(__file__).resolve().parents[2]/'Lagrangian/RunGeneralT3EFT1ScalarThresholdProbe.wl'
    wl=path.read_text()
    assert 'Mass->{Light,MSX}' in wl and 'Mass->{Light,MSY}' in wl
    assert 'Mass->{Heavy,MSX}' in wl and 'Mass->{Heavy,MSY}' in wl
    assert 'RemoveField[NewFermion]' in wl
    assert 'CheckAbort[CheckLagrangian[scalarEFT],$Aborted]' in wl
    assert 'Quiet@Check[CheckLagrangian[scalarEFT]' not in wl
    assert 'RunT3Matching[scalarEFT,5,0]' in wl
    assert 'RunT3Matching[scalarEFT,5,1]' in wl
    assert 'correction=Expand[loopMatched-treeMatched]' in wl
    assert '"physical_sequential_C5_computed"->False' in wl
    assert '"hard_pole_rge_certified"->False' in wl
    assert 'UsingFrontEnd[Catch[T3MEScalarThresholdPerform[]' in wl
    assert 'T3MECheckpoint["scalar-threshold-one-loop-Match-started"]' in wl
    assert 'HcSimplify[scalarEFT]' in wl
    assert 'GreensSimplify[scalarEFT]' in wl
    assert '"light_EFT1_hermiticity"' in wl
    assert '"heavy_EFT1_original_hermiticity"' in wl
    assert '"equivalence_rewrite_accepted"' in wl
    assert '"hermiticity_check_passed"->True' in wl
    assert 'scalarEFT=candidate' in wl
    assert 'PlusHc[scalarEFT]' not in wl
    assert 'T3MEPerform[] :=' not in wl


def test_missing_verified_first_stage(monkeypatch,tmp_path):
    with pytest.raises(ValueError,match='missing live tree EFT1 certificate'):
        task._verified_tree_input('B_-1','final',tmp_path/'tree',tmp_path/'loop')


def test_mock_wolfram_run(monkeypatch,tmp_path):
    root=tmp_path/'repo'
    (root/'Lagrangian').mkdir(parents=True)
    (root/'output').mkdir()
    wl=root/'Lagrangian/RunGeneralT3EFT1ScalarThresholdProbe.wl';wl.write_text('Matchete source')
    match=root/'Lagrangian/RunMatching.wl';match.write_text('RunT3Matching')
    params=root/'output/parameters.json';params.write_text('{}')
    mass=root/'output/mass.json';mass.write_text('{}')
    tcert=root/'output/tree.json';tcert.write_text('{}')
    lcert=root/'output/oneloop.json';lcert.write_text('{}')
    tree=root/'output/source_eft1_tree.inputform.txt';tree.write_text('0123456789'*4)
    for name,val in [('ROOT',root),('WOLFRAM',wl),('MATCHING_SOURCE',match)]:
        monkeypatch.setattr(task,name,val)
    monkeypatch.setattr(task,'prepare_model',lambda *a,**k:{'payload':{'unitary_U':{}},'parameters':params})
    monkeypatch.setattr(task,'_prior_certificate',lambda *a,**k:mass)
    monkeypatch.setattr(task,'_verified_tree_input',lambda *a,**k:(tree,tcert,lcert))
    monkeypatch.setattr(task.shutil,'which',lambda x:'wolframscript')
    def fake_run(cmd,**kwargs):
        assert cmd[-2]==str(tree)
        assert 'T3_EXPERIMENTAL_GENERAL_QUARTICS' not in kwargs['env']
        target=Path(cmd[-1]);stem=target.stem
        write_artifacts(target.parent,stem)
        target.write_text(json.dumps(certificate(stem=stem)))
        return SimpleNamespace(returncode=0,stdout='scalar success',stderr='')
    monkeypatch.setattr(task.subprocess,'run',fake_run)
    result=task.run_model('B_-1',outdir=root/'output/scalar')
    assert result['raw_scalar_llhh_term_count']==2
    assert result['physical_sequential_C5_computed'] is False


def test_failed_wolfram_exit_not_accepted(monkeypatch,tmp_path):
    root=tmp_path/'repo';(root/'Lagrangian').mkdir(parents=True);(root/'output').mkdir()
    wl=root/'Lagrangian/RunGeneralT3EFT1ScalarThresholdProbe.wl';wl.write_text('Matchete source')
    match=root/'Lagrangian/RunMatching.wl';match.write_text('RunT3Matching')
    parms=root/'output/parameters.json';parms.write_text('{}')
    mass=root/'output/mass.json';mass.write_text('{}')
    tcert=root/'output/tree.json';tcert.write_text('{}')
    lcert=root/'output/oneloop.json';lcert.write_text('{}')
    tree=root/'output/source_eft1_tree.inputform.txt';tree.write_text('0123456789'*4)
    for name,val in [('ROOT',root),('WOLFRAM',wl),('MATCHING_SOURCE',match)]:
        monkeypatch.setattr(task,name,val)
    monkeypatch.setattr(task,'prepare_model',lambda *a,**k:{'payload':{'unitary_U':{}},'parameters':parms})
    monkeypatch.setattr(task,'_prior_certificate',lambda *a,**k:mass)
    monkeypatch.setattr(task,'_verified_tree_input',lambda *a,**k:(tree,tcert,lcert))
    monkeypatch.setattr(task.shutil,'which',lambda x:'wolframscript')
    def mock_run(cmd,**kwargs):
        target=Path(cmd[-1]);stem=target.stem
        write_artifacts(target.parent,stem)
        target.write_text(json.dumps(certificate(stem=stem)))
        return SimpleNamespace(returncode=5,stdout='check failed',stderr='failed')
    monkeypatch.setattr(task.subprocess,'run',mock_run)
    with pytest.raises(ValueError,match='failed scalar-threshold matching'):
        task.run_model('B_-1',outdir=root/'output/scalar')


def test_canonical_transition_from_validated_tree(monkeypatch,tmp_path):
    src=tmp_path/'tree';src.mkdir()
    one=tmp_path/'one';one.mkdir()
    tstem=task.tree_stem('B_-1','final')
    lstem=task.one_loop_stem('B_-1','final')
    tree_report={'transition_representation':'EOMEFT',
        'transition_eom_eft1_file':tstem+'_eft1_transition_eom.inputform.txt'}
    (src/(tstem+'.json')).write_text(json.dumps(tree_report))
    (one/(lstem+'.json')).write_text('{}')
    transition=src/tree_report['transition_eom_eft1_file']
    transition.write_text('GammaCC Field[NewScalarX] Field[NewScalarY] '*50)
    monkeypatch.setattr(task,'validate_tree_report',lambda *a,**k:None)
    monkeypatch.setattr(task,'validate_one_loop_report',lambda *a,**k:{'matched_eft_file':tmp_path/'old_expanded_eft.txt'})
    result,_,_=task._verified_tree_input('B_-1','final',src,one)
    assert result==transition
    tree_report['transition_representation']='MatchedEFT'
    (src/(tstem+'.json')).write_text(json.dumps(tree_report))
    with pytest.raises(ValueError,match='canonical EOMEFT'):
        task._verified_tree_input('B_-1','final',src,one)


def test_diagnostic_hard_failures(tmp_path):
    stem=task._stem('B_-1','final')
    tree=make_tree(tmp_path)
    write_artifacts(tmp_path,stem)
    diag_path=tmp_path/(stem+'_hermiticity_diagnostics.json')
    diag=json.loads(diag_path.read_text())
    for bad in ({'final_hermiticity': False},
                {'final_ordinary_CheckLagrangian': 'False'},
                {'physical_sequential_C5_computed': True},
                {'equivalence_rewrite_accepted': 'HcSimplify'}):
        new=diag.copy();new.update(bad)
        diag_path.write_text(json.dumps(new))
        with pytest.raises(ValueError, match='Hermiticity diagnostic'):
            task.validate_report(certificate(),'B_-1',tmp_path,stem,tree)
    diag_path.write_text(json.dumps(diag))
    assert task.validate_report(certificate(),'B_-1',tmp_path,stem,tree)


def test_identity_rewrite_still_requires_strict_certificate(tmp_path):
    stem=task._stem('C_-1','final')
    tree=make_tree(tmp_path)
    write_artifacts(tmp_path,stem)
    diag_path=tmp_path/(stem+'_hermiticity_diagnostics.json')
    diag=json.loads(diag_path.read_text())
    diag['equivalence_rewrite_accepted']='GreensSimplify'
    diag_path.write_text(json.dumps(diag))
    c=certificate(model='C_-1',stem=stem)
    c['canonical_EOMEFT_hermiticity_rewrite']='GreensSimplify'
    assert task.validate_report(c,'C_-1',tmp_path,stem,tree)


def test_diagnostic_path_required(tmp_path):
    stem=task._stem('B_-1','final')
    tree=make_tree(tmp_path)
    write_artifacts(tmp_path,stem)
    (tmp_path/(stem+'_hermiticity_diagnostics.json')).unlink()
    with pytest.raises(ValueError,match='Hermiticity diagnostic not exported'):
        task.validate_report(certificate(),'B_-1',tmp_path,stem,tree)

@pytest.mark.parametrize('method', ['OriginalRelabel', 'ExpandOnly', 'ExpandRelabel'])
def test_identity_preserving_rewrites_have_strict_certificates(tmp_path, method):
    stem = task._stem('B_-1', 'final')
    tree = make_tree(tmp_path)
    write_artifacts(tmp_path, stem)
    p = tmp_path / (stem + '_hermiticity_diagnostics.json')
    payload = json.loads(p.read_text())
    payload['equivalence_rewrite_accepted'] = method
    p.write_text(json.dumps(payload))
    report = certificate()
    report['canonical_EOMEFT_hermiticity_rewrite'] = method
    task.validate_report(report, 'B_-1', tmp_path, stem, tree)


def test_wolfram_sector_audit_does_not_modify_eft():
    root = Path(__file__).resolve().parents[2]
    wl = (root / 'Lagrangian/RunGeneralT3EFT1ScalarThresholdProbe.wl').read_text()
    assert 'T3MEHermiticitySectorAudit[eft_, dest_String]' in wl
    assert 'scalar-threshold-hermiticity-sector-audit-started' in wl
    for fragment in ('"background_no_LLSS_or_BSM_quartics"',
                     '"background_plus_all_LLSS"',
                     '"background_plus_all_BSM_quartics"',
                     '"background_plus_XX_LLSS"',
                     '"background_plus_XY_LLSS"',
                     '"background_plus_YY_LLSS"',
                     '"background_plus_scalar_self_quartics"',
                     '"background_plus_Higgs_portal_quartics"',
                     '"partition_conserves_exact_terms"->True',
                     '"physical_sequential_C5_computed"->False'):
        assert fragment in wl
    assert 'scalarEFT,hermiticitySectorPath' in wl
    assert 'scalarEFT=Total' not in wl
    assert 'scalarEFT=PlusHc' not in wl
    assert '1/2*(scalarEFT' not in wl
    runner = (root/'scripts/RunGeneralT3EFT1ScalarThresholdProbe.py').read_text()
    assert "'_hermiticity_sector_audit.json'" in runner


def test_strict_renormalizable_reconstruction_acceptance(tmp_path):
    stem = task._stem('B_-1', 'final')
    tree = make_tree(tmp_path)
    write_artifacts(tmp_path, stem)
    c = certificate()
    c['canonical_EOMEFT_hermiticity_rewrite'] = 'VerifiedRenormalizableReconstruction'
    c['reconstruction_audit_file'] = stem + '_reconstruction_audit.json'
    diag = tmp_path / (stem + '_hermiticity_diagnostics.json')
    payload = json.loads(diag.read_text())
    payload['equivalence_rewrite_accepted'] = 'VerifiedRenormalizableReconstruction'
    diag.write_text(json.dumps(payload))
    proof = tmp_path / c['reconstruction_audit_file']
    proof.write_text(json.dumps({
        'status': 'VerifiedEquivalentReconstructionAccepted',
        'algebraic_equality_proven': True,
        'independent_matchete_check_passed': True,
        'accepted': True,
        'physical_sequential_C5_computed': False,
    }))
    task.validate_report(c, 'B_-1', tmp_path, stem, tree)
    for key, invalid in (
        ('algebraic_equality_proven', False),
        ('independent_matchete_check_passed', False),
        ('accepted', False),
        ('physical_sequential_C5_computed', True),
    ):
        data = json.loads(proof.read_text())
        data[key] = invalid
        proof.write_text(json.dumps(data))
        with pytest.raises(ValueError, match='unverified reconstruction'):
            task.validate_report(c, 'B_-1', tmp_path, stem, tree)
        data[key] = not invalid
        proof.write_text(json.dumps(data))
    proof.unlink()
    with pytest.raises(ValueError,match='audit missing'):
        task.validate_report(c, 'B_-1', tmp_path, stem, tree)


def test_wolfram_reconstruction_cannot_skip_equality_or_full_check():
    wl = (Path(__file__).resolve().parents[2]
          / 'Lagrangian/RunGeneralT3EFT1ScalarThresholdProbe.wl').read_text()
    assert 'T3MEVerifiedRenormalizableReconstruction[' in wl
    assert 'TrueQ[residual===0]' in wl
    assert 'TimeConstrained[Simplify[relabeled==0],15,False]' in wl
    assert 'renormFree=CheckAbort[FreeLag[NewScalarX,NewScalarY]' in wl
    assert 'Total[majorana]' in wl
    assert 'clean=exact && TrueQ[fullCheck]' in wl
    assert '"scalar_self_quartic_coupling_isolation"->qCouplingResults' in wl
    assert '"physical_sequential_C5_computed"->False' in wl
    assert 'scalarEFT=reconstruction["candidate"]' in wl
    assert 'scalarEFT=(scalarEFT+Hc' not in wl
