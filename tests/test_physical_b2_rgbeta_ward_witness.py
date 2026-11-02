"""Synthetic exact Ward, negative controls, source-provenance and promotion guards."""
from __future__ import annotations
import json
from pathlib import Path
import sys
import pytest
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'Numerical/lo_certification/physical_normalisation'))
from b2_gauge_coupling import AuditBlocked
from b2_rgbeta_ward_witness import (
    ANTI,FUND,EPS,all16_charge_warning,compute,delta,fit_tensor,
    formal_evidence,invariant_identity,literal_B_y1,literal_B_y2,
    literal_rgbeta_22_mix,mix_charge_support,norm2,repaired_22_mix,
    repaired_B_y1,source_provenance,ward_residual,
)


def test_mix_literal_fails_su2_exact():
    assert ward_residual(literal_rgbeta_22_mix,(FUND,FUND,FUND,ANTI))==(3,0,3)


def test_mixed_delta_epsilon_exact_singlet():
    assert ward_residual(repaired_22_mix,(FUND,FUND,FUND,ANTI))==(0,0,0)


def test_B_y1_literal_epsilon_wrong_bar_orientation():
    assert ward_residual(literal_B_y1,(ANTI,FUND))==(2,0,2)
    assert ward_residual(repaired_B_y1,(ANTI,FUND))==(0,0,0)


def test_B_y2_epsilon_correct_for_two_fundamentals():
    assert ward_residual(literal_B_y2,(FUND,FUND))==(0,0,0)


def test_wrong_delta_for_two_fundamentals_fails():
    assert ward_residual(lambda ab:delta(*ab),(FUND,FUND))!=(0,0,0)


def test_repaired_mix_Bose_symmetric():
    for h in range(2):
        for hh in range(2):
            for s in range(2):
                for sb in range(2):
                    assert repaired_22_mix((h,hh,s,sb))==repaired_22_mix((hh,h,s,sb))


def test_bad_charge_all_B_C_alphas():
    rows=all16_charge_warning()
    assert len(rows)==6
    assert set(r['alpha'] for r in rows)=={-3,-1,1}
    for r in rows:
        assert any(not x['electric_charge_preserved'] for x in r['literal_source_tensor'])
        assert all(x['electric_charge_preserved'] for x in r['explicit_Ward_invariant_candidate'])


def test_nonzero_H0H0_good_vs_bad():
    a=-1
    bad=mix_charge_support(a,literal_rgbeta_22_mix)
    good=mix_charge_support(a,repaired_22_mix)
    assert [(x['S1_weight_index_1based'],x['S2_weight_index_1based']) for x in bad]==[(2,2)]
    assert [(x['S1_weight_index_1based'],x['S2_weight_index_1based']) for x in good]==[(1,2)]


def test_exact_fit_accepts_phase_and_rejects_wrong_tensor():
    factor=sp.I*sp.sqrt(3)/2
    good=fit_tensor(repaired_22_mix,lambda ix:factor*repaired_22_mix(ix),4)
    assert good['is_proportional_exact'] and good['complex_ratio_to_chosen_candidate']==str(factor)
    bad=fit_tensor(repaired_22_mix,literal_rgbeta_22_mix,4)
    assert not bad['is_proportional_exact']
    assert bad['ratio_is_a_physical_RGBeta_coupling_conversion'] is False


def test_ward_summary_is_not_physical_certification():
    out=invariant_identity()
    assert len(out['tensor_checks'])==5
    assert not out['native_RGBeta_index_semantics_verified']
    assert not out['native_RGBeta_action_prefactors_verified']


def test_reject_missing_source(tmp_path):
    with pytest.raises(AuditBlocked):source_provenance(tmp_path)


def test_reject_ancestor_repo_as_project_root(tmp_path):
    import subprocess
    repo=tmp_path/'outer';repo.mkdir()
    subprocess.run(['git','init','-q',str(repo)],check=True)
    child=repo/'child';child.mkdir()
    with pytest.raises(AuditBlocked,match='exact requested repository root'):
        source_provenance(child)


def test_reject_wrong_formal_status():
    with pytest.raises(AuditBlocked):formal_evidence({
        'format':'b2-native-rgbeta-formal-invariant-probe-v1',
        'physical_source_certified':True})


def test_reject_fake_physical_registry(tmp_path):
    reg=tmp_path/'reg.json'
    reg.write_text(json.dumps({'format':'b2-rgbeta-extended-native-registry-v1',
                              'physical_source_certified':True}),encoding='utf-8')
    with pytest.raises(AuditBlocked):compute(tmp_path,None,reg,validate_sources=False)


def test_synthetic_computation_stays_false():
    obj=compute(Path('.'),None,None,validate_sources=False)
    assert obj['physical_source_certified'] is False
    assert obj['physical_crosspackage_complete'] is False
    assert obj['crosspackage_yukawa_ratios'] is None
    assert obj['crosspackage_mass_ratio'] is None
    assert obj['crosspackage_mixing_ratio'] is None
    assert not obj['source_bug_confirmed']


def test_conjugation_negative_control():
    # Wrongly classifying Bar[S1] as unbarred would reverse the conclusion.
    assert ward_residual(literal_B_y1,(FUND,FUND))==(0,0,0)
    assert ward_residual(repaired_B_y1,(FUND,FUND))!=(0,0,0)


def _mini_git_repo(tmp_path, monkeypatch, *, worktree_crlf=False):
    """Real mini Git repo, with a Git-*canonical* source blob expectation.

    On Windows, pathlib.write_text() can create CRLF bytes while `git add`
    stores LF. The old fixture incorrectly used git_blob(path.read_bytes())
    as the committed reference hash, making one test fail on Windows despite
    the source_provenance() verifier behaving correctly.
    """
    import subprocess
    import b2_rgbeta_ward_witness as witness
    from b2_gauge_coupling import git_blob

    root=tmp_path/'mini';root.mkdir()
    rel='RGE/running/rgbeta/T3RGBetaModel.wl'
    path=root/rel;path.parent.mkdir(parents=True)
    content=(
        'T3RGBetaMixInvariant[2, 2] := (\n'
        ' delS2[SU2L, a, #1, #2] * delS2[SU2L, a, #3, #4] &);\n'
        'T3RGBetaYukawaInvariant[2, 1] :=\n'
        ' (eps[SU2L @ fund, #1, #2] &);\n'
        'AddQuartic[lambdaT3, {H, H, S1, Bar @ S2}, Test -> 1];\n'
    ).encode('utf-8')
    # Write exact fixture bytes; do not depend on OS text-mode translation.
    path.write_bytes(content.replace(b'\n',b'\r\n') if worktree_crlf else content)

    def git(*args):
        return subprocess.run(['git','-C',str(root),*args],check=True,
                              capture_output=True,text=True).stdout.strip()
    git('init','-q')
    # Test source_provenance() in the same canonical LF / optional CRLF
    # conditions encountered with Git for Windows.
    git('config','core.autocrlf','true')
    git('add','--',rel)
    git('-c','user.email=test@example.invalid','-c','user.name=Test',
        'commit','-qm','fixture')
    committed_blob=git('rev-parse','HEAD:'+rel)
    assert git('hash-object','--path='+rel,'--',rel)==committed_blob
    if worktree_crlf:
        assert git_blob(path.read_bytes())!=committed_blob
    # Compare with actual Git's committed blob, not working-tree raw bytes.
    monkeypatch.setattr(witness,'SOURCE_BLOBS',{rel:committed_blob})
    return root,path,git


def test_git_source_guard_accepts_clean_crlf_worktree(tmp_path,monkeypatch):
    root,path,_git=_mini_git_repo(tmp_path,monkeypatch,worktree_crlf=True)
    assert b'\r\n' in path.read_bytes()
    report=source_provenance(root)
    assert report['all_six_source_blobs_unchanged']
    assert report['validated_git_blobs'][0]['validated_with_git_clean_filter']


def test_git_source_guard_accepts_new_head_when_source_blobs_same(tmp_path,monkeypatch):
    root,path,git=_mini_git_repo(tmp_path,monkeypatch)
    first=source_provenance(root)
    assert first['all_six_source_blobs_unchanged']  # fixture uses the same verification logic
    (root/'README').write_text('unrelated commit',encoding='utf-8')
    git('add','README')
    git('-c','user.email=test@example.invalid','-c','user.name=Test',
        'commit','-qm','new commit')
    later=source_provenance(root)
    assert later['repository_head_measured']!=first['repository_head_measured']
    assert later['validated_git_blobs'][0]['git_blob_sha']==first['validated_git_blobs'][0]['git_blob_sha']


def test_git_source_guard_rejects_real_worktree_change(tmp_path,monkeypatch):
    root,path,git=_mini_git_repo(tmp_path,monkeypatch)
    path.write_text(path.read_text()+'\n(* real source drift *)\n')
    with pytest.raises(AuditBlocked,match='changed/unverified'):
        source_provenance(root)


def test_git_source_guard_rejects_staged_change(tmp_path,monkeypatch):
    root,path,git=_mini_git_repo(tmp_path,monkeypatch)
    path.write_text(path.read_text()+'\n(* staged drift *)\n')
    git('add',str(path.relative_to(root)))
    with pytest.raises(AuditBlocked,match='changed/unverified'):
        source_provenance(root)


def test_git_source_guard_rejects_changed_head_tree(tmp_path,monkeypatch):
    root,path,git=_mini_git_repo(tmp_path,monkeypatch)
    path.write_text(path.read_text()+'\n(* committed drift *)\n')
    git('add',str(path.relative_to(root)))
    git('-c','user.email=test@example.invalid','-c','user.name=Test',
        'commit','-qm','source drift')
    with pytest.raises(AuditBlocked,match='changed/unverified'):
        source_provenance(root)


def test_formal_evidence_actual_BC_expression_contract():
    from b2_gauge_coupling import CLASSES
    records=[]
    form='(del[SU2L[fund], b2h, b2t]*del[SU2L[fund], b2h2, b2s])/2 + (del[SU2L[fund], b2h, b2s]*del[SU2L[fund], b2h2, b2t])/2'
    for family,(d1,d2,df,alphas) in CLASSES.items():
        records.append({
            'family':family,'status':'RegistryAndFormalInvariantCaptured',
            'dimensions_alpha':[d1,d2,df,alphas[0]],
            'physical_action_certified':False,
            'native_numeric_RGBeta_component_axes_certified':False,
            'mixing_invariant_formal':form if family in ('B','C') else 'pending',
            'registered_yukawas':[
                {'coupling':'y1','actual_registered_fields':'{Bar[S1], l, F}',
                 'formal_invariant_inputform':
                     'eps[SU2L[fund], b2a, b2b]' if family=='B'
                     else 'delS2[SU2L, b2c, b2a, b2b]'},
                {'coupling':'y2','actual_registered_fields':'{S2, l, F}',
                 'formal_invariant_inputform':'pending'}]})
    fake={'format':'b2-native-rgbeta-formal-invariant-probe-v1',
          'status':'FormalRegistriesCapturedPhysicalNormalisationPending',
          'physical_source_certified':False,
          'numerical_component_axes_verified':False,
          'full_RGBeta_action_normalisation_verified':False,
          'five_classes':records}
    assert formal_evidence(fake)['B_C_formal_delta_delta_confirmed']
    records[1]['mixing_invariant_formal']='eps[SU2L[fund],x,y]'
    with pytest.raises(AuditBlocked,match='formal expansion changed'):
        formal_evidence(fake)
