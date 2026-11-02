"""B2 hostile symbolic tests (run with pytest)."""
import json
import math
import sys
from pathlib import Path
import pytest
import sympy as sp

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from Numerical.lo_certification.physical_normalisation.b2_gauge_coupling import (
    AuditBlocked, source_rotations, mass_metrics, models, charge_channels,
    spurion_phase_ward, portal_singlet_embedding, scalar_A_mixing_offdiag,
    exact, git_blob, source_inventory, native_tensor,
)


def generator_record():
    i=sp.I;r=sp.sqrt(2)
    u=sp.Matrix([[-i/r,0,i/r],[0,1,0],[1/r,0,1/r]])
    def encode(mat):return [[str(v) for v in row] for row in mat.tolist()]
    return {'status':'generator_algebra_verified_CG_basis_unverified',
            'representations':[
                {'dimension':2,'status':'generator_algebra_verified',
                 'T3_in_original_matchete_basis':encode(sp.diag(sp.Rational(1,2),-sp.Rational(1,2))),
                 'U_columns_are_weight_eigenvectors':encode(sp.eye(2))},
                {'dimension':3,'status':'generator_algebra_verified',
                 'T3_in_original_matchete_basis':encode(u*sp.diag(1,0,-1)*u.H),
                 'U_columns_are_weight_eigenvectors':encode(u)}]}


def test_ward_and_basis():
    u=source_rotations(generator_record())
    evidence=mass_metrics(u)
    assert evidence['triplet_metric_in_source_weight_basis']==[['0','0','1'],['0','1','0'],['1','0','0']]
    assert evidence['native_RGBeta_epsilon_orientation_measured'] is False
    assert sp.Matrix([[0,1],[-1,0]]).T != sp.Matrix([[0,1],[-1,0]])


def test_incorrect_triplet_rotation_blocked():
    item=generator_record()
    item['representations'][1]['U_columns_are_weight_eigenvectors'][0][0]='1'
    with pytest.raises(AuditBlocked):source_rotations(item)


def test_allowed_weight_rephasing_is_not_error():
    item=generator_record()
    d=item['representations'][1]
    d['U_columns_are_weight_eigenvectors']=[[f'({z})*I' if j==0 and z!='0' else z
                                           for j,z in enumerate(row)]
                                           for row in d['U_columns_are_weight_eigenvectors']]
    u=source_rotations(item)
    assert mass_metrics(u)['native_RGBeta_triplet_del_axes_measured'] is False


def test_all_sixteen_charge_skeletons():
    items=models()
    assert len(items)==16
    assert [m['label'] for m in items if m['majorana']]==['B_m1','C_m1']
    for m in items:
        for ch in charge_channels(m):
            assert ch['Q']==str(sp.Rational(ch['Q']))
            assert ch['f']<=m['dF'] and ch['s1']<=m['dS1'] and ch['s2']<=m['dS2']


def test_spurion_and_swapped_heavy_field_and_wrong_HC():
    for majorana in (True,False):
        evidence=spurion_phase_ward(majorana)
        assert evidence['spurion_charges']==[0,0,0,0]
    ev=spurion_phase_ward(False)
    charges=lambda key:sp.Matrix(ev[key])
    y2wrong=sp.Matrix([0,-1,0,1])  # swapped F -> Fc
    assert charges('charges_y1')+y2wrong-charges('charges_MF')-charges('charges_lambdaT3') != sp.zeros(4,1)
    assert sum((charges(k) for k in ('charges_y1','charges_y2')),
               sp.zeros(4,1))-charges('charges_MF')+charges('charges_lambdaT3') != sp.zeros(4,1)


def test_complex_mixing_potential_HC_orientation():
    lam=sp.Integer(2)+3*sp.I
    v=sp.Integer(3)
    x,y=scalar_A_mixing_offdiag(lam,v,sp.Integer(-1))
    assert x == sp.Rational(9,2)*(2+3*sp.I)
    assert y == sp.conjugate(x) and y != x


def test_hermitian_portal_nullspace_not_unique():
    ident=sp.Matrix(4,4,lambda i,j:int(i in (0,3) and j in (0,3)))
    result=portal_singlet_embedding([ident,sp.I*sp.diag(1,0,0,0)],2)
    assert result['unique'] is False
    assert result['rank']==1
    assert result['nullspace_basis']==[['0','1']]
    assert result['canonical_pivot_coefficients']==['1','0']
    assert not portal_singlet_embedding([2*ident],2)['canonical_pivot_coefficients']==['1']
    assert portal_singlet_embedding([2*ident],2)['canonical_pivot_coefficients']==['1/2']


def test_false_unique_portal_rejected():
    ident=sp.eye(4)
    with pytest.raises(AuditBlocked):
        portal_singlet_embedding([ident,ident],2)


def test_untrusted_value_and_drift_detected(tmp_path):
    with pytest.raises(AuditBlocked):exact('__import__("os").system("true")')
    assert git_blob(b'abc') != git_blob(b'ABC')
    record=source_inventory(tmp_path)
    assert record['source_reference_matches'] is False
    assert record['current_git_head_measured'] is False


def test_explicit_native_tensor_orientation_and_conjugation():
    u={1:sp.eye(1),2:sp.eye(2),3:sp.eye(3)}
    v={'dimensions':[2,2,1], 'components':[{'indices':[1,2,1],'cg':'I'},
                                            {'indices':[2,1,1],'cg':'-I'}]}
    tensor=native_tensor(v,(True,True,False),u)
    assert tensor[0,1,0]==sp.I
    assert tensor[1,0,0]==-sp.I
    vbad={'dimensions':[2,2,1], 'components':[{'indices':[3,1,1],'cg':'1'}]}
    with pytest.raises(AuditBlocked):native_tensor(vbad,(True,True,False),u)


def test_no_false_certification_from_missing_or_stale_native(tmp_path):
    from Numerical.lo_certification.physical_normalisation.b2_gauge_coupling import (
        compute, prior_portal_dictionary,
    )
    cross, channels=compute(tmp_path)
    assert cross['physical_source_certified'] is False
    assert cross['physical_crosspackage_complete'] is False
    assert cross['cross_package_yukawa_ratios'] is None
    assert cross['prior_portal_evidence']['status']=='prior_report_not_available'
    assert channels['native_physical_source_certified_models']==0
    assert channels['model_count']==16
    assert all(v['conditional_native_CG_channels'] is None for v in channels['models'])


def test_sqrt2_internal_registry_is_not_action_factor():
    # Hypothetical RGBeta Weyl action: -1/2 * (l,F) permutations *
    # internal sqrt(2) / complex-to-real scalar sqrt(2) = -1.
    internal=sp.sqrt(2)
    physical_candidate=-sp.Rational(1,2)*2*internal/sp.sqrt(2)
    assert physical_candidate==-1
    assert physical_candidate != internal
    # This is algebra under explicit hypotheses, never native certification.


def _checkout_with_windows_line_endings(tmp_path):
    """A synthetic clean Git checkout whose raw CRLF blob != committed LF blob."""
    import subprocess
    def git(*args):
        result=subprocess.run(['git','-C',str(tmp_path),*args],
            capture_output=True, text=True, check=True)
        return result.stdout.strip()
    git('init','-q')
    git('config','user.name','B2 Fixture')
    git('config','user.email','b2@example.invalid')
    git('config','core.autocrlf','true')
    rel='Lagrangian/model/SU2Invariants.wl'
    target=tmp_path/rel
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(b'(* tested *)\nOne[] := 1;\n')
    git('add',rel)
    git('commit','-qm','synthetic canonical LF')
    head=git('rev-parse','HEAD')
    expected=git('rev-parse',f'HEAD:{rel}')
    # Reproduce a normal Windows checkout with changed raw bytes.
    target.write_bytes(target.read_bytes().replace(b'\n',b'\r\n'))
    assert git_blob(target.read_bytes()) != expected
    assert git('hash-object',f'--path={rel}',rel) == expected
    return rel,target,head,expected


def test_windows_crlf_git_canonical_source_is_accepted(tmp_path, monkeypatch):
    import Numerical.lo_certification.physical_normalisation.b2_gauge_coupling as b2
    rel,target,head,expected=_checkout_with_windows_line_endings(tmp_path)
    # Add a missing synthetic file so source_inventory only tests provenance;
    # no T3 WL source-expression assertions are possible in a tiny fixture.
    monkeypatch.setattr(b2, 'HEAD', head)
    monkeypatch.setattr(b2, 'SOURCE_BLOBS', {rel:expected,
                        'a_deliberately_missing_wl_source.wl':'0'*40})
    status=b2.source_inventory(tmp_path)
    item=status['source_files'][rel]
    assert item['state']=='reference_match'
    assert item['verification']=='path_specific_git_clean_filter'
    assert item['git_canonical_matches_reference'] is True
    assert item['git_path_is_clean_against_HEAD'] is True
    assert status['source_reference_matches'] is False  # missing second file
    assert status['current_git_head']==head
    assert status['current_git_head_measured'] is True
    assert item['raw_worktree_git_blob_sha'] != expected
    assert item['git_canonical_blob_sha'] == expected


def test_windows_crlf_real_drift_and_unknown_commit_are_blocked(tmp_path, monkeypatch):
    import Numerical.lo_certification.physical_normalisation.b2_gauge_coupling as b2
    rel,target,head,expected=_checkout_with_windows_line_endings(tmp_path)
    monkeypatch.setattr(b2,'SOURCE_BLOBS',{rel:expected,'missing.wl':'0'*40})
    monkeypatch.setattr(b2,'HEAD','a'*40)
    status=b2.source_inventory(tmp_path)['source_files'][rel]
    assert status['state']=='drift_or_unreviewed'  # HEAD not known
    assert status['git_canonical_matches_reference'] is False
    monkeypatch.setattr(b2,'HEAD',head)
    target.write_bytes(target.read_bytes().replace(b'One[] := 1;',b'One[] := 2;'))
    status=b2.source_inventory(tmp_path)['source_files'][rel]
    assert status['state']=='drift_or_unreviewed'
    assert status['git_canonical_matches_reference'] is False
    assert status['git_path_is_clean_against_HEAD'] is False
    assert status['git_canonical_blob_sha'] != expected


def test_untracked_file_cannot_impersonate_reference_git_blob(tmp_path, monkeypatch):
    import Numerical.lo_certification.physical_normalisation.b2_gauge_coupling as b2
    rel,target,head,expected=_checkout_with_windows_line_endings(tmp_path)
    target.unlink()
    target.write_bytes(b'changed\r\n')
    monkeypatch.setattr(b2,'HEAD',head)
    monkeypatch.setattr(b2,'SOURCE_BLOBS',{rel:expected,'missing.wl':'0'*40})
    status=b2.source_inventory(tmp_path)['source_files'][rel]
    assert status['state']=='drift_or_unreviewed'
    assert status['git_canonical_matches_reference'] is False
