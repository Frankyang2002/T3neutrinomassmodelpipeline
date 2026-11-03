"""Exact G2 pseudoreal-basis tests, including negative and promotion guards.

These tests do NOT establish the native RGBeta component convention.
"""
from __future__ import annotations

from itertools import product
from pathlib import Path
import sys

import pytest
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'Numerical/lo_certification/physical_normalisation'))
from b2_gauge_coupling import AuditBlocked
from b2_rgbeta_ward_witness import ANTI,FUND,EPS,ward_residual,repaired_22_mix
import b2_su2_intertwiner as G


def test_exact_pseudoreal_intertwiner():
    p=[sp.Matrix([[0,1],[1,0]])/2,
       sp.Matrix([[0,-sp.I],[sp.I,0]])/2,
       sp.diag(sp.Rational(1,2),-sp.Rational(1,2))]
    assert all(sp.simplify(EPS*t.conjugate()*EPS.inv()+t)==sp.zeros(2) for t in p)


def test_mix_literal_negative_and_alternative_zero():
    w=G.core_witness()
    assert w['literal_standard_BC_mix_ward_squared']==[3,0,3]
    assert w['alternative_formal_BC_mix_ward_squared']==[0,0,0]
    assert w['transported_standard_BC_mix_ward_squared']==[0,0,0]


def test_y1_literal_negative_and_alternative_zero():
    w=G.core_witness()
    assert w['literal_standard_B_y1_ward_squared']==[2,0,2]
    assert w['alternative_formal_B_y1_ward_squared']==[0,0,0]
    assert w['transported_standard_B_y1_ward_squared']==[0,0,0]


def test_transport_equals_candidate_all_16_components():
    assert all(sp.simplify(G.transported_mix(ix)-repaired_22_mix(ix))==0
               for ix in product(range(2),repeat=4))


def test_y1_transport_equals_negative_delta():
    assert all(G.transported_y1(ix)==-int(ix[0]==ix[1])
               for ix in product(range(2),repeat=2))


def test_transport_preserves_H_bose_symmetry():
    assert all(G.transported_mix((h,k,s,b))==G.transported_mix((k,h,s,b))
               for h,k,s,b in product(range(2),repeat=4))


def test_wrong_literal_mix_fails_standard_ward():
    assert ward_residual(G.literal_rgbeta_22_mix,(FUND,FUND,FUND,ANTI))!=(0,0,0)


def test_wrong_literal_y1_fails_standard_ward():
    assert ward_residual(G.literal_B_y1,(ANTI,FUND))!=(0,0,0)


def test_epsilon_transpose_changes_sign_not_ward():
    assert all(G.transported_mix(ix,EPS.T)==-G.transported_mix(ix)
               for ix in product(range(2),repeat=4))
    assert ward_residual(lambda ix:G.transported_mix(ix,EPS.T),
                         (FUND,FUND,FUND,ANTI))==(0,0,0)


def test_bad_intertwiner_blocks():
    with pytest.raises(AuditBlocked):G.core_witness(sp.zeros(2))
    with pytest.raises(AuditBlocked):G.core_witness(sp.eye(2))


def test_triplet_weight_metric_synthetic_matchete_axes():
    u=sp.Matrix([[-sp.sqrt(2)*sp.I/2,0,sp.sqrt(2)*sp.I/2],
                 [0,1,0],
                 [sp.sqrt(2)/2,0,sp.sqrt(2)/2]])
    result=G.triplet_witness({3:u})
    assert result['Matchete_triplet_weight_metric']==[
        ['0','0','1'],['0','1','0'],['1','0','0']]
    assert result['metric_equals_RGBeta_delS2_axes'] is False


def test_wrong_triplet_rotation_rejected():
    with pytest.raises(AuditBlocked):G.triplet_witness({3:sp.eye(3)})


def test_g1_and_formal_evidence_mandatory(tmp_path):
    with pytest.raises(AuditBlocked,match='missing required evidence'):
        G._json_file(tmp_path/'missing.json')


def test_reject_fake_physical_output_or_overwrite(tmp_path):
    out=tmp_path/'audit'
    item={'format':G.FORMAT,'physical_source_certified_models':0,
          'safe_to_integrate_into_production':False,
          'source_provenance':{'repository_head_measured':'a'*40},
          'g1_capture_reuse_provenance':{'capture_head':'a'*40,'same_head':True}}
    G.write_report(out,item)
    assert (out/'intertwiner_witness.json').is_file()
    with pytest.raises(AuditBlocked,match='refusing existing'):
        G.write_report(out,item)
    with pytest.raises(AuditBlocked,match='invalid G2'):
        G.write_report(tmp_path/'unsafe',dict(item,physical_source_certified_models=16))


def test_reject_incomplete_matchete_cg_evidence():
    with pytest.raises(AuditBlocked):G.evidence_guard({'cases':[]},{'status':'bad','count':16})


def test_physical_flags_stay_false_in_core():
    obj=G.core_witness()
    assert obj['epsilon_is_RGBeta_native_component_witness'] is False
    assert obj['physical_coupling_conversion'] is None


def test_all16_matchete_inventory_adversarial():
    from b2_gauge_coupling import models
    arr=[];cov=[]
    for m in models():
        rep={k:m[k] for k in ('dS1','dS2','dF','alpha')}
        arr.append({'model':rep,'vertices':{
            'T3Y1CG':{},'T3Y2CG':{},'T3MixCG':{}}})
        cov.append({'model':rep,'all_covariant':True,'vertices':{
            'T3Y1CG':{'covariant':True},
            'T3Y2CG':{'covariant':True},
            'T3MixCG':{'covariant':True}}})
    data={'cases':arr}
    consistency={'status':'cg_generator_covariance_verified',
                 'count':16,'models':cov}
    assert len(G.evidence_guard(data,consistency))==16
    cov[-1]['vertices']['T3MixCG']['covariant']=False
    with pytest.raises(AuditBlocked):G.evidence_guard(data,consistency)
    cov[-1]['vertices']['T3MixCG']['covariant']=True
    arr[-1]['model']=dict(arr[0]['model'])
    with pytest.raises(AuditBlocked):G.evidence_guard(data,consistency)


def test_six_bc_synthetic_mix_comparisons(monkeypatch):
    from b2_gauge_coupling import models
    from b2_rgbeta_ward_witness import literal_rgbeta_22_mix
    tensor=sp.MutableDenseNDimArray(
        [2*repaired_22_mix(ix) for ix in product(range(2),repeat=4)],
        (2,2,2,2))
    monkeypatch.setattr(G,'native_tensor',lambda *args:tensor)

    def independent_fit(ref,obs,n):
        num=sum(sp.conjugate(ref(ix))*obs(ix)
                for ix in product(range(2),repeat=n))
        den=sum(sp.conjugate(ref(ix))*ref(ix)
                for ix in product(range(2),repeat=n))
        ratio=sp.simplify(num/den)
        error=sp.simplify(sum(
            sp.conjugate(obs(ix)-ratio*ref(ix))*(obs(ix)-ratio*ref(ix))
            for ix in product(range(2),repeat=n)))
        return {'complex_ratio_to_chosen_candidate':str(ratio),
                'exact_residual_squared_norm':str(error),
                'is_proportional_exact':bool(error==0)}
    monkeypatch.setattr(G,'fit_tensor',independent_fit)
    cases={(m['dS1'],m['dS2'],m['dF'],m['alpha']):
           {'vertices':{'T3MixCG':{}}} for m in models()}
    rows=G.compare_b_c_mixing(cases,{})
    assert len(rows)==6
    assert all(r['Matchete_over_transported_candidate']=='2' for r in rows)
    assert all(r['actual_RGBeta_to_Matchete_coupling_ratio'] is None for r in rows)

    # Adversarial incompatible tensor: refusing a false conditional fit.
    bad=sp.MutableDenseNDimArray(
        [literal_rgbeta_22_mix(ix) for ix in product(range(2),repeat=4)],
        (2,2,2,2))
    monkeypatch.setattr(G,'native_tensor',lambda *args:bad)
    with pytest.raises(AuditBlocked,match='disagrees'):
        G.compare_b_c_mixing(cases,{})


def _fixture_committed_physics_sources(tmp_path, monkeypatch):
    """Real Git ancestry and six-file source guards, independent of user repo."""
    import subprocess
    from b2_gauge_coupling import SOURCE_BLOBS as ORIGINAL_SOURCE_BLOBS

    root = tmp_path / 'repo'
    root.mkdir()

    def git(*args):
        proc = subprocess.run(['git', '-C', str(root), *args],
                              capture_output=True, text=True, check=True)
        return proc.stdout.strip()

    git('init', '-q')
    git('config', 'core.autocrlf', 'false')
    pinned = {}
    for rel in ORIGINAL_SOURCE_BLOBS:
        file = root / rel
        file.parent.mkdir(parents=True, exist_ok=True)
        file.write_bytes((f'(* fixture source: {rel} *)\n').encode('utf-8'))
        git('add', '--', rel)
        pinned[rel] = git('hash-object', f'--path={rel}', '--', rel)
    git('-c', 'user.name=Fixture', '-c', 'user.email=fixture@example.invalid',
        'commit', '-qm', 'initial fixture source')
    captured = git('rev-parse', 'HEAD')
    monkeypatch.setattr(G, 'SOURCE_BLOBS', pinned)
    native = {
        'git_head_at_capture': captured,
        'git_source_blobs_at_capture': {
            rel: pinned[rel] for rel in G.ALLOWED_B2_PATHS
        },
    }
    return root, git, captured, pinned, native


def test_g1_reuse_accepts_real_descendant_unchanged_six_sources(tmp_path, monkeypatch):
    root, git, capture, pinned, native = _fixture_committed_physics_sources(tmp_path, monkeypatch)
    (root / 'UNRELATED').write_text('unrelated new work\n', encoding='utf-8')
    git('add', 'UNRELATED')
    git('-c', 'user.name=Fixture', '-c', 'user.email=fixture@example.invalid',
        'commit', '-qm', 'unrelated advancement')
    head = git('rev-parse', 'HEAD')
    assert capture != head
    result = G.g1_capture_reuse_guard(native, root, {
        'repository_head_measured': head,
        'all_six_source_blobs_unchanged': True,
    })
    assert result['capture_head'] == capture
    assert result['current_head'] == head
    assert result['capture_is_ancestor'] is True
    assert result['same_head'] is False
    assert result['all_six_source_blobs_verified_in_capture_tree'] is True
    assert result['g1_native_reran_at_current_head'] is False
    assert result['physical_source_certified'] is False


def test_g1_reuse_same_head_delegates_original_guard(tmp_path, monkeypatch):
    root, git, capture, pinned, native = _fixture_committed_physics_sources(tmp_path, monkeypatch)
    called = []
    monkeypatch.setattr(G, 'native_snapshot_guard',
                        lambda report, repo, source: called.append((report, repo, source)))
    provenance = {'repository_head_measured': capture,
                  'all_six_source_blobs_unchanged': True}
    result = G.g1_capture_reuse_guard(native, root, provenance)
    assert len(called) == 1
    assert called[0] == (native, root, provenance)
    assert result['same_head'] is True


def test_g1_reuse_rejects_nonancestor(tmp_path, monkeypatch):
    root, git, capture, pinned, native = _fixture_committed_physics_sources(tmp_path, monkeypatch)
    git('checkout', '--orphan', 'unrelated')
    for rel in pinned:
        (root / rel).write_bytes((f'(* fixture source: {rel} *)\n').encode('utf-8'))
    git('add', '--all')
    git('-c', 'user.name=Fixture', '-c', 'user.email=fixture@example.invalid',
        'commit', '-qm', 'other unrelated history')
    with pytest.raises(AuditBlocked, match='not a verified ancestor'):
        G.g1_capture_reuse_guard(native, root, {
            'repository_head_measured': git('rev-parse', 'HEAD'),
            'all_six_source_blobs_unchanged': True})


def test_g1_reuse_rejects_changed_capture_tree(tmp_path, monkeypatch):
    root, git, capture, pinned, native = _fixture_committed_physics_sources(tmp_path, monkeypatch)
    # Adversarial expected hash: the real capture tree cannot be rewritten.
    rel = next(iter(pinned))
    pinned[rel] = 'f' * 40
    (root/'UNRELATED').write_text('new\n', encoding='utf-8')
    git('add', 'UNRELATED')
    git('-c', 'user.name=Fixture', '-c', 'user.email=fixture@example.invalid',
        'commit', '-qm', 'unrelated')
    with pytest.raises(AuditBlocked, match='different/missing pinned source'):
        G.g1_capture_reuse_guard(native, root, {
            'repository_head_measured': git('rev-parse','HEAD'),
            'all_six_source_blobs_unchanged': True})


def test_g1_reuse_rejects_dirty_worktree(tmp_path, monkeypatch):
    root, git, capture, pinned, native = _fixture_committed_physics_sources(tmp_path, monkeypatch)
    (root/'UNRELATED').write_text('new\n', encoding='utf-8')
    git('add', 'UNRELATED')
    git('-c', 'user.name=Fixture', '-c', 'user.email=fixture@example.invalid',
        'commit', '-qm', 'unrelated')
    rel = next(iter(pinned))
    (root/rel).write_text('unexpected source drift\n',encoding='utf-8')
    with pytest.raises(AuditBlocked, match='source drift'):
        G.g1_capture_reuse_guard(native, root, {
            'repository_head_measured': git('rev-parse','HEAD'),
            'all_six_source_blobs_unchanged': True})


def test_g1_reuse_rejects_falsified_capture_inventory(tmp_path, monkeypatch):
    root, git, capture, pinned, native = _fixture_committed_physics_sources(tmp_path, monkeypatch)
    (root/'UNRELATED').write_text('new\n', encoding='utf-8')
    git('add', 'UNRELATED')
    git('-c', 'user.name=Fixture', '-c', 'user.email=fixture@example.invalid',
        'commit', '-qm', 'unrelated')
    rel = G.ALLOWED_B2_PATHS[0]
    native['git_source_blobs_at_capture'][rel] = 'e' * 40
    with pytest.raises(AuditBlocked, match='native recorded different source blob'):
        G.g1_capture_reuse_guard(native, root, {
            'repository_head_measured': git('rev-parse','HEAD'),
            'all_six_source_blobs_unchanged': True})


def test_g1_reuse_requires_six_source_current_verification(tmp_path, monkeypatch):
    root, git, capture, pinned, native = _fixture_committed_physics_sources(tmp_path, monkeypatch)
    with pytest.raises(AuditBlocked, match='current six-source provenance'):
        G.g1_capture_reuse_guard(native, root, {
            'repository_head_measured': capture,
            'all_six_source_blobs_unchanged': False})


def test_core_witness_uses_json_primitive_ward_values():
    import json
    core = G.core_witness()
    encoded = json.dumps(core)
    loaded = json.loads(encoded)
    assert loaded['literal_standard_BC_mix_ward_squared'] == [3,0,3]
    assert loaded['transported_standard_BC_mix_ward_squared'] == [0,0,0]
    assert loaded['epsilon_is_RGBeta_native_component_witness'] is False
