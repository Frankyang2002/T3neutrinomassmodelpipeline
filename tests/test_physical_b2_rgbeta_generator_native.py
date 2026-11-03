"""B2 G2 native symbolic generator probe: adversarial isolated Python tests.

These tests mock the native engine. Passing them does NOT establish numerical
RGBeta components, physical action factors, or physical certifications.
"""
from __future__ import annotations

import copy
import importlib.util
from pathlib import Path
import json
import subprocess
import sys
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'Numerical/lo_certification/physical_normalisation'))
import b2_rgbeta_generator_native as B
from b2_gauge_coupling import AuditBlocked, SOURCE_BLOBS, CLASSES


def proven():
    return {
        'repository_head_measured': 'a'*40,
        'all_six_source_blobs_unchanged': True,
        'validated_git_blobs': [
            {'path':p,'git_blob_sha':v,'validated_with_git_clean_filter':True}
            for p,v in SOURCE_BLOBS.items()
        ],
    }


def gauge(d):
    return ('{U1Y[0]}' if d==1 else
            '{U1Y[0], SU2L[fund]}' if d==2 else
            '{U1Y[0], SU2L[S2]}')


def native():
    rows=[]
    for name,(d1,d2,df,alphas) in CLASSES.items():
        a=-1 if name in ('B','C') else 0
        rows.append({
            'family':name,'dimensions_alpha':[d1,d2,df,a],
            'registered':True,'scalar_representations_match_registration':True,
            'fermion_representations_match_registration':True,
            'all_scalar_field_and_conjugate_indices_present':True,
            'yukawa_field_orders_match':True,
            'S1_gauge_reps':gauge(d1),
            'S2_gauge_reps':gauge(d2),
            'F_gauge_reps':gauge(df),
            'field_index_ids':{'S1':10,'Bar[S1]':11,'S2':12,'Bar[S2]':13},
            'native_package_field_index_map_available':True,
            'physical_action_certified':False,'numeric_SU2_axes_certified':False,
        })
    return {
        'format':B.FORMAT,'status':B.STATUS,
        'git_head_at_capture':'a'*40,
        'git_source_blobs_at_capture':{
            p:SOURCE_BLOBS[p] for p in B.ALLOWED_B2_PATHS},
        'git_provenance_capture_method':'Python six-source clean-filter guard',
        'rgbeta_package_file':'C:/Users/x/RGBeta/Kernel/init.m',
        'field_index_map_symbol_names':[
            'RGBeta`FieldsAndCouplings`PackagePrivate`$fieldIndexMap'],
        'field_index_map_unique_package_scope':True,
        'Tscal_symbol_names':['RGBeta`FieldsAndCouplings`PackagePrivate`Tscal'],
        'Tscal_downvalues_inputform':'{HoldPattern[Tscal[A_,a_,b_]] :> '
            'Module[{$fieldIndexMap,tGen,Bar[scal]}, tGen[rep,A,a,b] + Bar[scal]]}',
        'Tscal_rule_has_conjugate_scalar_blocks':True,
        'tGen_bar_transpose_identity_native':True,
        'eps_sample_remains_symbolic':True,
        'delS2_sample_remains_symbolic':True,
        'classes':rows,
        'RGBeta_numeric_eps_components_certified':False,
        'RGBeta_numeric_delS2_components_certified':False,
        'RGBeta_to_Matchete_basis_identified':False,
        'physical_coupling_dictionary_certified':False,
        'physical_source_certified_models':0,
        'safe_to_integrate_into_production':False,
        'EFT_performed':False,
    }


def test_five_native_families_validate_without_physical_promotion():
    v=B.validate_data(native(),proven())
    assert v['families_registered']==5
    assert v['numeric_SU2_component_axes_identified'] is False
    assert v['B_C_index_orientation_resolved_from_native_numeric_axes'] is False
    assert v['physical_source_certified_models']==0


@pytest.mark.parametrize('field', [
    'physical_source_certified_models','physical_coupling_dictionary_certified',
    'RGBeta_to_Matchete_basis_identified','RGBeta_numeric_eps_components_certified',
    'RGBeta_numeric_delS2_components_certified',
    'safe_to_integrate_into_production','EFT_performed'])
def test_reject_false_physical_claim(field):
    n=native();n[field]=1 if field=='physical_source_certified_models' else True
    with pytest.raises(AuditBlocked):B.validate_data(n,proven())


@pytest.mark.parametrize('field,value', [
    ('tGen_bar_transpose_identity_native',False),
    ('Tscal_rule_has_conjugate_scalar_blocks',False),
    ('eps_sample_remains_symbolic',False),
    ('delS2_sample_remains_symbolic',False),
    ('Tscal_symbol_names',[]),
    ('Tscal_downvalues_inputform','unresolved'),
    ('git_head_at_capture','b'*40),
    ('status','BLOCKED_NATIVE_SYMBOLIC_GENERATOR_WITNESS'),
])
def test_reject_incomplete_native_measurement(field,value):
    n=native();n[field]=value
    with pytest.raises(AuditBlocked):B.validate_data(n,proven())


def test_missing_or_duplicate_models_rejected():
    n=native();n['classes'].pop()
    with pytest.raises(AuditBlocked):B.validate_data(n,proven())
    n=native();n['classes'][1]=copy.deepcopy(n['classes'][0])
    with pytest.raises(AuditBlocked):B.validate_data(n,proven())


def test_source_inventory_mismatch_rejected():
    n=native();n['git_source_blobs_at_capture'][B.ALLOWED_B2_PATHS[0]]='b'*40
    with pytest.raises(AuditBlocked):B.validate_data(n,proven())


def test_unbarred_and_barred_field_index_negative_controls():
    n=native();n['classes'][1]['field_index_ids']['Bar[S1]']=10
    with pytest.raises(AuditBlocked):B.validate_data(n,proven())
    n=native();n['classes'][1]['field_index_ids']['Bar[S1]']=-1
    with pytest.raises(AuditBlocked):B.validate_data(n,proven())


def test_field_registered_SU2_dimensions_negative_control():
    n=native();n['classes'][1]['S1_gauge_reps']=gauge(3)
    with pytest.raises(AuditBlocked):B.validate_data(n,proven())
    n=native();n['classes'][0]['S1_gauge_reps']=gauge(2)
    with pytest.raises(AuditBlocked):B.validate_data(n,proven())


def test_source_guard_must_cover_six_verified_paths():
    p=proven();p['validated_git_blobs'].pop()
    with pytest.raises(AuditBlocked):B.child_environment(p,{})
    p=proven();p['validated_git_blobs'][0]['validated_with_git_clean_filter']=False
    with pytest.raises(AuditBlocked):B.child_environment(p,{})
    p=proven();p['all_six_source_blobs_unchanged']=False
    with pytest.raises(AuditBlocked):B.child_environment(p,{})


def test_env_discards_inherited_injected_provenance():
    parent={x:'f'*40 for x in ('T3_B2_PROVENANCE_HEAD',*B.ENV_NAMES.values())}
    before=parent.copy();env=B.child_environment(proven(),parent)
    assert parent==before
    assert env['T3_B2_PROVENANCE_HEAD']=='a'*40
    assert all(env[v]==SOURCE_BLOBS[k] for k,v in B.ENV_NAMES.items())


def test_restricted_flag_blocks_native_probe(tmp_path):
    root=tmp_path/'root';(root/B.SCRIPT).parent.mkdir(parents=True)
    (root/B.SCRIPT).write_text('(* fixture only *)')
    out=root/B.OUTPUT_DEFAULT
    with pytest.raises(AuditBlocked,match='non-restricted'):
        B.run_native(root,out,parent_env={'T3_EXPERIMENTAL_GENERAL_FERMIONS':'1'})


def test_existing_output_rejected_before_native(tmp_path):
    root=tmp_path/'root';(root/B.SCRIPT).parent.mkdir(parents=True)
    (root/B.SCRIPT).write_text('(* fixture only *)')
    out=root/B.OUTPUT_DEFAULT;out.parent.mkdir(parents=True)
    out.write_text('{}')
    with pytest.raises(AuditBlocked,match='existing'):
        B.run_native(root,out,parent_env={})


def test_runner_provenance_and_native_validator_with_mocked_wolfram(tmp_path,monkeypatch):
    root=tmp_path/'root';(root/B.SCRIPT).parent.mkdir(parents=True)
    (root/B.SCRIPT).write_text('(* fixture only *)')
    out=root/B.OUTPUT_DEFAULT
    monkeypatch.setattr(B,'source_provenance',lambda r:proven())
    seen={}
    def fake_run(args,**kwargs):
        seen.update({'args':args,**kwargs})
        out.parent.mkdir(parents=True)
        out.write_text(json.dumps(native()))
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(B.subprocess,'run',fake_run)
    v=B.run_native(root,out,wolfram='fake',parent_env={'T3_B2_PROVENANCE_HEAD':'f'*40})
    assert seen['args']==['fake','-file',str(root/B.SCRIPT),str(out)]
    assert seen['env']['T3_B2_PROVENANCE_HEAD']=='a'*40
    assert v['physical_source_certified_models']==0
    assert len(v['native_output_sha256'])==64


def test_native_nonzero_exit_is_blocking(tmp_path,monkeypatch):
    root=tmp_path/'root';(root/B.SCRIPT).parent.mkdir(parents=True)
    (root/B.SCRIPT).write_text('(* fixture only *)')
    monkeypatch.setattr(B,'source_provenance',lambda r:proven())
    monkeypatch.setattr(B.subprocess,'run',lambda *a,**k:SimpleNamespace(returncode=7))
    with pytest.raises(AuditBlocked,match='exited 7'):
        B.run_native(root,root/B.OUTPUT_DEFAULT,parent_env={})


def test_wolfram_source_checks_guard_against_git_runprocess():
    script=(ROOT/B.SCRIPT).read_text(encoding='utf-8')
    assert 'RunProcess[' not in script
    assert 'Needs["RGBeta`"]' in script
    assert 'DownValues[symbol]' in script
    assert 'tGen[Bar[SU2L[fund]]' in script
    for key in ('T3_B2_PROVENANCE_HEAD',*B.ENV_NAMES.values()):
        assert f'Environment["{key}"]' in script
    assert '"physical_source_certified_models"->0' in script


def test_package_script_import_and_entrypoint():
    p=ROOT/'scripts/RunPhysicalB2RGBetaGeneratorNative.py'
    spec=importlib.util.spec_from_file_location('B2G2Runner',p)
    assert spec and spec.loader
    m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    assert m.ROOT==ROOT
    assert callable(m.main)


@pytest.mark.parametrize('bad_names', [[], ['Global`$fieldIndexMap'],
    ['RGBeta`FieldsAndCouplings`PackagePrivate`$fieldIndexMap',
     'RGBeta`Other`$fieldIndexMap']])
def test_unresolved_or_ambiguous_private_field_index_symbol_rejected(bad_names):
    n=native();n['field_index_map_symbol_names']=bad_names
    with pytest.raises(AuditBlocked,match='uniquely resolved'):
        B.validate_data(n,proven())


def test_unmeasured_private_field_index_map_rejected():
    n=native();n['field_index_map_unique_package_scope']=False
    with pytest.raises(AuditBlocked,match='uniquely resolved'):
        B.validate_data(n,proven())
    n=native();n['classes'][2]['native_package_field_index_map_available']=False
    with pytest.raises(AuditBlocked,match='field-index registration failed'):
        B.validate_data(n,proven())


def test_wolfram_reads_native_private_field_index_registry_not_global():
    script=(ROOT/B.SCRIPT).read_text(encoding='utf-8')
    assert 'Names["RGBeta`*`$fieldIndexMap"]' in script
    assert 'ToExpression[First[names]]' in script
    assert 'indexMap = B6NativeIndexMap[];' in script
    assert 'Lookup[indexMap,"All",<||>]' in script
    assert 'Lookup[$fieldIndexMap' not in script
    assert '"native_package_field_index_map_available"' in script
    assert '"field_index_map_symbol_names"' in script
