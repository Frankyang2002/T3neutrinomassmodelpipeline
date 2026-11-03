"""Adversarial tests for B2's native internal-vs-physical-action boundary."""
from __future__ import annotations

import copy
from pathlib import Path
import sys
import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'Numerical/lo_certification/physical_normalisation'))
import b2_rgbeta_internal_action as B
from b2_gauge_coupling import AuditBlocked, CLASSES, SOURCE_BLOBS


def fixture():
    families=[]
    for family,(d1,d2,df,alphas) in CLASSES.items():
        alpha=-1 if family in ('B','C') else 0
        maj=alpha==-1 and df%2==1
        rows=[]
        for name in B.OPERATOR_ORDER:
            exp=B.EXPECTED_RATIOS[name]
            kind='yukawa' if name in ('y1','y2') else 'mass' if name=='MF' else 'quartic'
            rows.append({'coupling':name,'kind':kind,'status':'NativeRegistryRatioMeasured',
                         'registered_fields':B.expected_fields(family,name,maj),
                         'actual_expected_fields_match':True,
                         'chirality':'not_applicable' if kind=='quartic' else 'Right',
                         'self_conjugate':(name!='lambdaT3') if kind=='quartic' else None,
                         'native_exact_ratio_and_metadata_verified':True,
                         'physical_Lagrangian_coefficient_verified':False,
                         'expected_internal_ratio_inputform':exp,
                         'sample_Coupling_over_Matrix_inputform':exp,
                         'sample_CouplingBar_over_Matrix_inputform':exp,
                         'sample_Coupling_inputform':(
                             f'{exp}*{name}' if kind=='quartic' else
                             f'{exp}*Matrix[Bar[{name}]][b4p,b4r]'),
                         'sample_CouplingBar_inputform':(
                             f'{exp}*Bar[{name}]' if name=='lambdaT3' else
                             f'{exp}*{name}' if kind=='quartic' else
                             f'{exp}*Matrix[{name}][b4p,b4r]'),
                         'registered_invariant_inputform':'eps[SU2L[fund], b4a, b4b]',
                         'unique_arrangements_count':1})
        families.append({'family':family,'dimensions_alpha':[d1,d2,df,alpha],
                         'majorana':maj,'status':'AllSixRegistryRatiosMeasured',
                         'operators':rows,'physical_action_certified':False})
    return {'format':B.FORMAT,
            'status':'NativeInternalTensorPrefactorsCapturedPhysicalActionPending',
            'internal_prefactor_is_physical_action_factor':False,
            'numeric_component_axes_verified':False,
            'physical_source_certified':False,
            'git_head_at_capture':'a'*40,
            'git_source_blobs_at_capture':{r:SOURCE_BLOBS[r] for r in B.ALLOWED_B2_PATHS},
            'classes':families}


def change(which, model='A', coupling='y1'):
    x=fixture()
    row=next(v for v in x['classes'] if v['family']==model)
    op=next(v for v in row['operators'] if v['coupling']==coupling)
    which(x,row,op)
    return x


def test_valid_all_five_and_sixteen():
    x=B.validate_native_data(fixture())
    assert x['families_checked']==5 and x['native_registry_ratios_measured']==30
    assert len(B.expand_all16(x))==16
    assert not x['physical_source_certified']
    assert all(not r['physical_source_certified'] for r in B.expand_all16(x))


def test_internal_factors_not_promoted_to_action():
    x=B.validate_native_data(fixture())
    for cls in x['families'].values():
        for op in cls['operator_internal_conventions']:
            assert op['physical_action_factor_from_this_ratio'] is None
            assert not op['physical_action_verified']
    assert x['physical_coupling_ratios_deducible_from_registry_only'] is False


@pytest.mark.parametrize('coupling,bad_ratio', [
    ('y1','2'),('y2','1'),('MF','Sqrt[2]'),
    ('lambdaT3','12'),('lambdaH1','3'),('lambdaH2','1/2'),
])
def test_wrong_internal_factor_blocked(coupling,bad_ratio):
    x=change(lambda x,r,o:o.update(sample_Coupling_over_Matrix_inputform=bad_ratio),coupling=coupling)
    with pytest.raises(AuditBlocked):B.validate_native_data(x)


def test_wrong_conjugate_ratio_blocked():
    x=change(lambda x,r,o:o.update(sample_CouplingBar_over_Matrix_inputform='-6'),coupling='lambdaT3')
    with pytest.raises(AuditBlocked):B.validate_native_data(x)


def test_wrong_majorana_fermion_blocked():
    x=change(lambda x,r,o:o.update(registered_fields='{S2, l, Fc}'),model='B',coupling='y2')
    with pytest.raises(AuditBlocked):B.validate_native_data(x)


def test_wrong_dirac_mass_order_blocked():
    x=change(lambda x,r,o:o.update(registered_fields='{Fc, F}'),model='E',coupling='MF')
    with pytest.raises(AuditBlocked):B.validate_native_data(x)


def test_wrong_mix_hc_flag_blocked():
    x=change(lambda x,r,o:o.update(self_conjugate=True),coupling='lambdaT3')
    with pytest.raises(AuditBlocked):B.validate_native_data(x)


def test_wrong_yukawa_chirality_blocked():
    x=change(lambda x,r,o:o.update(chirality='Left'),coupling='y1')
    with pytest.raises(AuditBlocked):B.validate_native_data(x)


def test_wrong_portal_selfconjugacy_blocked():
    x=change(lambda x,r,o:o.update(self_conjugate=False),coupling='lambdaH2')
    with pytest.raises(AuditBlocked):B.validate_native_data(x)


def test_wrong_reference_factor_blocked():
    x=change(lambda x,r,o:o.update(expected_internal_ratio_inputform='1'),coupling='lambdaH1')
    with pytest.raises(AuditBlocked):B.validate_native_data(x)


def test_missing_explicit_native_sample_blocked():
    x=change(lambda x,r,o:o.update(sample_Coupling_inputform='$Failed'),coupling='y1')
    with pytest.raises(AuditBlocked):B.validate_native_data(x)


def test_false_source_certificate_blocked():
    x=fixture();x['physical_source_certified']=True
    with pytest.raises(AuditBlocked):B.validate_native_data(x)


def test_false_action_certificate_blocked():
    x=change(lambda x,r,o:o.update(physical_Lagrangian_coefficient_verified=True))
    with pytest.raises(AuditBlocked):B.validate_native_data(x)


def test_missing_family_blocked():
    x=fixture();x['classes'].pop()
    with pytest.raises(AuditBlocked):B.validate_native_data(x)


def test_duplicate_family_blocked():
    x=fixture();x['classes'][1]=copy.deepcopy(x['classes'][0])
    with pytest.raises(AuditBlocked):B.validate_native_data(x)


def test_snapshot_head_mismatch_blocked(monkeypatch,tmp_path):
    x=fixture()
    provenance={'repository_head_measured':'b'*40}
    with pytest.raises(AuditBlocked,match='different Git HEAD'):
        B.native_snapshot_guard(x,tmp_path,provenance)


def test_snapshot_source_drift_blocked(monkeypatch,tmp_path):
    x=fixture();x['git_source_blobs_at_capture'][B.ALLOWED_B2_PATHS[0]]='0'*40
    provenance={'repository_head_measured':'a'*40}
    with pytest.raises(AuditBlocked,match='stale/mismatched'):
        B.native_snapshot_guard(x,tmp_path,provenance)


def test_snapshot_rejects_missing_source_in_inventory(tmp_path):
    x=fixture();del x['git_source_blobs_at_capture'][B.ALLOWED_B2_PATHS[-1]]
    with pytest.raises(AuditBlocked,match='inventory'):
        B.native_snapshot_guard(x,tmp_path,{'repository_head_measured':'a'*40})


def test_snapshot_rejects_changed_live_clean_hash(monkeypatch,tmp_path):
    x=fixture()
    monkeypatch.setattr(B,'git_result',lambda root,*args:(0,'f'*40))
    with pytest.raises(AuditBlocked,match='does not match local'):
        B.native_snapshot_guard(x,tmp_path,{'repository_head_measured':'a'*40})


def test_output_refuses_overwrite(tmp_path):
    output=tmp_path/'audit';output.mkdir()
    with pytest.raises(AuditBlocked,match='refusing existing'):
        B.write_outputs(output,{'format':B.FORMAT})


def test_physical_action_is_never_identified_from_sqrt2():
    # Direct negative control for mistaking the internally stored real-scalar
    # tensor normalization for a physical Lagrangian coupling ratio.
    x=B.validate_native_data(fixture())
    for family in x['families'].values():
        assert [o['internal_tensor_prefactor'] for o in
                family['operator_internal_conventions']]==[
                    'Sqrt[2]','Sqrt[2]','2','6','6','6']
        assert not family['physical_action_certified']


@pytest.mark.parametrize('coupling', ('lambdaT3', 'lambdaH1', 'lambdaH2'))
def test_native_scalar_quartic_samples_accepted(coupling):
    # Actual RGBeta v1.2.0 InputForm from the user's native witness.
    data=fixture()
    for cls in data['classes']:
        op=next(r for r in cls['operators'] if r['coupling']==coupling)
        assert 'Matrix[' not in op['sample_Coupling_inputform']
        assert 'Matrix[' not in op['sample_CouplingBar_inputform']
    assert B.validate_native_data(data)['native_registry_ratios_measured']==30


@pytest.mark.parametrize('coupling,bad_original,bad_bar', [
    ('lambdaT3','6*Bar[lambdaT3]','6*lambdaT3'),
    ('lambdaT3','6*lambdaH1','6*Bar[lambdaT3]'),
    ('lambdaH1','6*lambdaH1','6*Bar[lambdaH1]'),
    ('lambdaH2','6*Matrix[lambdaH2][]','6*lambdaH2'),
    ('y1','Sqrt[2]*y1','Sqrt[2]*Matrix[y1][b4p,b4r]'),
    ('MF','2*Matrix[MF][b4p,b4r]','2*Matrix[Bar[MF]][b4p,b4r]'),
])
def test_native_sample_orientation_and_shape_mutations_blocked(
    coupling,bad_original,bad_bar
):
    data=change(lambda x,r,o:o.update(
        sample_Coupling_inputform=bad_original,
        sample_CouplingBar_inputform=bad_bar),coupling=coupling)
    with pytest.raises(AuditBlocked,match='unexpected native coupling sample'):
        B.validate_native_data(data)
