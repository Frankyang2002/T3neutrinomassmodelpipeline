"""Contract tests for the continuous RG -> C5 -> low-energy adapter."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from scripts import RunT3FullScalarNeutrino as m


def test_classification_is_exhaustive():
    assert len(m.MODELS)==16
    assert len(m.MIXED)==2
    assert len(m.SELF)==6
    assert len(set(m.MODELS)-m.MIXED-m.SELF)==8
    assert not (m.MIXED & m.SELF)


def test_complex_json_strict():
    assert m.complex_pair({'real':2,'imag':-1},'mu')==2-1j
    with pytest.raises(ValueError,match='complex real/imag'):
        m.complex_pair({'real':2},'mu')
    with pytest.raises(ValueError,match='nonfinite'):
        m.complex_pair({'real':float('inf'),'imag':0},'mu')


def test_majorana_c5_is_symmetric():
    c=np.eye(3,dtype=complex)*1e-14
    assert np.array_equal(m.validated_c5(c),c)
    c[0,1]=1e-13
    with pytest.raises(ValueError,match='symmetric'):
        m.validated_c5(c)
    with pytest.raises(ValueError,match='3x3'):
        m.validated_c5(np.eye(2))


def test_complex_matrix_roundtrip():
    original=np.diag([1+2j,3-4j,5+0j])
    assert np.array_equal(m.complex_matrix(m.matrix_payload(original),'matrix'),original)


def test_upstream_requires_uncertified_continuous_running(tmp_path):
    source=tmp_path/'source.json'
    row={'format':m.full.FORMAT,'status':m.full.STATUS,'model':'B_-1',
         'physical_matchete_rgbeta_mapping_certified':False,
         'full_sequential_neutrino_mass_prediction':False,
         'fermion_threshold_wilson_matching_performed':False,
         'scalar_threshold_wilson_matching_performed':False,
         'scalar_endpoint_physically_admissible':True,
         'scalar_endpoint_spectrum':{'tachyonic':False},
         'quadratic_kind':'mixed',
         'settings':{'scalar_endpoint_gev':1e3,'fermion_threshold_gev':m.full.s30.STOP}}
    source.write_text(json.dumps(row),encoding='utf-8')
    assert m.verify_continuous('B_-1',source)['model']=='B_-1'
    row['full_sequential_neutrino_mass_prediction']=True
    source.write_text(json.dumps(row),encoding='utf-8')
    with pytest.raises(ValueError,match='certification changed'):
        m.verify_continuous('B_-1',source)


def _fixture(tmp_path,monkeypatch,model):
    monkeypatch.setattr(m,'ROOT',tmp_path)
    folder=tmp_path/'running'/m.full.model_tag(model)
    folder.mkdir(parents=True)
    (folder/'continuous_trajectory.json').write_text('{}')
    (folder/'eft1.json').write_text('{}')
    matching=tmp_path/'matching.txt'
    matching.write_text('source expression')
    monkeypatch.setattr(m,'verify_continuous',lambda model,path:{
        'saved_candidate_eft1_endpoint':'eft1.json',
        'uv_fermion_threshold_state':{},'eft1_scalar_endpoint_state':{},
        'quadratic_kind':('mixed' if model in m.MIXED else
                          'self_s1' if model in m.SELF else 'none')})
    monkeypatch.setattr(m,'decode_uv',lambda obj:'uv')
    monkeypatch.setattr(m,'decode_eft',lambda obj:'eft')
    monkeypatch.setattr(m,'low_energy',lambda model,c5,eft,**kw:{
        'takagi_masses_ev_descending':[0.05,0.01,0.001],
        'oscillation_observables_normal_ordering':{'sin2_theta13':0.02},
        'sm_weinberg_rge_performed':True})
    return folder,matching


def test_mixed_uses_actual_mass_basis_matching(tmp_path,monkeypatch):
    folder,matching=_fixture(tmp_path,monkeypatch,'B_-1')
    called=[]
    def match(model,uv,endpoint_dir):
        called.append((model,uv,endpoint_dir))
        return np.eye(3)*1e-12,matching,'mass-eigenstate-one-scale-direct-A5',{}
    monkeypatch.setattr(m,'matched_bc_c5',match)
    monkeypatch.setattr(m,'matched_nonmixed_c5',lambda *args:pytest.fail('wrong kernel'))
    data=m.compute_one('B_-1',source_root=folder.parent,destination_root=tmp_path/'dest',resume=False)
    assert called==[('B_-1','uv',folder)]
    assert data['nonzero_bilinear_in_matching'] is True
    assert data['additional_quadratic_c5_omitted'] is False
    assert data['uv_evolved_yukawas_used'] is True
    assert data['full_general_scalar_physical_c5_certified'] is False


def test_self_bilinear_is_never_mislabeled_as_full_matching(tmp_path,monkeypatch):
    folder,matching=_fixture(tmp_path,monkeypatch,'A_+0')
    monkeypatch.setattr(m,'matched_nonmixed_c5',lambda *args:(np.eye(3)*1e-12,matching,'reduced-hierarchical'))
    monkeypatch.setattr(m,'matched_bc_c5',lambda *args:pytest.fail('wrong kernel'))
    data=m.compute_one('A_+0',source_root=folder.parent,destination_root=tmp_path/'dest',resume=False)
    assert data['status']=='RestrictedKernelControlSelfBilinearNotMatched'
    assert data['additional_quadratic_c5_omitted'] is True
    assert data['nonzero_bilinear_in_matching'] is False
    assert data['full_sequential_threshold_matching_performed'] is False


def test_unmixed_model_uses_existing_original_c5_kernel(tmp_path,monkeypatch):
    folder,matching=_fixture(tmp_path,monkeypatch,'B_+1')
    monkeypatch.setattr(m,'matched_nonmixed_c5',lambda *args:(np.eye(3)*1e-12,matching,'hierarchical'))
    data=m.compute_one('B_+1',source_root=folder.parent,destination_root=tmp_path/'dest',resume=False)
    assert data['status']=='CandidateNoAdditionalQuadraticMatchedC5'
    assert data['additional_quadratic_c5_omitted'] is False
    assert data['nonzero_bilinear_in_matching'] is False


def test_resume_requires_same_input_hash(tmp_path,monkeypatch):
    folder,matching=_fixture(tmp_path,monkeypatch,'B_+1')
    monkeypatch.setattr(m,'matched_nonmixed_c5',lambda *args:(np.eye(3)*1e-12,matching,'hierarchical'))
    target=tmp_path/'dest'
    initial=m.compute_one('B_+1',source_root=folder.parent,destination_root=target,resume=False)
    assert m.compute_one('B_+1',source_root=folder.parent,destination_root=target,resume=True)==initial
    with pytest.raises(ValueError,match='output exists'):
        m.compute_one('B_+1',source_root=folder.parent,destination_root=target,resume=False)
    (folder/'continuous_trajectory.json').write_text('{"modified":true}')
    with pytest.raises(ValueError,match='sources changed'):
        m.compute_one('B_+1',source_root=folder.parent,destination_root=target,resume=True)
