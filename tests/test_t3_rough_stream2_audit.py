"""Small local algebra/contract tests; do not substitute for native Windows run."""
from __future__ import annotations
from pathlib import Path
import json
import numpy as np
import pytest
from scripts.AuditT3RoughStream2 import (audit_model, slug, sha256,
                                          CONVENTION, NAMES, W_FORMAT)


def fixture(tmp_path):
    root=tmp_path/'A_m4';root.mkdir()
    canonical=tmp_path/'canonical_input.json'
    ids={'matching_context_id':'stream1-real-id','matched_seed_id':'stream1-seed',
         'lagrangian_source_id':'real-source-id','scalar_basis_id':'A_m4:canonical'}
    canonical.write_text(json.dumps(ids))
    sh=sha256(canonical)
    c0=np.zeros((45,45,12,12),complex)
    c0[36,36,4,4]=1.5e-8+2e-9j
    c1=c0.copy();c1[36,36,2,2]=5e-12/4
    z=np.zeros_like(c0)
    data={'C_psi2phi2_ms':c1,'C_psi2phi2_tree_mf':c0,
          'C5_EFT1_ms':np.diag([5e-12,0,0]).astype(complex),
          'C5_EFT1_mf':np.zeros((3,3),complex),
          'C_LLSS_ms':c1[np.ix_(range(36,42),range(36,42),range(4,12),range(4,12))],
          'C_LLSS_mf':c0[np.ix_(range(36,42),range(36,42),range(4,12),range(4,12))],
          'beta_C_at_MF_16pi2':z,'beta_C_at_MS_16pi2':z,
          'trajectory_mu_gev':np.array([1e5,1e3]),
          'trajectory_C5_norm_gev_inverse':np.array([0,5e-12]),
          'trajectory_full_norm_gev_inverse':np.array([np.linalg.norm(c0),np.linalg.norm(c1)])}
    for endpoint in ('MF','MS'):
        for name in NAMES:
            data[f'beta_C5_{name}_at_{endpoint}_16pi2']=np.zeros((3,3),complex)
    np.savez_compressed(root/'wilson_at_MS.npz',**data)
    meta={'format':W_FORMAT,'tensor_npz':'wilson_at_MS.npz','model':'A_-4',
          'tensor_npz_sha256':sha256(root/'wilson_at_MS.npz'),
          'canonical_wilson_convention':CONVENTION,'scheme':'MSbar','units':'GeV',
          'synthetic_fixture':False,'wilson_evolved_to_ms':True,
          'includes_eft1_weinberg':True,'legacy_direct_llss_to_c5_log_added':False,
          'solver':{'method':'DOP853','atol':1e-21,
                    'backends':['existing_repository_matrix_free_engine']},
          'adm_variant':'fonseca_eq22_six_term','strict_1loop_endpoint_at_ms':False,
          'mu_f_gev':1e5,'mu_s_gev':1e3,
          'lepton_indices':[[36,37],[38,39],[40,41]],
          'higgs_source_indices':[0,1,2,3],
          'heavy_source_indices':list(range(4,12)),
          'C5_generated_from_LLSS_only_boundary':True,
          'C5_EFT1_ms_norm_gev_inverse':5e-12,
          'LLSS_mf_norm_gev_inverse':float(np.linalg.norm(data['C_LLSS_mf'])),
          'non_LLSS_tree_mf_norm_gev_inverse':0,
          'input':{'file_sha256':{'canonical_input_json':
                                {'path':str(canonical),'sha256':sh}}},
          **ids}
    (root/'wilson.json').write_text(json.dumps(meta))
    return root,meta


def test_slug_and_endpoint(tmp_path):
    assert slug('A_-4')=='A_m4' and slug('D_+2')=='D_p2'
    path,_=fixture(tmp_path)
    report=audit_model(path,tmp_path,require_existing_engine=True)
    assert report['generated_from_pure_LLSS']
    assert report['C5_MS_norm_GeV_inv']==pytest.approx(5e-12)
    assert report['n_weyl']==45


def test_checksums_reject_tampering(tmp_path):
    path,_=fixture(tmp_path)
    with (path/'wilson_at_MS.npz').open('ab') as f:f.write(b'tampered')
    with pytest.raises(ValueError,match='SHA256'):
        audit_model(path,tmp_path)


def test_requires_source_provenance(tmp_path):
    path,meta=fixture(tmp_path)
    meta['matching_context_id']='different-seed'
    (path/'wilson.json').write_text(json.dumps(meta))
    with pytest.raises(ValueError,match='inherited unchanged'):
        audit_model(path,tmp_path)


def test_denies_non_native_backend_when_requested(tmp_path):
    path,meta=fixture(tmp_path)
    meta['solver']['backends']=['portable_formula_mirror_no_repository_engine']
    (path/'wilson.json').write_text(json.dumps(meta))
    with pytest.raises(ValueError,match='repository-native'):
        audit_model(path,tmp_path,require_existing_engine=True)
