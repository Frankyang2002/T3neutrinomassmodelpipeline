"""Independent synthetic-only adapter contract regression (no physics outputs).

Runs without the C scalar matcher, SMEFT ODE or external Matchete data. It
specifically guards against losing the model, units, scheme and ordered real
scalar labels when writing C-owned normalized endpoint manifests.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from Numerical.parallel_t3.threshold.models import ALL_MODELS, n_real, scalar_labels
from Numerical.parallel_t3.threshold.pipeline import B_FORMAT, W_FORMAT, Q_CONVENTION, W_CONVENTION
from Numerical.parallel_t3.threshold.rough_adapter import AdapterBlocked, prepare_pair


def _producer_pair(root: Path, model: str, *, declare_labels: bool = False):
    n=n_real(model)
    folder=root/model.replace('_-','_m').replace('_+','_p')
    folder.mkdir(parents=True,exist_ok=True)
    sid=f'synthetic:manifest-roundtrip:{model}'
    common={'model':model,'matching_context_id':sid,'matched_seed_id':sid,
            'lagrangian_source_id':'synthetic:source_only','mu_f_gev':1e5,
            'mu_s_gev':1e3,'scheme':'MSbar','units':'GeV','synthetic_fixture':True}
    w={'format':W_FORMAT,'tensor_npz':'wilson_at_MS.npz',**common,
       'canonical_wilson_convention':W_CONVENTION,
       'wilson_evolved_to_ms':True,'includes_eft1_weinberg':True,
       'legacy_direct_llss_to_c5_log_added':False,
       'scalar_basis_id':f'{model}:canonical_interleaved',
       'lepton_indices':[[0,1],[2,3],[4,5]],
       'higgs_source_indices':list(range(4)),
       'heavy_source_indices':list(range(4,4+n))}
    b={'format':B_FORMAT,'tensor_npz':'background.npz',**common,
       'quartic_potential_convention':Q_CONVENTION,
       'quartic_inventory_status':'full_tensor',
       'wilson_scalar_basis_id':f'{model}:canonical_interleaved',
       'background_scalar_basis_id':f'{model}:canonical_interleaved',
       'single_common_scalar_threshold':True,
       'scalar_mass_origin':'synthetic_algebra_only',
       'sm_couplings_at_ms':{'gY':0.36,'g2':0.64,'g3':1.1,'lambdaH':0.25},
       'quartic_channel_names':[]}
    if declare_labels:
        w['scalar_source_real_labels']=scalar_labels(model)
        b['scalar_background_real_labels']=scalar_labels(model)
    # Nonzero, finite source: test numerical data are preserved, never rederived.
    C=np.zeros((8,8,4+n,4+n),complex)
    C[0,0,2,2]=3e-15+1e-15j
    C[0,0,4,4]=1e-12+2e-13j
    mass=np.diag(np.linspace(1.4e6,2.6e6,n))
    q=np.zeros((4+n,)*4,float)
    q[0,0,0,0]=0.6
    ye=np.diag([2e-6,6e-4,1e-2]).astype(complex)
    yu=np.diag([1e-5,4e-3,0.9]).astype(complex)
    yd=np.diag([3e-5,5e-4,2e-2]).astype(complex)
    np.savez_compressed(folder/'wilson_at_MS.npz',C_psi2phi2_ms=C)
    np.savez_compressed(folder/'background.npz',
                        mass_squared_gev2_ms=mass,
                        R_wilson_from_background=np.eye(n),
                        quartic_real_full_ms=q,ye=ye,yu=yu,yd=yd)
    for meta,name,npz in ((w,'wilson.json','wilson_at_MS.npz'),
                          (b,'background.json','background.npz')):
        meta['tensor_npz_sha256']=hashlib.sha256((folder/npz).read_bytes()).hexdigest()
        (folder/name).write_text(json.dumps(meta),encoding='utf-8')
    return folder/'wilson.json',folder/'background.json'


@pytest.mark.parametrize('model',ALL_MODELS)
def test_canonical_all16_without_labels_uses_declared_basis_only(tmp_path,model):
    w,b=_producer_pair(tmp_path/'producer',model)
    prepared=prepare_pair(model,w,b,tmp_path/'adapted'/model,allow_synthetic=True)
    wa=json.loads(prepared.wilson.read_text())
    ba=json.loads(prepared.background.read_text())
    for meta in (wa,ba):
        assert meta['model']==model
        assert meta['scheme']=='MSbar'
        assert meta['units']=='GeV'
        assert meta['mu_f_gev']==1e5 and meta['mu_s_gev']==1e3
        assert meta['synthetic_fixture'] is True
    assert wa['scalar_source_real_labels']==scalar_labels(model)
    assert ba['scalar_background_real_labels']==scalar_labels(model)
    with np.load(w.parent/'wilson_at_MS.npz',allow_pickle=False) as original, \
         np.load(prepared.wilson.parent/'wilson_at_MS.npz',allow_pickle=False) as adapted:
        np.testing.assert_array_equal(original['C_psi2phi2_ms'],adapted['C_psi2phi2_ms'])
    again=prepare_pair(model,prepared.wilson,prepared.background,
                       tmp_path/'repeat'/model,allow_synthetic=True)
    assert json.loads(again.background.read_text())['model']==model


def test_unknown_source_basis_does_not_get_invented_labels(tmp_path):
    w,b=_producer_pair(tmp_path/'producer','B_-1')
    wm=json.loads(w.read_text());bm=json.loads(b.read_text())
    wm['scalar_basis_id']='unverified_X_Y_frame'
    bm['wilson_scalar_basis_id']='unverified_X_Y_frame'
    w.write_text(json.dumps(wm));b.write_text(json.dumps(bm))
    with pytest.raises(AdapterBlocked,match='scalar_source_real_labels'):
        prepare_pair('B_-1',w,b,tmp_path/'adapted',allow_synthetic=True)


def test_genuine_run_rejects_synthetic_even_when_model_complete(tmp_path):
    w,b=_producer_pair(tmp_path/'producer','A_-4',declare_labels=True)
    with pytest.raises(AdapterBlocked,match='Synthetic input forbidden'):
        prepare_pair('A_-4',w,b,tmp_path/'adapted')
    assert not (tmp_path/'adapted'/'wilson.json').exists()


def test_native_B_permutation_and_explicit_full_scalar_labels(tmp_path):
    w,b=_producer_pair(tmp_path/'producer','A_-4',declare_labels=False)
    meta=json.loads(b.read_text())
    with np.load(b.parent/'background.npz',allow_pickle=False) as data:
        arrays={name:data[name].copy() for name in data.files}
    original=arrays['quartic_real_full_ms'].copy()
    labels=['H[1,R]','H[1,I]','H[2,R]','H[2,I]']+scalar_labels('A_-4')
    order=[4,0,5,1,2,3]+list(range(6,len(labels)))
    meta.pop('format');meta.pop('tensor_npz_sha256')
    meta['schema']='t3-parallel-B-background-v1'
    meta['scalar_full_real_labels']=[labels[i] for i in order]
    for axis in range(4):
        arrays['quartic_real_full_ms']=np.take(arrays['quartic_real_full_ms'],order,axis=axis)
    np.savez_compressed(b.parent/'background.npz',**arrays)
    b.write_text(json.dumps(meta),encoding='utf-8')
    result=prepare_pair('A_-4',w,b,tmp_path/'adapted',allow_synthetic=True)
    with np.load(result.background.parent/'background.npz',allow_pickle=False) as data:
        np.testing.assert_array_equal(original,data['quartic_real_full_ms'])
    bmeta=json.loads(result.background.read_text())
    assert bmeta['scalar_background_real_labels']==scalar_labels('A_-4')
    assert bmeta['scalar_labels_origin']=='B_full_real_scalar_labels'
