"""C v2: real-endpoint readiness and fixed-order/RGE-improved comparison.

All A/B tensors here are expressly synthetic *contract fixtures*, not physics.
The preflight/CLI must never automatically treat them as real endpoints.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from Numerical.parallel_t3.threshold.fixtures import make_fixture
from Numerical.parallel_t3.threshold.models import ALL_MODELS, n_real
from Numerical.parallel_t3.threshold.preflight import PREFLIGHT_FORMAT, preflight_pair

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT/'scripts/RunT3ParallelCRealIntegration.py'


def _update_npz(path: Path, **updates):
    with np.load(path,allow_pickle=False) as z:
        arrays={key:np.array(z[key]) for key in z.files}
    arrays.update(updates)
    np.savez_compressed(path,**arrays)
    meta=path.parent/'wilson.json' if path.name.startswith('wilson') else path.parent/'background.json'
    data=json.loads(meta.read_text(encoding='utf-8'))
    data['tensor_npz_sha256']=hashlib.sha256(path.read_bytes()).hexdigest()
    meta.write_text(json.dumps(data),encoding='utf-8')


def test_missing_A_B_reports_blocked_and_does_not_create_tensors(tmp_path):
    a=tmp_path/'A_m4'/'wilson.json'
    b=tmp_path/'A_m4'/'background.json'
    pf=preflight_pair(a,b,expected_model='A_-4')
    assert pf['format']==PREFLIGHT_FORMAT
    assert pf['status']=='blocked'
    assert len(pf['blockers'])==2
    assert not a.exists() and not b.exists()
    json.dumps(pf,allow_nan=False)


@pytest.mark.parametrize('model',ALL_MODELS)
def test_all16_contract_preflight_strict_without_writes(tmp_path,model):
    a,b=make_fixture(model,tmp_path)
    hashes=(hashlib.sha256(a.read_bytes()).hexdigest(),
            hashlib.sha256(b.read_bytes()).hexdigest())
    pf=preflight_pair(a,b,mode='strict_one_loop',expected_model=model,
                      permit_synthetic=True)
    assert pf['status']=='ready_conditional',pf['blockers']
    result=pf['checks']['single_MS_finite_bubble_match_candidate']['details']
    assert result['eigenstates_integrated_out_at_single_MS']==n_real(model)
    assert result['unordered_eigenstate_pairs']==n_real(model)*(n_real(model)+1)//2
    assert result['already_evolved_C5_kept_once']
    assert not pf['real_physical_prediction_certified']
    assert pf['input_kind']=='synthetic_contract_probe'
    json.dumps(pf,allow_nan=False)
    assert hashes==(hashlib.sha256(a.read_bytes()).hexdigest(),
                    hashlib.sha256(b.read_bytes()).hexdigest())


def test_synthetic_rejected_without_optin(tmp_path):
    a,b=make_fixture('A_-4',tmp_path)
    pf=preflight_pair(a,b,mode='rge_improved')
    assert pf['status']=='blocked'
    assert 'synthetic fixture forbidden' in str(pf['blockers'])


def test_full_Wilson_symmetry_error_found_before_run(tmp_path):
    a,b=make_fixture('B_-3',tmp_path)
    with np.load(a.parent/'wilson_at_MS.npz',allow_pickle=False) as z:
        full=np.array(z['C_psi2phi2_ms'])
    full[0,1,5,5]=5e-13+3e-13j
    _update_npz(a.parent/'wilson_at_MS.npz',C_psi2phi2_ms=full)
    pf=preflight_pair(a,b,permit_synthetic=True)
    assert pf['status']=='blocked'
    assert 'full fermion/scalar pair symmetry failed' in str(pf['blockers'])


def test_stale_source_pair_rejected_before_any_matching(tmp_path):
    a,b=make_fixture('A_-4',tmp_path)
    meta=json.loads(b.read_text(encoding='utf-8'))
    meta['lagrangian_source_id']='STALE_PRE_PORTAL'
    b.write_text(json.dumps(meta),encoding='utf-8')
    pf=preflight_pair(a,b,permit_synthetic=True)
    assert pf['status']=='blocked'
    assert 'lagrangian_source_id' in str(pf['blockers'])


def test_missing_quartic_channel_detected(tmp_path):
    a,b=make_fixture('C_-1',tmp_path)
    p=b.parent/'background.npz'
    with np.load(p,allow_pickle=False) as z:
        bad=np.array(z['quartic_real_channels_ms'])
    bad[0]*=.5
    _update_npz(p,quartic_real_channels_ms=bad)
    pf=preflight_pair(a,b,permit_synthetic=True)
    assert pf['status']=='blocked'
    assert 'channels omit/double-count' in str(pf['blockers'])


def test_positive_mass_squared_preflight_guard(tmp_path):
    a,b=make_fixture('D_+0',tmp_path)
    p=b.parent/'background.npz'
    with np.load(p,allow_pickle=False) as z:
        m=np.array(z['mass_squared_gev2_ms'])
    m[0,0]=-1.0e8
    _update_npz(p,mass_squared_gev2_ms=m)
    pf=preflight_pair(a,b,permit_synthetic=True)
    assert pf['status']=='blocked'
    assert 'Tachyonic' in str(pf['blockers'])


def test_double_count_direct_log_rejected(tmp_path):
    a,b=make_fixture('A_-4',tmp_path)
    meta=json.loads(a.read_text(encoding='utf-8'))
    meta['legacy_direct_llss_to_c5_log_added']=True
    a.write_text(json.dumps(meta),encoding='utf-8')
    pf=preflight_pair(a,b,permit_synthetic=True)
    assert pf['status']=='blocked'
    assert 'log' in str(pf['blockers'])


def test_script_run_both_modes_and_compare(tmp_path):
    a,b=make_fixture('A_-4',tmp_path/'source')
    output=tmp_path/'output'
    command=[sys.executable,str(SCRIPT),'--wilson',str(a),
             '--background',str(b),'--mode','both','--models','A_-4',
             '--permit-synthetic','--match-only','--run','--output',str(output)]
    proc=subprocess.run(command,cwd=ROOT,text=True,capture_output=True)
    assert proc.returncode==0,(proc.stdout,proc.stderr)
    report=json.loads((output/'real_integration_summary.json').read_text())
    assert report['counts']['all_requested_modes_ready_or_completed']==1
    diff=json.loads((output/'A_m4'/'counting_comparison.json').read_text())
    assert diff['no_direct_LLSS_to_C5_log_added_twice']
    assert (diff['scalar_bubble_rge_improved_minus_strict']
            ['difference_frobenius_gev_inverse']>1e-21)
    assert (output/'A_m4'/'strict_one_loop'/'manifest.json').exists()
    assert (output/'A_m4'/'rge_improved'/'manifest.json').exists()


def test_script_missing_real_inputs_blocks_without_fixture_fallback(tmp_path):
    output=tmp_path/'out'
    proc=subprocess.run([sys.executable,str(SCRIPT),
            '--wilson-root',str(tmp_path/'not_A'),
            '--background-root',str(tmp_path/'not_B'),
            '--models','A_-4','--output',str(output)],
            cwd=ROOT,text=True,capture_output=True)
    assert proc.returncode==1
    report=json.loads((output/'real_integration_summary.json').read_text())
    assert report['counts']['blocked']==1
    assert list(output.rglob('C5_*.npz'))==[]


def test_strict_requires_independently_truncated_A_endpoint(tmp_path):
    a,b=make_fixture('A_-4',tmp_path)
    p=a.parent/'wilson_at_MS.npz'
    with np.load(p,allow_pickle=False) as z:
        arrays={key:np.array(z[key]) for key in z.files
                if key!='C_psi2phi2_ms_strict_1loop'}
    np.savez_compressed(p,**arrays)
    manifest=json.loads(a.read_text(encoding='utf-8'))
    manifest['tensor_npz_sha256']=hashlib.sha256(p.read_bytes()).hexdigest()
    a.write_text(json.dumps(manifest),encoding='utf-8')
    strict=preflight_pair(a,b,mode='strict_one_loop',permit_synthetic=True)
    improved=preflight_pair(a,b,mode='rge_improved',permit_synthetic=True)
    assert strict['status']=='blocked'
    assert 'one-loop-truncated' in str(strict['blockers'])
    assert improved['status']=='ready_conditional'
