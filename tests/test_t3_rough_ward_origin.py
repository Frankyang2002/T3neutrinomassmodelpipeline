"""Independent scalar-coordinate charge/Ward tests (synthetic algebra only)."""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from hashlib import sha256

import numpy as np
import pytest

from Numerical.parallel_t3.threshold.matching import master_weinberg_template
from Numerical.parallel_t3.threshold.models import n_real
from Numerical.parallel_t3.threshold.pipeline import W_FORMAT
from Numerical.parallel_t3.threshold.rough_origin import (
    source_charge_analysis,real_to_complex_llhh,complex_to_real_llhh,
    complex_higgs_coordinate_map,audit_pair)

ROOT=Path(__file__).resolve().parents[1]


def _pure():
    c5=np.array([[.3+.2j,.07,0],[.07,.6-.1j,.04],[0,.04,.15]],complex)*1e-15
    return np.einsum('pq,ijAB->pqijAB',c5,master_weinberg_template(),optimize=True)


def _charge_example(kind):
    z=real_to_complex_llhh(_pure())
    if kind=='noninvariant_HH':
        z[0,0,0,0,0,0]+=1e-14
    if kind=='wrong_charge':
        z[0,0,0,0,0,2]+=1e-14
        z[0,0,0,0,2,0]+=1e-14
        z[0,0,0,0,2,2]+=2e-14
    return complex_to_real_llhh(z)


def test_unitary_complex_higgs_map_and_exact_roundtrip():
    M=complex_higgs_coordinate_map()
    np.testing.assert_allclose(M.conj().T@M,np.eye(4),rtol=0,atol=1e-14)
    a=_pure()
    np.testing.assert_allclose(complex_to_real_llhh(real_to_complex_llhh(a)),a,rtol=0,atol=1e-29)


def test_pure_weinberg_is_only_holomorphic_and_gauge_singlet():
    a=source_charge_analysis(_pure())
    assert a['status']=='consistent_with_unique_Weinberg_singlet'
    assert a['charge_sector_fractions']['HH']==pytest.approx(1)
    assert a['charge_sector_fractions']['H_Hbar_plus_Hbar_H']<1e-14
    assert a['charge_sector_fractions']['Hbar_Hbar']<1e-14
    assert a['allowed_HH_singlet_residual_relative']<1e-12
    assert a['direct_U1_Ward_relative']<1e-12
    assert max(a['allowed_HH_SU2_Ward_relative'].values())<1e-12


def test_bad_su2_can_preserve_hypercharge_exactly():
    x=source_charge_analysis(_charge_example('noninvariant_HH'))
    assert x['status']=='correct_HH_hypercharge_but_broken_SU2_tensor'
    assert x['direct_U1_Ward_relative']<1e-12
    assert x['charge_sector_fractions']['HH']>1-1e-12
    assert x['allowed_HH_singlet_residual_relative']>.5
    assert max(x['allowed_HH_SU2_Ward_relative'].values())>.5


def test_wrong_hypercharge_is_independent_of_su2_breakdown():
    x=source_charge_analysis(_charge_example('wrong_charge'))
    assert x['status']=='wrong_hypercharge_sector_in_LLHH'
    assert x['direct_U1_Ward_relative']>.5
    assert x['charge_sector_fractions']['H_Hbar_plus_Hbar_H']>.1
    assert x['charge_sector_fractions']['Hbar_Hbar']>.1
    assert x['direct_U1_Ward_relative']==pytest.approx(x['independently_reconstructed_U1_Ward_relative'],abs=1e-13)


def test_zero_LLHH_reports_no_relative_source_information():
    r=source_charge_analysis(np.zeros((3,3,2,2,4,4),complex))
    assert r['status']=='zero_LLHH_no_relative_source_information'
    assert r['direct_U1_Ward_relative']==0


def _write_source(tmp_path,model='A_-4',status='noninvariant_HH'):
    n=n_real(model)
    full=np.zeros((6,6,4+n,4+n),complex)
    tree=np.zeros_like(full)
    good=_pure();bad=_charge_example(status)
    for p in range(3):
        for q in range(3):
            for i in range(2):
                for j in range(2):
                    full[2*p+i,2*q+j,:4,:4]=bad[p,q,i,j]
                    tree[2*p+i,2*q+j,:4,:4]=good[p,q,i,j]
    loc=tmp_path/model.replace('_-','_m').replace('_+','_p')
    loc.mkdir(parents=True,exist_ok=True)
    npz=loc/'wilson_at_MS.npz'
    np.savez_compressed(npz,C_psi2phi2_ms=full,C_psi2phi2_tree_mf=tree)
    meta={'format':W_FORMAT,'model':model,'tensor_npz':npz.name,
          'tensor_npz_sha256':sha256(npz.read_bytes()).hexdigest(),
          'lepton_indices':[[0,1],[2,3],[4,5]],
          'higgs_source_indices':[0,1,2,3],
          'heavy_source_indices':list(range(4,4+n)),
          'scalar_basis_id':f'{model}:canonical_interleaved',
          'synthetic_fixture':True}
    manifest=loc/'wilson.json'
    manifest.write_text(json.dumps(meta),encoding='utf8')
    return manifest


@pytest.mark.parametrize('model',['A_-4','B_-1','E_+2'])
def test_source_origin_comparison_uses_saved_tree_not_guessed_boundary(tmp_path,model):
    p=_write_source(tmp_path,model)
    d=audit_pair(p,model=model)
    assert d['status']=='SourceWardOriginDiagnosticOnly'
    assert d['MS_full_evolved_LLHH']['status']=='correct_HH_hypercharge_but_broken_SU2_tensor'
    assert d['MF_orderzero_LLHH']['status']=='consistent_with_unique_Weinberg_singlet'
    assert 'arose_after_saved_tree_endpoint' in d['origin_interpretation_NOT_proof']


def test_hypercharge_sector_gauge_break_and_provenance_hash(tmp_path):
    p=_write_source(tmp_path,'B_-1','wrong_charge')
    d=audit_pair(p,model='B_-1')
    assert d['MS_full_evolved_LLHH']['status']=='wrong_hypercharge_sector_in_LLHH'
    assert d['source_sha256']['npz']==sha256((p.parent/'wilson_at_MS.npz').read_bytes()).hexdigest()
    np.testing.assert_allclose(d['MS_full_evolved_LLHH']['direct_U1_Ward_relative'],
        d['MS_full_evolved_LLHH']['independently_reconstructed_U1_Ward_relative'])
    with pytest.raises(ValueError,match='model mismatch|model mismatch'):
        audit_pair(p,model='C_-1')


def test_cli_writes_csv_and_does_not_modify_source_files(tmp_path):
    base=tmp_path/'results';folder=base/'_adapted_inputs'
    p=_write_source(folder,'A_-4')
    before=(p.read_bytes(),(p.parent/'wilson_at_MS.npz').read_bytes())
    cmd=[sys.executable,str(ROOT/'scripts/AuditT3RoughWardOrigin.py'),
         '--models','A_-4','--root',str(base)]
    done=subprocess.run(cmd,cwd=ROOT,text=True,capture_output=True)
    assert done.returncode==0,(done.stdout,done.stderr)
    out=base/'_diagnostics'/'origin_v5'
    assert (out/'origin_summary.csv').is_file()
    document=json.loads((out/'A_m4_origin.json').read_text())
    assert document['source_sha256']['npz']==sha256(before[1]).hexdigest()
    assert p.read_bytes()==before[0]
    assert (p.parent/'wilson_at_MS.npz').read_bytes()==before[1]
