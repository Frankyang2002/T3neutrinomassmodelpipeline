"""Batch 4: analytical identities, negative controls and 16-model read-only fixture."""
from __future__ import annotations

from hashlib import sha256
import json
import math
import os
from pathlib import Path

import numpy as np
import pytest

from Numerical.lo_certification.direct_uv.native_amplitude_contract import model_dictionary
from Numerical.lo_patterns_16 import flavour_geometry as x


def pack(z):
    z=np.asarray(z,dtype=complex)
    return {'real':z.real.tolist(),'imag':z.imag.tolist()}


def records(*hs):
    return [{'matrix_ev':np.asarray(h,dtype=complex),
             'mass_gev':float(10000*(i+1)),'multiplicity':1}
            for i,h in enumerate(hs)]


def put(path:Path,obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,allow_nan=False,indent=2)+'\n',encoding='utf-8')


def sha(path):return sha256(path.read_bytes()).hexdigest()


def test_exact_mode_factorization_complex_interference():
    one=np.diag([1.+.3j,2.,3.])
    two=np.diag([-.98-.3j,.5,1.])
    got=x.geometry(one+two,records(one,two),'synthetic')
    tiny=got['modes'][0]
    assert tiny['mass_ev']==pytest.approx(.02,rel=1e-12)
    assert tiny['coherence_fraction']<.02
    assert tiny['mechanism_diagnostic']=='heavy_eigenspace_cancellation_dominant'
    assert tiny['light_mass_fraction_of_group_strength']==pytest.approx(
        tiny['geometric_projection_fraction']*tiny['coherence_fraction'],abs=1e-12)
    assert max(m['product_identity_absolute_error'] for m in got['modes'])<1e-12


def test_small_lightest_mass_from_flavour_projection_not_cancellation():
    one=np.diag([1e-7,1.,2.]);two=np.diag([2e-7,2.,1.])
    got=x.geometry(one+two,records(one,two),'synthetic')
    lo=got['modes'][0]
    assert lo['mass_ev']==pytest.approx(3e-7,rel=1e-10)
    assert lo['coherence_fraction']==pytest.approx(1.,rel=1e-10)
    assert lo['geometric_projection_fraction']<1e-6
    assert lo['mechanism_diagnostic']=='mode_geometry_dominant'


def test_nonzero_singular_masses_of_complex_symmetric_matrix():
    H=np.array([[2+.7j,1.3-.1j,0],[1.3-.1j,1.-.2j,.4],[0,.4,3.]],complex)
    got=x.geometry(H,records(H),'synthetic')
    assert np.allclose(got['masses_ev'],np.sort(np.linalg.svd(H,compute_uv=False)))
    assert all(m['coherence_fraction']==pytest.approx(1.) for m in got['modes'])
    assert got['determinant_absolute_ev3']==pytest.approx(abs(np.linalg.det(H)),rel=1e-12)
    assert np.allclose([p['amplitude_ev']['imag'] for p in got['projections']],0,atol=1e-12)


def test_unitary_flavour_rephasing_invariant():
    H=np.array([[1+1j,.4-.1j,.5],[.4-.1j,2+0j,.3j],[.5,.3j,3]],complex)
    K=np.array([[.8,0,.2],[0,-.4j,.5],[.2,.5,.5]],complex)
    basis=np.diag(np.exp(1j*np.array([.8,-1.2,.3])))
    a=x.geometry(H+K,records(H,K),'native')
    b=x.geometry(basis.T@(H+K)@basis,
                 records(basis.T@H@basis,basis.T@K@basis),'rephased')
    assert np.allclose(a['masses_ev'],b['masses_ev'],rtol=1e-12)
    for i in range(3):
        for k in ('geometric_projection_fraction','coherence_fraction'):
            assert a['modes'][i][k]==pytest.approx(b['modes'][i][k],abs=1e-12)


def test_heavy_group_sum_is_not_sum_of_individual_mass_spectra():
    H=np.diag([1.,0,0]);K=np.diag([-1.,2.,3.])
    g=x.geometry(H+K,records(H,K),'synthetic')
    assert g['masses_ev']==pytest.approx([0.,2.,3.],abs=1e-12)
    assert not np.allclose(g['masses_ev'],x.sv(H)+x.sv(K))


def test_group_decomposition_rejects_missing_or_rescaled_group():
    H=np.diag([1.,2.,3.]);K=np.diag([.1,.2,.3])
    with pytest.raises(x.FlavourError,match='sum does not reproduce'):
        x.geometry(H+K,records(H),'synthetic')
    with pytest.raises(x.FlavourError,match='sum does not reproduce'):
        x.geometry(H+K,records(H,2*K),'synthetic')


def test_group_rejects_nonsymmetric_and_nonfinite():
    H=np.eye(3,dtype=complex)
    bad=H.copy();bad[0,1]=.3
    with pytest.raises(x.FlavourError,match='not symmetric'):
        x.geometry(H,records(bad,H-bad),'synthetic')
    with pytest.raises(x.FlavourError,match='invalid mass matrix'):
        x.geometry(H*np.nan,records(H),'synthetic')
    with pytest.raises(x.FlavourError,match='nonfinite'):
        x.matrix(pack(np.full((3,3),np.nan)),'synthetic')


def test_determinant_suppressed_without_small_frobenius_norm():
    a=np.diag([1e-5,1,2.]);g=x.geometry(a,records(a),'synthetic')
    assert g['matrix_frobenius_ev']>2.
    assert g['determinant_scaled_by_frobenius_cubed']<1e-5
    assert g['masses_ev'][0]/g['masses_ev'][-1]==pytest.approx(5e-6)


def test_loo_mass_changes_are_not_matrix_part_eigenmasses():
    A=np.diag([.1,1,3]);B=np.diag([1,1,1])
    g=x.geometry(A+B,records(A,B),'synthetic')
    assert g['leave_one_group_out'][0]['m1_without_ev']==pytest.approx(1.)
    assert g['leave_one_group_out'][1]['m1_without_ev']==pytest.approx(.1)
    assert g['leave_one_group_out'][0]['m1_change_ev']==pytest.approx(-.1)


def test_single_group_coherence_is_unity():
    A=np.diag([.2,1,2]);g=x.geometry(A,records(A),'synthetic')
    assert all(m['coherence_fraction']==pytest.approx(1) for m in g['modes'])
    assert g['fermion_eigenspace_coherence']==pytest.approx(1)


def mock_project(tmp_path,monkeypatch):
    """Self-consistent synthetic 16-model fixture, not a physics benchmark."""
    root=tmp_path/'repo';root.mkdir()
    old=root/x.STUDY;prev=root/x.BATCH02;third=root/x.BATCH03
    allrows=[];r2=[];r3=[];fingerprints={}
    for j,item in enumerate(model_dictionary()):
        label=item['model']
        e1=np.diag([.0005+.00003*j,.008+.0001*j,.04+.0005*j]);
        e2=np.diag([.0001+.00001*j,.004+.0001*j,.012+.0001*j]);
        e3=np.diag([.0002,.001,.004]);
        if j==8:
            e1[0,0]=-.04; e2[0,0]=.04 # C+1: synthetic heavy-state cancellation
        parts=[e1,e2,e3];full=sum(parts)
        m=np.sort(np.linalg.svd(full,compute_uv=False))
        filename=root/('Numerical/benchmarks/fitted_reference/'+label+'.json')
        inputfile=root/('output/general_observable/lo_direct_one_loop/provisional_inputs/'+label+'.json')
        put(filename,{'model':label,'synthetic':True})
        FM=np.diag([10000.,20000.,30000.]);put(inputfile,{'fermion_mass_gev':pack(FM)})
        source={'benchmark':str(filename.relative_to(root)).replace('\\','/'),
                'benchmark_sha256':sha(filename),
                'input':str(inputfile.relative_to(root)).replace('\\','/'),
                'input_sha256':sha(inputfile)}
        saved={'model':label,'status':'conditional_uv_evaluated',
               'physical_normalisation_certified':False,'source':source,
               'mnu_ev':pack(full),'takagi_masses_ev_sorted':m.tolist()}
        p=old/'models'/f'{label}.json';put(p,saved)
        B={'model':label,'benchmark_sha256':sha(filename),'component_input_sha256':sha(inputfile)}
        allrows.append(B);r2.append(B.copy())
        nh=[{'fermion_mass_gev':10000.*(i+1),'eigenspace_multiplicity':1,
             'matrix_ev':pack(H),'matrix_norm_ev':float(np.linalg.norm(H))}
            for i,H in enumerate(parts)]
        groupcoh=float(np.linalg.norm(full)/sum(np.linalg.norm(h) for h in parts))
        modelrow={'model':label,'fermion_type':item['fermion_type'],
                  'physical_normalisation_certified':False,
                  'yukawa_frobenius_product':1.5e-8,'lambdaT3_abs':.1,
                  'heavy_mass_groups':nh,'heavy_mass_group_coherence':groupcoh}
        r3.append(modelrow)
        per={'model':label,'source_SHA256':source,'summary':modelrow,
             'physical_normalisation_certified':False,
             'observed_source_conditional_mnu_gev':pack(full*1e-9),
             'exact_closed_kernel_mnu_gev':pack(full*1e-9)}
        p3=third/'models'/f'{label}.json';put(p3,per)
        fingerprints[label]={'model_detail_sha256':sha(p),
                             'benchmark_sha256':sha(filename),
                             'weyl_input_sha256':sha(inputfile)}
    summary={'format':'t3-conditional-pattern-study-all16-v1',
             'counts':{'conditional_evaluated':16,'blocked':0,'physical_source_certified':0},
             'physical_normalisation_certified':False,'models':allrows,
             'provenance':{'revision':{'local_head':'synthetic'},
                           'source_audit':{'source_files':{}},'source_audit_sha256':{}}}
    put(old/'all16_summary.json',summary)
    b2={'format':'t3-chatA-readonly-batch02-numerical-patterns-v1',
        'physical_normalisation_certified':False,
        'input_run_sha256':sha(old/'all16_summary.json'),
        'counts':{'models':16,'source_scan_rows':288,'physical_source_certified':0},'models':r2}
    put(prev/'batch02_analysis.json',b2)
    b3={'format':'t3-chatA-batch03-analytic-scalar-kernel-v1',
        'physical_normalisation_certified':False,'safe_for_physical_source_certification':False,
        'source_summary_sha256':sha(old/'all16_summary.json'),
        'model_count':16,'models':r3}
    put(third/'batch03_analytic_explanation.json',b3)
    man={'source_summary_sha256':sha(old/'all16_summary.json'),
         'batch02_analysis_sha256':sha(prev/'batch02_analysis.json'),
         'physical_source_certified':False,'models':fingerprints}
    put(third/'SOURCE_MANIFEST.json',man)
    # This test fixture explicitly lacks the production source-audit artefacts.
    # The production hash guard has independent tests in Batch 2 and remains
    # enabled in the deployed workflow.
    monkeypatch.setattr(x,'_source_hash_guard',lambda *a,**k:None)
    return root


def test_synthetic_all16_complete_and_original_files_unchanged(tmp_path,monkeypatch):
    root=mock_project(tmp_path,monkeypatch)
    before={str(p):sha(p) for p in root.rglob('*') if p.is_file()}
    result=x.analyse(root)
    assert result['model_count']==16
    assert result['physical_normalisation_certified'] is False
    assert len(result['models'])==16
    assert (root/x.TARGET/'BATCH04_FLAVOUR_REPORT.md').is_file()
    assert (root/x.TARGET/'batch04_mode_projections.csv').is_file()
    assert all(row['matrix_reconstruction_relative_error']<1e-12 for row in result['models'])
    assert any(row['lightest_mode']['mechanism_diagnostic']=='heavy_eigenspace_cancellation_dominant'
               for row in result['models'])
    assert all(sha(Path(p))==h for p,h in before.items())


def test_fails_closed_on_modified_batch3_eigenspace(tmp_path,monkeypatch):
    root=mock_project(tmp_path,monkeypatch)
    p=root/x.BATCH03/'batch03_analytic_explanation.json';data=json.loads(p.read_text())
    data['models'][0]['heavy_mass_groups'][0]['matrix_ev']['real'][0][0]*=2
    put(p,data)
    with pytest.raises(x.FlavourError,match='summary drift|norm drift'):
        x.analyse(root)
    assert not (root/x.TARGET).exists()


def test_fails_closed_on_changed_frozen_input(tmp_path,monkeypatch):
    root=mock_project(tmp_path,monkeypatch)
    p=root/'output/general_observable/lo_direct_one_loop/provisional_inputs/A_-4.json'
    raw=json.loads(p.read_text());raw['fermion_mass_gev']['real'][0][0]+=100
    put(p,raw)
    # On the real checkout _source_hash_guard checks this SHA before decoding.
    # Our minimal fixture bypasses that external guard; spectral mass mismatch
    # must nevertheless block the modified input.
    with pytest.raises(x.FlavourError,match='do not match source'):
        x.analyse(root)
    assert not (root/x.TARGET).exists()


def test_fails_closed_on_modified_batch3_manifest(tmp_path,monkeypatch):
    root=mock_project(tmp_path,monkeypatch)
    p=root/x.BATCH03/'SOURCE_MANIFEST.json';raw=json.loads(p.read_text())
    raw['models']['C_-1']['model_detail_sha256']='0'*64;put(p,raw)
    with pytest.raises(x.FlavourError,match='provenance conflict'):
        x.analyse(root)
    assert not (root/x.TARGET).exists()


def test_fails_closed_on_false_certification(tmp_path,monkeypatch):
    root=mock_project(tmp_path,monkeypatch)
    p=root/x.BATCH03/'batch03_analytic_explanation.json';raw=json.loads(p.read_text())
    raw['physical_normalisation_certified']=True;put(p,raw)
    with pytest.raises(x.FlavourError,match='physical certification'):
        x.analyse(root)
    assert not (root/x.TARGET).exists()


def test_fails_closed_on_existing_output(tmp_path,monkeypatch):
    root=mock_project(tmp_path,monkeypatch)
    x.analyse(root)
    with pytest.raises(x.FlavourError,match='already populated'):
        x.analyse(root)


def test_fails_closed_on_output_outside_owned_tree(tmp_path,monkeypatch):
    root=mock_project(tmp_path,monkeypatch)
    with pytest.raises(x.FlavourError,match='outside Chat A output tree'):
        x.analyse(root,output_dir=root/'some_elsewhere')


@pytest.mark.skipif(os.getenv('RUN_LO_PATTERN_INTEGRATION')!='1',
                    reason='Opt-in: requires local full frozen Batch-1/2/3 outputs')
def test_actual_local_stored_all16_analytical_integration(tmp_path):
    root=Path(__file__).resolve().parents[1]
    dest=root/x.BASE/'analysis_batch04_integration_test_01'
    # Integration is not run automatically because it writes a complete report.
    # Allow an empty directory only; reruns must choose a fresh directory.
    result=x.analyse(root,output_dir=dest)
    assert result['model_count']==16
    assert max(r['matrix_reconstruction_relative_error'] for r in result['models'])<1e-7
