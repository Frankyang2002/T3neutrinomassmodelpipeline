"""Non-Weinberg source diagnostics: no suppressed or invented LLHH components."""
from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from Numerical.parallel_t3.threshold.matching import master_weinberg_template
from Numerical.parallel_t3.threshold.rough_projection import diagnose_llhh, diagnose_pair
from Numerical.parallel_t3.threshold.fixtures import make_fixture
from Numerical.parallel_t3.threshold.rough_adapter import prepare_pair
from Numerical.parallel_t3.threshold.models import ALL_MODELS

ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT/'scripts'/'RunT3RoughAll16.py'


def _amp(c5):
    return np.einsum('pq,ijAB->pqijAB',c5,master_weinberg_template(),optimize=True)


def test_exact_Weinberg_tensor_has_zero_projection_and_ward_residual():
    c5=np.array([[1+.25j,2,0],[2,3-.12j,0],[0,0,4]],complex)*1e-15
    d=diagnose_llhh(_amp(c5),master_weinberg_template(),declared_c5=c5)
    assert d['pure_Weinberg_at_tolerance_1e-8']
    assert d['relative_Weinberg_purity_residual']<1e-13
    assert d['C5_component_vs_LS_relative_to_larger']<1e-13
    assert all(x['relative_to_LLHH']<1e-13 for x in d['ward_SU2xU1'].values())
    assert d['stream2_separately_declared_C5']['relative_to_C_projector']<1e-13
    json.dumps(d,allow_nan=False)


def test_nonWeinberg_sources_are_recorded_with_ward_violations():
    c5=np.eye(3,dtype=complex)*1e-15
    a=_amp(c5)
    # Scalar and weak-leg symmetric, but NOT gauge singlet.
    a[0,0,0,0,0,0]+=9e-15
    d=diagnose_llhh(a,master_weinberg_template(),declared_c5=c5)
    assert d['relative_Weinberg_purity_residual']>.85
    assert d['ward_SU2xU1']['SU2_T1']['relative_to_LLHH']>.1
    assert d['ward_SU2xU1']['U1_Y']['relative_to_LLHH']>.1
    assert d['largest_nonWeinberg_components'][0]['absolute_difference_gev_inverse']>0
    assert not d['pure_Weinberg_at_tolerance_1e-8']


def test_higgs_conjugation_is_reported_as_diagnostic_only():
    a=_amp(np.eye(3,dtype=complex)*1e-14).conjugate()
    d=diagnose_llhh(a,master_weinberg_template())
    assert d['relative_Weinberg_purity_residual']>.1
    candidates=d['alternative_coordinate_hypotheses_DO_NOT_APPLY']
    assert candidates[0]['hypothesis_only']=='H_imaginary_components_conjugated'
    assert candidates[0]['relative_residual']<1e-13
    assert d['ward_SU2xU1']['U1_Y']['relative_to_LLHH']>.1


def test_source_neutral_component_vs_full_projection_is_exposed():
    c5=np.eye(3,dtype=complex)*1e-14
    a=_amp(c5)
    a[0,0,0,0,2,2]+=c5[0,0]*.5
    d=diagnose_llhh(a,master_weinberg_template(),declared_c5=c5)
    assert d['C5_component_vs_LS_relative_to_larger']>.01
    assert d['stream2_separately_declared_C5']['relative_to_neutral_component']>.01


def test_failed_run_still_persists_projection_diagnostic(tmp_path):
    a,b=make_fixture('A_-4',tmp_path/'source')
    with np.load(a.parent/'wilson_at_MS.npz',allow_pickle=False) as z:
        arrays={k:np.asarray(z[k]).copy() for k in z.files if k!='C5_EFT1_ms'}
    # This is intentionally invalid synthetic A: violates SU2 but leaves the
    # already-defined full Wilson pair symmetries exactly intact.
    arrays['C_psi2phi2_ms'][0,0,0,0]+=1e-13
    path=a.parent/'wilson_at_MS.npz'
    np.savez_compressed(path,**arrays)
    m=json.loads(a.read_text())
    m['tensor_npz_sha256']=sha256(path.read_bytes()).hexdigest()
    a.write_text(json.dumps(m),encoding='utf-8')
    pair=prepare_pair('A_-4',a,b,tmp_path/'adapted',allow_synthetic=True)
    diagnostic=diagnose_pair(pair.wilson,output=tmp_path/'report.json',model='A_-4')
    assert diagnostic['relative_Weinberg_purity_residual']>.5
    assert json.loads((tmp_path/'report.json').read_text())['status']=='SOURCE_DIAGNOSTIC_NOT_MATCHING'


def test_runner_does_not_convert_bad_weinberg_to_genuine_mass(tmp_path):
    a,b=make_fixture('A_-4',tmp_path/'source'/'A_m4')
    # Place the same canonical synthetic pair where the all16 runner looks.
    background=tmp_path/'B'/'A_m4';background.mkdir(parents=True)
    wilson=tmp_path/'A'/'A_m4';wilson.mkdir(parents=True)
    import shutil
    shutil.copyfile(a, wilson/'wilson.json')
    shutil.copyfile(a.parent/'wilson_at_MS.npz',wilson/'wilson_at_MS.npz')
    shutil.copyfile(b, background/'background.json')
    shutil.copyfile(b.parent/'background.npz',background/'background.npz')
    a=wilson/'wilson.json'
    with np.load(wilson/'wilson_at_MS.npz',allow_pickle=False) as z:
        arrays={k:np.asarray(z[k]).copy() for k in z.files if k!='C5_EFT1_ms'}
    arrays['C_psi2phi2_ms'][0,0,0,0]+=1e-13
    np.savez_compressed(wilson/'wilson_at_MS.npz',**arrays)
    wm=json.loads(a.read_text())
    wm['tensor_npz_sha256']=sha256((wilson/'wilson_at_MS.npz').read_bytes()).hexdigest()
    a.write_text(json.dumps(wm),encoding='utf-8')
    output=tmp_path/'results'
    command=[sys.executable,str(RUNNER),'--models','A_-4',
             '--background-root',str(tmp_path/'B'),
             '--wilson-root',str(tmp_path/'A'),
             '--output',str(output),'--allow-synthetic-smoke-test','--overwrite']
    result=subprocess.run(command,cwd=ROOT,capture_output=True,text=True)
    assert result.returncode==1,(result.stdout,result.stderr)
    info=json.loads((output/'summary.json').read_text())
    assert info['counts']['blocked']==1
    assert info['models'][0]['stage']=='existing_C_preflight'
    assert info['models'][0]['weinberg_purity_relative']>.5
    assert (output/'_diagnostics'/'A_m4_llhh_projection.json').is_file()
    assert not (output/'A_m4'/'C5_smeft_low.npz').exists()


def test_explicit_nonweinberg_projection_never_counts_as_real(tmp_path):
    a,b=make_fixture('A_-4',tmp_path/'source'/'A_m4')
    import shutil
    wilson=tmp_path/'A'/'A_m4';background=tmp_path/'B'/'A_m4'
    wilson.mkdir(parents=True);background.mkdir(parents=True)
    for src, dst in ((a,wilson/'wilson.json'),
                     (a.parent/'wilson_at_MS.npz',wilson/'wilson_at_MS.npz'),
                     (b,background/'background.json'),
                     (b.parent/'background.npz',background/'background.npz')):
        shutil.copyfile(src,dst)
    with np.load(wilson/'wilson_at_MS.npz',allow_pickle=False) as z:
        arrays={k:np.asarray(z[k]).copy() for k in z.files if k!='C5_EFT1_ms'}
    arrays['C_psi2phi2_ms'][0,0,0,0]+=1e-13
    np.savez_compressed(wilson/'wilson_at_MS.npz',**arrays)
    w=wilson/'wilson.json';m=json.loads(w.read_text())
    m['tensor_npz_sha256']=sha256((wilson/'wilson_at_MS.npz').read_bytes()).hexdigest()
    w.write_text(json.dumps(m),encoding='utf-8')
    output=tmp_path/'out'
    command=[sys.executable,str(RUNNER),'--models','A_-4',
             '--background-root',str(tmp_path/'B'),
             '--wilson-root',str(tmp_path/'A'),
             '--output',str(output),'--allow-synthetic-smoke-test',
             '--allow-projection-residual','--overwrite']
    proc=subprocess.run(command,cwd=ROOT,capture_output=True,text=True)
    assert proc.returncode==0,(proc.stdout,proc.stderr)
    summary=json.loads((output/'summary.json').read_text())
    assert summary['counts']['successful_real']==0
    assert summary['counts']['successful_synthetic_smoke']==0
    assert summary['counts']['nonweinberg_projected_diagnostics']==1
    manifest=json.loads((output/'A_m4'/'manifest.json').read_text())
    assert manifest['status']=='NonWeinbergProjectedDiagnostic'
    assert manifest['benchmark_approximation']['nonweinberg_prediction_not_valid_SMEFT']
    assert manifest['matching_residuals']['already_evolved_C5']>.5


def test_diagnostic_reports_orderzero_MF_witness_without_new_running(tmp_path):
    a,b=make_fixture('A_-4',tmp_path/'source')
    pair=prepare_pair('A_-4',a,b,tmp_path/'adapted',allow_synthetic=True)
    report=diagnose_pair(pair.wilson,output=tmp_path/'diag.json',model='A_-4')
    witness=report['source_MF_orderzero_witness']
    assert witness['availability']=='C_psi2phi2_tree_mf'
    assert witness['source_LLHH_norm_gev_inverse']>0
    assert witness['relative_Weinberg_purity_residual']<1e-10
    assert report['diagnostic_summary']['max_gauge_ward_relative']<1e-10


def test_projection_only_adapter_preserves_full_input_and_reports_redundant_c5(tmp_path):
    from Numerical.parallel_t3.threshold.rough_diagnostic_input import make_projection_only_pair
    from Numerical.parallel_t3.threshold.pipeline import _load_pair,W_FORMAT
    a,b=make_fixture('A_-4',tmp_path/'source')
    with np.load(a.parent/'wilson_at_MS.npz',allow_pickle=False) as z:
        arrays={k:np.asarray(z[k]).copy() for k in z.files}
    # Contaminate only a symmetric flavour, scalar and weak component.
    arrays['C_psi2phi2_ms'][0,0,2,2]+=1e-13
    np.savez_compressed(a.parent/'wilson_at_MS.npz',**arrays)
    wm=json.loads(a.read_text())
    wm['tensor_npz_sha256']=sha256((a.parent/'wilson_at_MS.npz').read_bytes()).hexdigest()
    a.write_text(json.dumps(wm),encoding='utf-8')
    pair=prepare_pair('A_-4',a,b,tmp_path/'adapted',allow_synthetic=True)
    original_sha=sha256(pair.wilson.parent.joinpath('wilson_at_MS.npz').read_bytes()).hexdigest()
    report=diagnose_pair(pair.wilson,output=tmp_path/'diag.json',model='A_-4')
    projected, audit=make_projection_only_pair('A_-4',pair,tmp_path/'projection',
        observed_purity_residual=report['relative_Weinberg_purity_residual'])
    assert audit['performed'] and audit['relative_disagreement']>0.01
    assert sha256(pair.wilson.parent.joinpath('wilson_at_MS.npz').read_bytes()).hexdigest()==original_sha
    orig,oa,_,_=_load_pair(pair.wilson,W_FORMAT)
    copied,ca,_,_=_load_pair(projected.wilson,W_FORMAT)
    assert 'C5_EFT1_ms' in oa and 'C5_EFT1_ms' not in ca
    np.testing.assert_array_equal(oa['C_psi2phi2_ms'],ca['C_psi2phi2_ms'])
    assert projected.background==pair.background
    assert copied['nonweinberg_projection_diagnostic_only']
    assert (tmp_path/'projection'/'redundant_C5_diagnostic.json').is_file()
    preserved, no_change=make_projection_only_pair('A_-4',pair,tmp_path/'none',
        observed_purity_residual=0)
    assert preserved is pair and not no_change['performed']


def test_optin_nonweinberg_redundant_c5_is_not_reported_as_genuine(tmp_path):
    from Numerical.parallel_t3.threshold.rough_diagnostic_input import make_projection_only_pair
    from Numerical.parallel_t3.threshold.rough_numerics import run_rough_case
    a,b=make_fixture('A_-4',tmp_path/'source')
    with np.load(a.parent/'wilson_at_MS.npz',allow_pickle=False) as z:
        arrays={k:np.asarray(z[k]).copy() for k in z.files}
    arrays['C_psi2phi2_ms'][0,0,2,2]+=1e-13
    np.savez_compressed(a.parent/'wilson_at_MS.npz',**arrays)
    wm=json.loads(a.read_text())
    wm['tensor_npz_sha256']=sha256((a.parent/'wilson_at_MS.npz').read_bytes()).hexdigest()
    a.write_text(json.dumps(wm),encoding='utf-8')
    pair=prepare_pair('A_-4',a,b,tmp_path/'adapted',allow_synthetic=True)
    report=diagnose_pair(pair.wilson,output=tmp_path/'diagnosis.json',model='A_-4')
    probe, audit=make_projection_only_pair('A_-4',pair,tmp_path/'projection',
        observed_purity_residual=report['relative_Weinberg_purity_residual'])
    assert audit['performed']
    # Test uses the actual existing C solver and SM+Weinberg evolution.
    result=run_rough_case('A_-4',probe,tmp_path/'results',allow_synthetic=True,
                          allow_projection_residual=True)
    assert result['status']=='NonWeinbergProjectedDiagnostic'
    meta=json.loads((tmp_path/'results'/'manifest.json').read_text())
    assert meta['input_producer_sha256']['projected_diagnostic_input_only']['performed']
    assert meta['benchmark_approximation']['nonweinberg_prediction_not_valid_SMEFT']


def test_flattened_ward_report_exposes_numbers_not_nested_objects():
    from scripts.InspectT3RoughWardSummary import flatten_ward
    c5=np.eye(3,dtype=complex)*1e-15
    d=diagnose_llhh(_amp(c5),master_weinberg_template())
    d['model']='A_-4'
    flat=flatten_ward(d)
    assert flat['model']=='A_-4'
    assert flat['Weinberg_residual']<1e-12
    assert flat['max_Ward']<1e-12
    assert flat['MF_witness_kind']=='not_present_in_older_report'
    assert all(isinstance(flat[key],float) for key in ('SU2_T1','SU2_T2','SU2_T3','U1_Y'))


def test_overwrite_removes_stale_projected_result_before_strict_block(tmp_path):
    out=tmp_path/'out'
    old=out/'A_m4'
    old.mkdir(parents=True)
    (old/'manifest.json').write_text(json.dumps({
        'model':'A_-4','status':'NonWeinbergProjectedDiagnostic'}))
    (old/'C5_smeft_low.npz').write_bytes(b'old diagnostic file, not a valid current result')
    command=[sys.executable,str(RUNNER),'--models','A_-4',
             '--background-root',str(tmp_path/'missing_B'),
             '--wilson-root',str(tmp_path/'missing_A'),
             '--output',str(out),'--overwrite']
    proc=subprocess.run(command,cwd=ROOT,text=True,capture_output=True)
    assert proc.returncode==1,(proc.stdout,proc.stderr)
    assert not old.exists()
    status=json.loads((out/'summary.json').read_text())
    assert status['counts']['blocked']==1
    assert status['models'][0]['stale_results_removed_for_overwrite'] is True
