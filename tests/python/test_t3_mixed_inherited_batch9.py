"""Offline source-coupling and fail-closed tests for Batch 9 (no Wolfram)."""
from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import sympy as sp

from scripts import RunT3MixedInheritedPoleBatch9 as task


def inp():
    y1=np.diag([0.1+0.2j,0.2+0.1j,0.3+0.2j]).astype(complex)
    y2=np.diag([0.4-0.1j,0.5-0.1j,0.6-0.1j]).astype(complex)
    return {'mf':np.array([100_000.,110_000.,120_000.]),
            'y1':y1,'y2':y2,'U':np.eye(2,dtype=complex),
            'couplings':{'lambdaT3':1.0+0.1j,
                'lambdaT3ExtraQ020020I1':0j,
                'lambdaT3ExtraQ020200I1':0j}}


def test_matrix_orientation_transpose_not_dagger():
    p=np.array([1+2j,3-1j,4j])
    x=task.coefficient_matrix(p)
    assert x.shape==(2,2) and np.array_equal(x,x.T)
    assert not np.array_equal(x,x.conj().T)
    assert x[0,1]==p[1] and x[0,0]==2*p[0]


def test_f_source_sum_three_masses_and_log_scale():
    v={'pole':np.array([1+0j,2+1j,3+0j]),
       'log':np.array([1+0j,2+1j,3+0j]),
       'constant':np.array([1+0j,2+1j,3+0j])}
    f,f0,record=task.f_inherited_matrices(inp(),v,8_000.)
    _,f1,_=task.f_inherited_matrices(inp(),v,8_000.*np.exp(0.023))
    assert len(record)==3
    assert task.relerr((f1-f0)/(2*0.023),f)<1e-10
    assert not np.allclose(f,f.conj().T,atol=0,rtol=0)
    assert np.linalg.norm(f-f.T)<1e-20


def test_normalization_identifies_exact_half_factor_failure():
    # Given P_F = -2 P_s, Batch9 v10 reported 0.5 rather than a
    # vanishing total pole. Only the source-established A5 -> C5 conversion
    # can remove this artifact, not an arbitrary best-fit scale.
    f=np.eye(3,dtype=complex)*(2e-9+3e-9j)
    scalar_ordered=-f/2
    assert abs(task.assess_poles(scalar_ordered,f)[
        'common_pole_residual_relative']-.5)<1e-14
    assert task.assess_poles(2*scalar_ordered,f)[
        'common_pole_cancellation_verified_in_raw_projected_sectors']


def test_source_extractor_convention_is_checked_not_assumed(tmp_path):
    path=tmp_path/'Lagrangian'/'RunMatching.wl'
    path.parent.mkdir(parents=True)
    doc="""(* Matchete's ordered coefficient A5; C5^(ij) = A5^(ij) + A5^(ji) *)
ExtractWeinbergCoefficient[eft_] := Module[{},
 holomorphicTerms = Select[terms, !FreeQ[#, Proj[-1]] &];
 coefficient = CompactWeinbergCoefficient[holomorphicTerms];
 <|"Coefficient" -> coefficient,
 "CanonicalCoefficient" -> canonicalCoefficient|>]
"""
    path.write_text(doc,encoding='utf-8')
    report=task.verify_scalar_extraction_convention(tmp_path)
    assert report['factor_from_saved_batch8_Asym_to_C5']==2
    assert report['ordered_matchete_extractor_returns_A5']
    path.write_text(doc.replace('C5^(ij) = A5^(ij) + A5^(ji)',
                                'C5^(ij) = A5^(ij)'),encoding='utf-8')
    with pytest.raises(task.PoleGateError,match='Unsupported'):
        task.verify_scalar_extraction_convention(tmp_path)
    path.unlink()
    with pytest.raises(task.PoleGateError,match='without'):
        task.verify_scalar_extraction_convention(tmp_path)


def test_raw_pole_success_and_wrong_sign_rejected():
    p=np.eye(3,dtype=complex)*(2e-9+3e-9j)
    assert task.assess_poles(p,-p)['common_pole_cancellation_verified_in_raw_projected_sectors']
    assert not task.assess_poles(p,p)['common_pole_cancellation_verified_in_raw_projected_sectors']
    assert not task.assess_poles(p,-0.97*p)['common_pole_cancellation_verified_in_raw_projected_sectors']


def test_source_polynomial_values_preserve_complex_bars():
    U=sp.Symbol('T3MEU11'); b=sp.Symbol('bar_T3MEU11')
    a=sp.Symbol('lambdaT3');ba=sp.Symbol('bar_lambdaT3')
    source={'pole':(U*b*a, U**2*ba, ba),
            'log':(U*b*a,U**2*ba,ba),
            'constant':(U*b*a,U**2*ba,ba)}
    result=task.evaluate_polynomials(source,inp())
    assert result['pole'][0]==1.0+0.1j
    assert result['pole'][1]==1.0-0.1j
    assert result['pole'][2]==1.0-0.1j


def test_unknown_spurion_blocks():
    z=sp.Symbol('NOT_NATIVE')
    source={name:(z,sp.Integer(1),sp.Integer(1))
            for name in ('pole','log','constant')}
    with pytest.raises(task.PoleGateError,match='Unbound'):
        task.evaluate_polynomials(source,inp())


def test_missing_or_guessed_quartic_map_blocks():
    bad=inp();bad['couplings']={'lambdaT3':1}
    source={name:(sp.Integer(1),)*3 for name in ('pole','log','constant')}
    with pytest.raises(task.PoleGateError,match='quartics'):
        task.evaluate_polynomials(source,bad)


def test_missing_archive_blocks(tmp_path):
    with pytest.raises(task.PoleGateError,match='ONE-LOOP'):
        task.load_native_f_polynomials(tmp_path,tmp_path/'no_archive.inputform.txt','B_-1')


def test_decode_rejects_wrong_shape_and_non_symmetric():
    with pytest.raises(task.PoleGateError,match='shape'):
        task.decode_matrix({'real':[[1]],'imag':[[0]]},'test')
    a=np.eye(3,dtype=complex);a[1,0]=1
    with pytest.raises(task.PoleGateError,match='not symmetric'):
        task.decode_matrix(task.encode_complex(a),'test')


def test_triple_matrix_requires_valid_mf():
    source={name:np.ones(3,dtype=complex) for name in ('pole','log','constant')}
    x=inp();x['mf'][2]=-1
    with pytest.raises(task.PoleGateError,match='Bad MF'):
        task.f_inherited_matrices(x,source,1_000.)


def test_expected_archived_f_source_path():
    assert task.archived_source_path(Path('archive'),'C_-1').name == (
        'C_m1_final_eft1_fermion_first_oneloop_raw_hard_llhh_sector.inputform.txt')
    with pytest.raises(task.PoleGateError):
        task.archived_source_path(Path('archive'),'D_+2')


def test_readonly_native_gates(monkeypatch,tmp_path):
    r=tmp_path/'hard'/'mixed_hard_08'; f=tmp_path/'fixed'/'mixed_fixed_08'
    d=r/'B_m1';d.mkdir(parents=True)
    (r/'component_gate.json').write_text(json.dumps({
        'component_equivalence_verified':True,'full_pole_cancellation_certified':False,
        'physical_C5_certified':False}))
    (d/'hard_scalar_report.json').write_text(json.dumps({
        'component_equivalent_native_reconstruction_used':True,
        'full_EFT1_CheckLagrangian_passed':True,
        'independent_three_mass_native_matching':False}))
    (d/'hard_scalar_reconstruction_audit.json').write_text(json.dumps({
        'component_equivalence_proven':True,'independent_matchete_check_passed':True}))
    source=d/'hard_c5_candidate.txt';source.write_text('native-source')
    rr=tmp_path/'running.json';pp=tmp_path/'params.json'
    rr.write_text('{"model":"B_-1"}');pp.write_text('{}')
    t=f/'B_m1';t.mkdir(parents=True)
    obj={'model':'B_-1','status':'NativeHardThreeFlavorPoleSeparated_NotInheritedPoleCancelled',
         'quartic_map_is_conditional':True,'physical_C5_certified':False,
         'checks':{'three_distinct_MF_lift_inferred_from_common_MF_symbol':True},
         'hard_scalar_source_sha256':task.digest(source),
         'source_fingerprints':{
             'hard_candidate_sha256':task.digest(source),
             'source_fingerprints':{'running':task.digest(rr),'parameters':task.digest(pp)}}}
    (t/'result.json').write_text(json.dumps(obj))
    archive=tmp_path/'archive';archive.mkdir()
    task.archived_source_path(archive,'B_-1').write_text('original-F-archive')
    monkeypatch.setattr(task.native,'OUTPUT',tmp_path/'hard')
    monkeypatch.setattr(task.batch6,'OUTPUT',tmp_path/'fixed')
    monkeypatch.setattr(task.batch6,'verify_hard',lambda *_:(source,{}))
    monkeypatch.setattr(task.native,'paths',lambda *_:{'running':rr,'parameters':pp})
    assert task.verify_inputs('B_-1','mixed_hard_08','mixed_fixed_08',archive)['archive'].is_file()
    # A changed scalar source must not be accepted as the completed batch8 source.
    source.write_text('changed')
    with pytest.raises(task.PoleGateError,match='source changed'):
        task.verify_inputs('B_-1','mixed_hard_08','mixed_fixed_08',archive)
    source.write_text('native-source')
    (d/'hard_scalar_reconstruction_audit.json').write_text(json.dumps({
        'component_equivalence_proven':True,'independent_matchete_check_passed':False}))
    with pytest.raises(task.PoleGateError,match='not independently checked'):
        task.verify_inputs('B_-1','mixed_hard_08','mixed_fixed_08',archive)


def test_no_f_rematch_in_new_script():
    text=Path(task.__file__).read_text()
    assert 'Match[tree' not in text
    assert 'subprocess.run' not in text
    assert 'fermion_Match_recomputed' in text
    assert 'physical_C5_certified' in text

@pytest.mark.parametrize('wrong_pole',[False,True])
def test_calculation_exports_combined_only_after_independent_pole_gate(monkeypatch,tmp_path,wrong_pole):
    x=inp(); x['muF']=150_000.
    v={name:np.asarray([1+0.2j,2-0.3j,1+0.4j])
       for name in ('pole','log','constant')}
    fp,ff,_=task.f_inherited_matrices(x,v,8000.)
    # Batch8 stores ordered scalar Asym; the required normalised scalar
    # C5 is 2*Asym. The F projector is ALREADY half-normalized C5.
    hp=0.5*fp if wrong_pole else -0.5*fp
    h0=np.eye(3,dtype=complex)*(1e-9+1e-9j)
    def hard(*args,muS,**kwargs):
        fin=h0-fp*np.log(muS/8000.) if not wrong_pole else h0
        return {'hard_pole_residue_c5_gev_inverse':task.encode_complex(hp),
                'hard_finite_c5_gev_inverse':task.encode_complex(fin)}
    source=tmp_path/'source.txt';source.write_text('test')
    archive=tmp_path/'f_archive.txt';archive.write_text('original-legacy-source')
    monkeypatch.setattr(task,'load_native_f_polynomials',lambda *a,**kw:{
        'poly':{name:tuple(map(sp.sympify, v[name])) for name in v},
        'archive_path':str(archive),'archive_sha256':task.digest(archive),
        'terms':{'pole':68,'log':68,'constant':68},
        'F_finite_constant_equals_pole':True})
    # Numerical inputs and projectors are source-frozen in the real run.
    monkeypatch.setattr(task,'numerical_inputs',lambda *a,**k:x)
    monkeypatch.setattr(task,'verify_scalar_extraction_convention',lambda *_: {
        'factor_from_saved_batch8_Asym_to_C5':2,
        'ordered_matchete_extractor_returns_A5':True})
    monkeypatch.setattr(task,'evaluate_hard',hard)
    prior={'muS_gev':8000.,'hard_pole_residue_c5_gev_inverse':task.encode_complex(hp),
           'hard_finite_c5_gev_inverse':task.encode_complex(h0)}
    evidence={'hard_source':source,'fixed_result':prior,'archive':archive,
              'running':{},'parameters':{},'native_gate_sha256':'x',
              'fixed_result_sha256':'y',
              'archive_evidence':{'provenance':'OriginalFermionFirstOneLoopArchive',
                'source_sha256':task.digest(archive)}}
    result=task.calculate_one('B_-1',evidence,tmp_path)
    assert result['physical_C5_certified'] is False
    if wrong_pole:
        assert result['status']=='ProjectedRawPoleOrScaleMismatch_NotCombined'
        assert 'conditional_common_scale_sum_c5_gev_inverse' not in result
    else:
        assert result['status']==task.STATUS
        assert result['source_pole_check']['common_pole_cancellation_verified_in_raw_projected_sectors']
        assert result['common_scale_finite_matrix_relative_change']<2e-8
        assert len(result['conditional_common_scale_sum_singular_values_ev'])==3
        assert task.relerr(task.decode_matrix(result[
            'hard_scalar_pole_residue_c5_gev_inverse'],'normalised'),2*hp)<1e-13
        assert abs(result['source_pole_check'][
            'ordered_batch8_unconverted_relative_residual']-0.5)<1e-13
        assert task.relerr(task.decode_matrix(result[
            'conditional_raw_hard_finite_c5_gev_inverse'],'C5'),2*h0)<1e-13


@pytest.mark.parametrize('model,expected', [('B_-1',20),('C_-1',24)])
def test_certified_native_recovery_passes_when_source_and_report_hashes_match(tmp_path, model, expected, monkeypatch):
    d=task.recovered_native_f_archive_paths(tmp_path,model)
    d['source'].parent.mkdir(parents=True)
    d['source'].write_text('one loop test source; not a Matchete expression',encoding='utf-8')
    stem=task.native.tag(model)
    report={
        'format':'t3-mass-eigen-fermion-first-eft1-oneloop-v1',
        'status':'FermionFirstOneLoopRawSectorsExtractedNotRenormalized',
        'model':model,'matching_approach':'FermionFirstTreeAndOneLoopScalarsRemainLight',
        'matching_EFT_order':5,'matching_loop_order':1,'source_quartic_count':17,
        'source_yukawa_count':2,'integrated_heavy_fields':['NewFermion'],
        'scalar_light_fields':['NewScalarX','NewScalarY'],'llss_term_count':expected,
        'raw_hard_llhh_candidate_term_count':204,'raw_hard_llhh_file':d['source'].name,
        'tree_matching_completed':True,'one_loop_matching_completed':True,
        'loop_difference_fermion_absent':True,
        'raw_hard_llhh_scheme_renormalized':False,
        'raw_hard_llhh_pole_rge_certified':False,'physical_C5_computed':False}
    d['report'].write_text(json.dumps(report),encoding='utf-8')
    inspection={
        'status':'MATCHeTE_RAW_MIXED_LOOP_LLHH_EXTRACTED_NOT_MSbar',
        'model':model,'driver_type':'mixed_BC','native_exit_code':0,
        'one_loop_F_native_run_returned':True,'hard_sector_renormalized':False,
        'explicit_hard_LLHH_candidate_terms':204,'physical_wilson_coefficient':None,
        'raw_hard_llhh_file':d['source'].name,
        'raw_hard_llhh_sha256':task.digest(d['source']),
        'native_report_sha256':task.digest(d['report'])}
    d['inspection'].write_text(json.dumps(inspection),encoding='utf-8')
    monkeypatch.setitem(task.KNOWN_NATIVE_F_SHA256,model,task.digest(d['source']))
    evidence=task.inspect_recovered_native_archive(tmp_path,model)
    original_archive_root=task.ARCHIVE
    monkeypatch.setattr(task,'ARCHIVE',tmp_path/'legacy_archive')
    monkeypatch.setattr(task,'CERTIFIED_NATIVE_F_ARCHIVE',tmp_path)
    recovered=task.resolve_f_first_archive(task.ARCHIVE,model)
    assert recovered['archive']==d['source']
    assert evidence['source_sha256']==task.digest(d['source'])
    assert evidence['provenance']=='ExistingNativeSpinorCertificationFOneLoop'
    # A corrupt source is rejected even if the JSON report still claims success.
    d['source'].write_text('tampered',encoding='utf-8')
    with pytest.raises(task.PoleGateError,match='SHA256'):
        task.inspect_recovered_native_archive(tmp_path,model)
    d['source'].write_text('one loop test source; not a Matchete expression',encoding='utf-8')
    # A report claiming renormalized matching cannot be silently treated as raw.
    report['raw_hard_llhh_scheme_renormalized']=True
    d['report'].write_text(json.dumps(report),encoding='utf-8')
    inspection['native_report_sha256']=task.digest(d['report'])
    d['inspection'].write_text(json.dumps(inspection),encoding='utf-8')
    with pytest.raises(task.PoleGateError,match='report'):
        task.inspect_recovered_native_archive(tmp_path,model)


def test_recovery_resolver_prioritizes_primary_and_refuses_explicit_root(monkeypatch,tmp_path):
    archive=tmp_path/'original';archive.mkdir()
    primary=task.archived_source_path(archive,'B_-1')
    primary.write_text('archive',encoding='utf-8')
    assert task.resolve_f_first_archive(archive,'B_-1')['archive']==primary
    primary.unlink()
    with pytest.raises(task.PoleGateError,match='explicitly selected'):
        task.resolve_f_first_archive(archive,'B_-1')


def test_native_archive_no_fermion_rematch_in_run():
    source=Path(task.__file__).read_text()
    assert 'RunT3Matching[' not in source
    assert 'wolframscript' not in source.lower()
    assert 'source_one_loop_pole_rge_certified' in source
