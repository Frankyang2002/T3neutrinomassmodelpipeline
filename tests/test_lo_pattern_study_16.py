"""Restricted T3 pattern-study tests. Synthetic tests do NOT certify native physics."""
from __future__ import annotations

import copy
import json
import math
import os
from pathlib import Path
import tempfile

import numpy as np
import pytest

from Numerical.direct_one_loop.core import encoded, evaluate, b0_finite, takagi
from Numerical.lo_certification.direct_uv.native_uv_loop_reference import (
    independent_reference, decode, rephase_input,
)
from Numerical.lo_certification.direct_uv.native_amplitude_contract import model_dictionary
from Numerical.lo_patterns_16.pattern_study import (
    StudyError, _simple_expr, _scans, _compare, transformed_input,
    audit_weight_lookup, _validate_dictionary, run_study,
)

ROOT = Path(__file__).resolve().parents[1]


def dirac_toy():
    M=np.array([[0,150],[150,0]],complex)
    A=np.array([[40000,2000],[2000,52000]],complex)
    y1=np.zeros((3,2,2),complex);y2=np.zeros_like(y1)
    y1[:,0,0]=[.1,.2,.3]
    y2[:,1,1]=[.2,.1,-.4]
    return {'model':'synthetic','fermion_mass_gev':encoded(M),
            'scalar_A_gev2':encoded(A),'scalar_B_gev2':encoded(np.zeros((2,2))),
            'renormalization_scale_gev':1e7,
            'vertices':[{'name':'y1_HC','conjugate':encoded(y1)},
                        {'name':'y2_HC','field':encoded(y2)}]}


def ma_toy():
    y=.13; a=42000.; b=700.; mf=100.
    return {'model':'synthetic_Majorana','fermion_mass_gev':encoded([[mf]]),
            'scalar_A_gev2':encoded([[a]]),'scalar_B_gev2':encoded([[b]]),
            'renormalization_scale_gev':1e7,
            'vertices':[{'name':'y_HC','field':encoded(np.array([[[y]]]*3))}]}


def mat(data):
    return decode(evaluate(data,retain_individual_terms=False)['total']['complex_symmetric_matrix_gev'],'mass',(3,3))


def test_16_catalogue_two_Majoranas_and_six_additional_Yukawa_options():
    x=model_dictionary()
    assert len(x)==len(set(z['model'] for z in x))==16
    assert [z['model'] for z in x if z['fermion_type']=='Majorana']==['B_-1','C_-1']
    assert sum(z['extra_symmetry_allowed_yukawa'] is not None for z in x)==6


def test_exact_audit_weight_complex_phase():
    assert abs(_simple_expr('-sqrt(2)*3**(3/4)*I/3')+1j*math.sqrt(2)*3**(.75)/3)<1e-14
    assert abs(_simple_expr('2*sqrt(6)/3')-2*math.sqrt(6)/3)<1e-14
    with pytest.raises(StudyError):_simple_expr("__import__('os').system('echo bad')")


def test_source_group_weight_perturbation_rejected():
    models={m['model']:m for m in model_dictionary()}
    audit_path=ROOT/'output/general_observable/lo_direct_one_loop/source_audits/charge_resolved_loop_weights.json'
    if not audit_path.is_file():pytest.skip('Source audit unavailable in this checkout')
    audit=json.loads(audit_path.read_text(encoding='utf-8-sig'))
    assert len(audit_weight_lookup(audit,models))==16
    modified=copy.deepcopy(audit)
    modified['contraction']='y1 * y2 * T3MixCG' # wrong conjugation convention
    with pytest.raises(StudyError,match='contraction'):audit_weight_lookup(modified,models)
    modified=copy.deepcopy(audit)
    modified['models'][0]['model']['dF']=99 # source drift
    with pytest.raises(StudyError,match='representation'):audit_weight_lookup(modified,models)


def test_gaussian_majorana_one_loop_factor_and_no_double_count():
    d=ma_toy();ref=independent_reference(d)
    z=decode(ref['complex_symmetric_matrix_gev'],'mass',(3,3))
    a=42000.;b=700.;mf=100.;mu=1e7;y=.13
    predicted=y*y*mf/(32*math.pi**2)*(
        b0_finite(mf*mf,a+b,mu**2)-b0_finite(mf*mf,a-b,mu**2))
    assert np.allclose(z,predicted*np.ones((3,3)),rtol=1e-9,atol=1e-15)
    assert abs(predicted)>1e-10
    assert not np.allclose(z,2*predicted*np.ones((3,3)),rtol=1e-5,atol=1e-10)
    assert _compare(evaluate(d),ref)['relative_matrix_error']<1e-8


def test_yukawa_exact_linear_and_lambda_small_mixing():
    d=dirac_toy();m0=mat(d)
    assert np.linalg.norm(m0)>1e-12
    for mode in ('y1','y2'):
        for t in (-1.,0.,.5,2.):
            m=mat(transformed_input(d,mode,t,1))
            assert np.linalg.norm(m-t*m0)<2e-10*np.linalg.norm(m0)
    z=mat(transformed_input(d,'lambdaT3',0.,1))
    assert np.linalg.norm(z)<1e-9*np.linalg.norm(m0)
    minus=mat(transformed_input(d,'lambdaT3',-1.,1))
    assert np.linalg.norm(minus+m0)<2e-10*np.linalg.norm(m0)
    small=copy.deepcopy(d)
    a=decode(small['scalar_A_gev2'],'A',(2,2));a[0,1]=a[1,0]=20.
    small['scalar_A_gev2']=encoded(a)
    ms=mat(small)
    mhalf=mat(transformed_input(small,'lambdaT3',.5,1))
    assert np.linalg.norm(mhalf-.5*ms)/np.linalg.norm(ms)<1e-5


def test_pole_free_and_takagi_masses_are_singular_values():
    d=dirac_toy();e=evaluate(d,retain_individual_terms=True)
    assert e['diagnostics']['uv_finite']
    assert e['diagnostics']['neutrino_mass_matrix_symmetric']
    m=decode(e['total']['complex_symmetric_matrix_gev'],'m',(3,3))
    assert np.linalg.norm(m-m.T)<1e-10*np.linalg.norm(m)
    s=np.sort(np.linalg.svd(m,compute_uv=False))*1e9
    reported=np.asarray(e['total']['neutrino_masses_ev'],dtype=float)
    # The synthetic Dirac toy has only two independent flavour-Yukawa vectors,
    # so its 3x3 neutrino mass is exactly rank <= 2.  LAPACK's SVD can return
    # a tiny positive singular value for an algebraically zero mode, whereas
    # the production Takagi routine deliberately identifies that mode as zero.
    # Use a matrix-scale floating-point floor ONLY for this expected zero,
    # not an indiscriminate large absolute tolerance on nonzero masses.
    zero_floor_ev=32*np.finfo(float).eps*m.shape[0]*s[-1]
    assert s[0] <= zero_floor_ev
    assert reported.shape==(3,)
    assert np.all(reported>=0.) and np.all(np.diff(reported)>=0.)
    assert reported[0] <= zero_floor_ev
    np.testing.assert_allclose(s[1:],reported[1:],rtol=1e-8,atol=1e-12)
    U,masses,res=takagi(m)
    assert res<1e-8 and masses[0]*1e9 <= zero_floor_ev
    np.testing.assert_allclose(masses[1:]*1e9,s[1:],rtol=1e-8,atol=1e-12)
    assert len(e['individual_terms'])>0


def test_full_rank_takagi_matches_all_singular_values_without_zero_floor():
    # Independent full-rank control: all three masses must be resolved.
    m=np.array([[3.+.2j,.12+.31j,-.04],
                [.12+.31j,2.-.5j,.27+.1j],
                [-.04,.27+.1j,1.+.3j]],dtype=complex)
    U,masses,res=takagi(m)
    s=np.sort(np.linalg.svd(m,compute_uv=False))
    assert s[0] > .5
    assert res<1e-8
    np.testing.assert_allclose(masses,s,rtol=1e-11,atol=1e-12)


def test_rephase_invariance_and_negative_conjugation_control():
    d=dirac_toy()
    original=decode(independent_reference(d)['complex_symmetric_matrix_gev'],'m',(3,3))
    rotated=rephase_input(d,np.exp(1j*np.array([.23,-.47])),np.exp(1j*np.array([.52,-.13])))
    new=decode(independent_reference(rotated)['complex_symmetric_matrix_gev'],'m',(3,3))
    assert np.linalg.norm(original-new)/np.linalg.norm(original)<1e-8
    bad=copy.deepcopy(d)
    bad['vertices'][0]['field']=bad['vertices'][0].pop('conjugate')
    changed=decode(independent_reference(bad)['complex_symmetric_matrix_gev'],'m',(3,3))
    assert np.linalg.norm(original-changed)>1e-8*np.linalg.norm(original)


def test_scalar_tachyon_and_nonfinite_negative_controls():
    d=dirac_toy()
    bad=transformed_input(d,'lambdaT3',1e9,1)
    with pytest.raises((ValueError,RuntimeError)):
        evaluate(bad)
    bad=copy.deepcopy(d)
    bad['fermion_mass_gev']['real'][0][1]=float('nan')
    with pytest.raises(ValueError):evaluate(bad)


def test_scan_report_includes_sign_flip_zero_and_mass_ratio():
    d=dirac_toy();ans=evaluate(d,retain_individual_terms=False)
    rows=_scans(d,ans,1,'synthetic')
    assert len(rows)==18 # 3*5 coupling points + 3 heavy ratios
    assert all(x['status']=='conditional_uv_scanned' for x in rows)
    assert {x['t'] for x in rows if x['mode']=='lambdaT3'}=={-1.,0.,.5,1.,2.}
    assert len([r for r in rows if r['mode']=='heavy_mass'])==3


def test_false_certificate_never_promoted(tmp_path):
    assert _validate_dictionary(None,tmp_path)['physical_source_certified'] is False
    path=tmp_path/'false_certificate.json'
    path.write_text(json.dumps({'physical_source_certified':True}))
    with pytest.raises(StudyError,match='cannot be promoted'):
        _validate_dictionary(path,tmp_path)
    path.write_text(json.dumps({'physical_source_certified':False}))
    with pytest.raises(StudyError,match='schema'):
        _validate_dictionary(path,tmp_path)


def test_output_tree_protected(tmp_path):
    with pytest.raises(StudyError,match='NEW directory'):
        run_study(tmp_path,tmp_path/'data_unrelated')


def test_genuine_frozen_16_on_local_checkout():
    """Opt-in only. Not 16 synthetic models and NOT physically certified."""
    if os.environ.get('RUN_LO_PATTERN_INTEGRATION')!='1':
        pytest.skip('Set RUN_LO_PATTERN_INTEGRATION=1 for genuine frozen input integration')
    owned=ROOT/'output/general_observable/lo_pattern_study16'
    owned.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='pytest_pattern16_',dir=owned) as tdir:
        # Remove empty dir: driver refuses to reuse any nonempty directory.
        report=run_study(ROOT,Path(tdir),do_scans=False,do_plots=False)
        assert len(report['models'])==16
        assert report['counts']=={'catalogue':16,'conditional_evaluated':16,'blocked':0,
                                    'physical_source_certified':0,'scan_passes':0,'scan_blocked':0}
        assert all(row['benchmark_sha256'] and row['component_input_sha256'] for row in report['models'])
        for row in report['models']:
            assert row['m1_ev'] is not None and row['m2_ev'] is not None and row['m3_ev'] is not None


def test_degenerate_internal_fermion_basis_permutation():
    d=dirac_toy()
    reference=mat(d)
    # Two identical Takagi singular values M=150; swap the Weyl source basis.
    permutation=np.array([1,0])
    changed=copy.deepcopy(d)
    M=decode(d['fermion_mass_gev'],'M',(2,2))
    changed['fermion_mass_gev']=encoded(M[np.ix_(permutation,permutation)])
    for row in changed['vertices']:
        for key in ('field','conjugate'):
            if key in row:
                original=decode(row[key],'vertex',(3,2,2))
                row[key]=encoded(original[:,permutation,:])
    new=mat(changed)
    assert np.linalg.norm(new-reference)/np.linalg.norm(reference)<1e-8


def test_msbar_scale_cancels_and_nondegenerate_singular_values():
    d=dirac_toy()
    default=mat(d)
    changed=copy.deepcopy(d)
    changed['renormalization_scale_gev']=2.3e5
    new=mat(changed)
    assert np.linalg.norm(new-default)/np.linalg.norm(default)<1e-8
    assert np.all(np.isfinite(np.linalg.svd(new,compute_uv=False)))


def test_exact_cg_weight_failure_wrong_factor_and_wrong_charge():
    models={m['model']:m for m in model_dictionary()}
    audit_path=ROOT/'output/general_observable/lo_direct_one_loop/source_audits/charge_resolved_loop_weights.json'
    if not audit_path.exists():pytest.skip('Source audit unavailable')
    audit=json.loads(audit_path.read_text(encoding='utf-8-sig'))
    modified=copy.deepcopy(audit)
    modified['models'][0]['channels'][0]['complex_weight']='2*sqrt(3)/3'
    with pytest.raises(StudyError,match='product/complex phase'):
        audit_weight_lookup(modified,models)
    modified=copy.deepcopy(audit)
    modified['models'][0]['channels'][0]['charge']='123'
    with pytest.raises(StudyError,match='charge mismatch'):
        audit_weight_lookup(modified,models)
