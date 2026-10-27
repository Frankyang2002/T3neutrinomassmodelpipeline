"""Offline checks on certified two-threshold symbolic export and numeric converter."""
from __future__ import annotations

import json
from pathlib import Path
import shutil
import zipfile

import pytest

from scripts.RunGeneralT3TwoThresholdNumericalC5 import (
    H_PATTERNS, EPS, OUT_DIR, ASSEMBLY_DIR, ERROR,
    project_expression, project_holomorphic_term, projection_report,
    numerical_evaluation, stem,
)
from scripts.RunGeneralT3TwoThresholdPoleAssembly import split_terms

HERE=Path(__file__).resolve().parents[2]


@pytest.mark.parametrize('model', ['B_-1','C_-1'])
def test_real_uploaded_projection(model,tmp_path):
    report, coef = projection_report(model, ASSEMBLY_DIR, tmp_path)
    assert report['PL_holomorphic_term_count']==116
    assert report['PR_hermitian_conjugate_term_count']==116
    assert report['upstream_pole_cancellation_verified'] is True
    assert report['physical_C5_computed'] is False
    assert report['full_EFT1_running_performed'] is False
    text=coef.read_text()
    assert 'Coupling[y1' in text and 'Coupling[y2' in text
    assert 'GammaCC' not in text and 'Field[H' not in text and 'CG[' not in text
    assert len(split_terms(text))==116


@pytest.mark.parametrize('model,orientations',[
    ('B_-1',{'ordinary':116}),
    ('C_-1',{'ordinary':34,'double_epsilon_swap':82}),
])
def test_exact_operator_orientation_counts(model, orientations):
    p=ASSEMBLY_DIR/(stem(model)+'_two_threshold_raw_llhh_fixed_muF_unprojected.inputform.txt')
    result=project_expression(p.read_text(),model)
    assert result['orientations']==orientations


def first_PL(model):
    p=ASSEMBLY_DIR/(stem(model)+'_two_threshold_raw_llhh_fixed_muF_unprojected.inputform.txt')
    return next(t for _,t in split_terms(p.read_text()) if 'Proj[-1]' in t)


@pytest.mark.parametrize('tamper',[
    lambda x:x.replace('Proj[-1]','Proj[1]',1),
    lambda x:x.replace('GammaCC','Gamma5',1),
    lambda x:x.replace('Field[H, Scalar','Field[NewScalarX, Scalar',1),
    lambda x:x.replace('CG[Bar[eps[SU2L]]','CG[eps[SU2L]',1),
    lambda x:x.replace('Index[d$$1, NFlavor]', 'Index[d$$2, NFlavor]',1),
    lambda x:x.replace('Index[d$$1, Flavor]', 'Index[d$$8, Flavor]',1),
    lambda x:x.replace('Coupling[MF, {}, 0]','Coupling[Unknown, {}, 0]',1),
])
def test_projector_refuses_corrupt_operator(tamper):
    raw=first_PL('B_-1')
    with pytest.raises(ValueError):project_holomorphic_term(tamper(raw))


def test_projection_rejects_corrupt_source(tmp_path):
    shutil.copytree(ASSEMBLY_DIR,tmp_path/'assembly')
    p=tmp_path/'assembly/B_m1_final_two_threshold_raw_llhh_fixed_muF_unprojected.inputform.txt'
    p.write_text(p.read_text().replace('Field[H, Scalar','Field[NewScalarX, Scalar',1))
    with pytest.raises(ValueError):projection_report('B_-1',tmp_path/'assembly',tmp_path/'out')


def test_projection_rejects_missing_pole_certificate(tmp_path):
    shutil.copytree(ASSEMBLY_DIR,tmp_path/'assembly')
    p=tmp_path/'assembly/B_m1_final_two_threshold_raw_llhh_audit.json'
    a=json.loads(p.read_text());a['matched_canonical_pole_term_count']=67;p.write_text(json.dumps(a))
    with pytest.raises(ValueError,match='certificate'):
        projection_report('B_-1',tmp_path/'assembly',tmp_path/'out')


def test_projection_refuses_wrong_model(tmp_path):
    with pytest.raises(ValueError):projection_report('A_-1',ASSEMBLY_DIR,tmp_path)


def test_numeric_requires_real_reference(tmp_path):
    report,p=projection_report('B_-1',ASSEMBLY_DIR,tmp_path)
    with pytest.raises(ValueError,match='frozen-aligned'):
        numerical_evaluation('B_-1',p.read_text(),tmp_path/'no-file.json',tmp_path,report)


def test_numeric_rejects_fabricated_certification(tmp_path):
    report,p=projection_report('B_-1',ASSEMBLY_DIR,tmp_path)
    ref=tmp_path/'fake.json'
    ref.write_text(json.dumps({'model':'B_-1','status':'Success','physical_c5_certified':True}))
    with pytest.raises(ValueError,match='Incompatible'):
        numerical_evaluation('B_-1',p.read_text(),ref,tmp_path,report)


@pytest.mark.parametrize('model',['B_-1','C_-1'])
def test_synthetic_aligned_benchmark_smoke(model,tmp_path):
    # Regression toy ONLY: no fitted-reference inputs or claimed physics.
    def c(x):return {'real':x,'imag':0}
    y1=[[c(0.005)],[c(0.003)],[c(0.004)]]
    y2=[[c(0.004)],[c(0.005)],[c(0.006)]]
    old={
        'model':model,
        'status':'FrozenScalarDiagonalDirectOneScaleC5NotPhysicalEFT1',
        'frozen_scalar_aligned_one_scale_scenario':True,
        'mass_basis_compatible_with_frozen_reference':True,
        'physical_c5_certified':False,'sequential_EFT1_performed':False,
        'source_mass_parameters_used_as_numeric_scalar_input':False,
        'fermion_masses_gev':[100000],
        'scalar_masses_gev':[1000,1142.857142857],
        'effective_scalar_basis_U':{'real':[[0,1],[1,0]],'imag':[[0,0],[0,0]]},
        'benchmark_config':{
            'format':'t3-frozen-scalar-diagonal-direct-C5-input-v1',
            'vev_gev':246,'y1':y1,'y2':y2,
            'quartics':{'lambdaT3':c(0.1),
               'lambdaT3ExtraQ020020I1':c(0),
               'lambdaT3ExtraQ020200I1':c(0)}},
        'takagi_singular_values_eV_descending':[1.0,0.5,0.1]}
    ref=tmp_path/'toy-regression-only.json';ref.write_text(json.dumps(old))
    report,p=projection_report(model,ASSEMBLY_DIR,tmp_path)
    result=numerical_evaluation(model,p.read_text(),ref,tmp_path,report)
    assert result['largest_takagi_mass_eV']>0
    assert result['renormalized_physical_C5_certified'] is False
    assert result['full_EFT1_running_performed'] is False
    assert len(result['mnu_gev'])==3
    assert len(result['takagi_singular_values_eV_descending'])==3


@pytest.mark.parametrize('model',['B_-1','C_-1'])
def test_conjugate_operator_tamper(model):
    from scripts.RunGeneralT3TwoThresholdNumericalC5 import check_conjugate_operator
    p=ASSEMBLY_DIR/(stem(model)+'_two_threshold_raw_llhh_fixed_muF_unprojected.inputform.txt')
    part=next(t for _,t in split_terms(p.read_text()) if 'Proj[1]' in t)
    assert check_conjugate_operator(part) in {'ordinary','double_epsilon_swap'}
    with pytest.raises(ValueError):
        check_conjugate_operator(part.replace('GammaCC','Gamma5',1))
    with pytest.raises(ValueError):
        check_conjugate_operator(part.replace('eps[SU2L]','eps[SU3c]',1))
