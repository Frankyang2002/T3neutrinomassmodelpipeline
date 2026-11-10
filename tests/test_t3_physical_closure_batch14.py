"""Standalone, source-format and mathematically independent three-gate checks."""
from __future__ import annotations

import importlib.util
from pathlib import Path
import json
import numpy as np
import pytest

PATH=Path(__file__).resolve().parents[1]/'Numerical/parallel_t3/complete/two_model_physical_closure_batch14.py'
spec=importlib.util.spec_from_file_location('batch14_audit',PATH)
audit=importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


def test_normalization_is_indeterminate_from_same_operator_logs():
    z=np.array([1.+2j,2.-1j,-3j])
    for factor in (1.,-1.,2.,-.5):
        result=audit.inferred_real_normalization(z,factor*z)
        assert result['factor']==pytest.approx(factor)
        assert result['residual']<1e-14
        assert result['absolute_gamma_CC_to_Weyl_factor_certified'] is False


def test_normalization_zero_log_cannot_be_fit():
    r=audit.inferred_real_normalization(np.zeros(4),np.ones(4))
    assert r['status']=='Unidentifiable_ZeroMasterLog'


def test_bad_matrix_decoding_fails_closed():
    z={'real':np.eye(3).tolist(),'imag':np.zeros((3,3)).tolist()}
    assert np.array_equal(audit.decoded_matrix(z),np.eye(3))
    with pytest.raises(ValueError,match='3x3'):
        audit.decoded_matrix({'real':[[1,2]],'imag':[[1,2]]})


def test_abelian_operator_ward_satisfied_then_violated():
    # Weyl charges q1=+1 and q2=-1, no charged scalars.
    t=np.diag([1.,-1.]).astype(complex)
    s=np.zeros((2,2),complex)
    c=np.zeros((2,2,2,2),complex)
    c[0,1,0,0]=c[1,0,0,0]=1.
    assert audit.operator_ward(c,t,s)==pytest.approx(0.)
    c[0,0,0,0]=1.
    assert audit.operator_ward(c,t,s)==pytest.approx(2.)


def test_abelian_scalar_quartic_ward_invariant():
    # Two real components with SO(2) generator i*antisymmetric.
    theta=1j*np.array([[0,1],[-1,0]],complex)
    delta=np.eye(2)
    lam=(np.einsum('ab,cd->abcd',delta,delta)
         +np.einsum('ac,bd->abcd',delta,delta)
         +np.einsum('ad,bc->abcd',delta,delta)).astype(complex)
    assert audit.scalar_quartic_ward(lam,theta)<1e-14


def test_abelian_yukawa_ward_invariant():
    t=np.diag([1.,-1.]).astype(complex)
    s=np.zeros((1,1),complex)
    y=np.zeros((2,2,1),complex)
    y[0,1,0]=y[1,0,0]=1
    assert audit.yukawa_ward(y,t,s)==pytest.approx(0.)
    y[0,0,0]=1
    assert audit.yukawa_ward(y,t,s)>0


def test_chirality_inventory_requires_three_each_A():
    entries=[]
    for c in ('Proj[1]','Proj[-1]'):
        entries += [{'sector':'CandidateLeptonNumberViolatingLLHH',
                     'raw_inputform':'NCM[GammaCC,'+c+']'}]*3
    root={'terms':entries}
    example={'model':'A_-4','master_candidate_EFT1_log_GeVinv':{
          'real':np.eye(3).tolist(),'imag':np.zeros((3,3)).tolist()},
             'native_B11_pole_derived_EFT1_log_GeVinv':{
          'real':np.eye(3).tolist(),'imag':np.zeros((3,3)).tolist()}}
    out=audit.audit_normalization(root,example)
    assert out['native_GammaCC_operator_chiralities']=={'PR':3,'PL':3}
    assert not out['gate_passed']


def test_source_matching_counts_not_just_LLHH():
    model='A_-4'
    sectors={'OtherEFT1Operator':81,
             'CandidateLeptonNumberViolatingLLScalarScalar':48,
             'CandidateLeptonNumberViolatingLLHH':6}
    ledger={'status':'UnrenormalisedOneLoopTermsEnumerated','sector_counts':sectors,
            'terms':[{'sector':k} for k,n in sectors.items() for _ in range(n)]}
    matrix=lambda x:{'real':(np.eye(3)*x).tolist(),'imag':np.zeros((3,3)).tolist()}
    b11={'source_hard_symmetric_pole_GeVinv':matrix(1),
         'source_soft_symmetric_pole_GeVinv':matrix(-1),
         'source_hard_symmetric_regular_GeVinv':matrix(2),
         'source_soft_symmetric_regular_GeVinv':matrix(3),
         'source_EFT1_scale_partition_log_GeVinv':matrix(4),
         'source_finite_sector_candidate_GeVinv':matrix(9),
         'native_pole_log_identity_relative':0.,
         'source_scale_partition_invariance_relative':0.}
    st={'A_physical_scalar_soft_Weinberg_terms':10}
    res=audit.audit_matching(ledger,b11,st,model,'hbar ϵ μbar2 T3EFT1W12')
    assert res['full_F_one_loop_operator_count']==135
    assert res['F_LNV_LLSS_one_loop_operator_count_not_dropped_from_ledger']==48
    assert not res['gate_passed']


def test_finite_matching_must_recompose():
    model='A_-4'
    sectors={'OtherEFT1Operator':81,
             'CandidateLeptonNumberViolatingLLScalarScalar':48,
             'CandidateLeptonNumberViolatingLLHH':6}
    ledger={'status':'UnrenormalisedOneLoopTermsEnumerated','sector_counts':sectors,
            'terms':[{'sector':k} for k,n in sectors.items() for _ in range(n)]}
    m=lambda x:{'real':(np.eye(3)*x).tolist(),'imag':np.zeros((3,3)).tolist()}
    b11={'source_hard_symmetric_pole_GeVinv':m(1),
         'source_soft_symmetric_pole_GeVinv':m(-1),
         'source_hard_symmetric_regular_GeVinv':m(2),
         'source_soft_symmetric_regular_GeVinv':m(3),
         'source_EFT1_scale_partition_log_GeVinv':m(4),
         'source_finite_sector_candidate_GeVinv':m(8),
         'native_pole_log_identity_relative':0.,
         'source_scale_partition_invariance_relative':0.}
    with pytest.raises(ValueError,match='recomposition'):
        audit.audit_matching(ledger,b11,{'A_physical_scalar_soft_Weinberg_terms':10},
                             model,'hbar ϵ μbar2 T3EFT1W12')


def test_missing_project_files_skip_integration():
    root=Path(__file__).resolve().parents[1]
    paths=audit.preflight_paths(root,'A_-4')
    if not all(p.exists() for p in paths):
        pytest.skip('Requires saved full native Batch11/13 outputs in Windows project')
    assert paths[0].is_file()
