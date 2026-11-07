"""Synthetic tests only: source-gated production comparison requires local repo."""
from __future__ import annotations
import importlib.util
import math
from pathlib import Path
import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[2]
SPEC=importlib.util.spec_from_file_location('t3_batch10_v3',ROOT/'scripts/CompareT3BCBatch9ToTriangle.py')
assert SPEC and SPEC.loader
b=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(b)

def test_majorana_transpose_matrix_roundtrip():
    x=np.array([[2j,1+2j,0],[1+2j,-4j,2],[0,2,1j]],complex)
    assert b.symmetry_gate(x,'test')==0
    assert b.relative_frobenius(x,b.decode_matrix(b.encode_matrix(x),'test'))==0
    with pytest.raises(b.ComparisonError):b.symmetry_gate(np.triu(x),'not symmetric')
    assert b.relative_frobenius(x,2*x)>.4

def test_independent_log_intercept_removes_exact_t_minus_two_correction():
    # Test includes both t^-2 constant and logarithm, the source of the v2 failure.
    a=np.array([[1,2j,3],[2j,4,5j],[3,5j,6]],complex)
    c=np.eye(3,dtype=complex)
    d=np.full((3,3),0.05j)
    d=(d+d.T)/2
    e=np.eye(3,dtype=complex)*(-0.09+0.3j)
    f=lambda t:a+c*math.log(t)+(d+e*math.log(t))/(t*t)
    values={float(t):f(float(t)) for t in (128,256,512,1024,2048)}
    a1,b1=b.fit_log_plus_power((128.,256.,512.,1024.),values)
    a2,b2=b.fit_log_plus_power((256.,512.,1024.,2048.),values)
    assert b.relative_frobenius(a1,a)<1e-11
    assert b.relative_frobenius(a2,a)<1e-11
    assert b.relative_frobenius(b2,c)<1e-11
    assert b.relative_frobenius(b1,b2)<1e-11


def test_subleading_log_t4_converges_not_fitted_to_target():
    a=np.eye(3,dtype=complex)*(0.7+0.2j)
    slope=np.eye(3,dtype=complex)*(-0.3+0.1j)
    correction=np.eye(3,dtype=complex)*(2+3j)
    remnant=np.eye(3,dtype=complex)*(4-5j)
    fun=lambda t:a+slope*math.log(t)+correction*math.log(t)/(t*t) + remnant*math.log(t)/(t**4)
    values={float(t):fun(float(t)) for t in (64,128,256,512,1024)}
    x1,s1=b.fit_log_plus_power((64.,128.,256.,512.),values)
    x2,s2=b.fit_log_plus_power((128.,256.,512.,1024.),values)
    assert b.relative_frobenius(x2,a)<b.relative_frobenius(x1,a)
    assert b.relative_frobenius(s2,slope)<b.relative_frobenius(s1,slope)
    assert b.relative_frobenius(x2,a)<1e-6

def test_source_ordered_A5_to_C5_no_factor_fit():
    masses=np.array([100.,110.,120.])
    y1=np.array([[1+1j,2,3],[4,5+1j,6],[7,8,9]],complex)/10
    y2=np.array([[1j,2j,3j],[2,3,4],[1,2,3]],complex)/20
    inp={'ms':np.array([10.,12.]),'U':np.eye(2),
         'couplings':{k:1. for k in b.Q_NAMES},'y1':y1,'y2':y2}
    def ordered(M,MX,MY,U,Q):
        return (1/M,2j/M,3/M)
    actual=b.assemble_direct(masses,inp,ordered)
    expected=np.zeros((3,3),complex)
    for k,M in enumerate(masses):
        L=np.column_stack((y1[:,k].conj(),y2[:,k].conj()))
        K=np.array([[2,2j],[2j,6]],complex)/M
        expected+=L@K@L.T
    np.testing.assert_allclose(actual,expected,atol=1e-17,rtol=1e-13)
    assert b.relative_frobenius(actual,2*expected)>.4

def test_asymptotic_matrix_from_synthetic_coefficient():
    masses=np.array([100.,110.,120.]); y1=np.eye(3,dtype=complex)
    inp={'mf':masses,'ms':np.array([10.,12.]),'U':np.eye(2,dtype=complex),
         'couplings':{k:1. for k in b.Q_NAMES},'y1':y1,'y2':np.eye(3)*0.1j}
    def ordered(M,MX,MY,U,Q):
        return ((1+math.log(M/100))/M, 0j,2/M)
    full=b.assemble_direct(masses,inp,ordered)
    result=b.leading_from_direct(inp,ordered,128)
    np.testing.assert_allclose(result['leading'],full,rtol=1e-11,atol=1e-12)
    assert result['asymptotic_matrix_convergence_relative']<1e-11
    assert b.symmetry_gate(result['leading'],'leading')<1e-12
