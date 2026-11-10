"""Batch 4 native PL Matchete parser, mixing, quartic and strict replacement tests.

Every test uses a source-style Laurent constant kernel; no false assertion
about physical hard/soft subtraction or Matchete-to-Weyl signs is made.
"""
import math
import numpy as np
import pytest

from Numerical.parallel_t3.wilson.native_finite_batch4 import (
    parse_native_kernel,parse_scalar_number,matchete_quartics,evaluate_kernels,
    MF,LOOP)

# Exact form of a tracked A(-4) native PL InputForm constant kernel.
A4 = (
    '+ (2*hbar*Bar[Coupling[y1, {Index[d$$1, Flavor], Index[d$$1, NFlavor]}, 0]]'
    '*Bar[Coupling[y2, {Index[d$$2, Flavor], Index[d$$1, NFlavor]}, 0]]'
    '*CG[Bar[eps[SU2L]], {Bar[Index[d$$1, SU2L[fund]]], Bar[Index[d$$2, SU2L[fund]]]}]'
    '*CG[Bar[eps[SU2L]], {Bar[Index[d$$3, SU2L[fund]]], Bar[Index[d$$4, SU2L[fund]]]}]'
    '*Coupling[lambdaT3, {}, 0]'
    '*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]'
    '*Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]'
    '*NCM[Transp[Field[l, Fermion, {Index[d$$4, SU2L[fund]], Index[d$$2, Flavor]}, {}]], '
    'DiracProduct[GammaCC, Proj[-1]], '
    'Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}]])'
    '/(Sqrt[3]*Coupling[MF, {Index[d$$1, NFlavor]}, 0])'
)
B_MIX = (
    '- (hbar*Bar[Coupling[lambdaT3ExtraQ020200I1, {}, 0]]'
    '*Bar[Coupling[T3MEU11, {}, 0]]^2'
    '*Bar[Coupling[y1, {Index[d$$1, Flavor], Index[d$$1, NFlavor]}, 0]]'
    '*Bar[Coupling[y1, {Index[d$$2, Flavor], Index[d$$1, NFlavor]}, 0]]'
    '*CG[Bar[eps[SU2L]], {Bar[Index[d$$1, SU2L[fund]]], Bar[Index[d$$2, SU2L[fund]]]}]'
    '*CG[Bar[eps[SU2L]], {Bar[Index[d$$3, SU2L[fund]]], Bar[Index[d$$4, SU2L[fund]]]}]'
    '*Coupling[T3MEU11, {}, 0]^2'
    '*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]'
    '*Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]'
    '*NCM[Transp[Field[l, Fermion, {Index[d$$4, SU2L[fund]], Index[d$$2, Flavor]}, {}]], '
    'DiracProduct[GammaCC, Proj[-1]], '
    'Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}]])'
    '/Coupling[MF, {}, 0]'
)
E_EXTRA=A4.replace('+ (2*hbar','- (2*Sqrt[2/3]*hbar').replace(
    '*Coupling[lambdaT3, {}, 0]',
    '*Bar[Coupling[lambdaT3ExtraQ020101I1, {}, 0]]').replace(
    'Bar[Coupling[y2,','Bar[Coupling[y2Conj,').replace(
    '/(Sqrt[3]*Coupling[MF,','/Coupling[MF,')[:-1]


def test_source_pl_narrow_grammar():
    a=parse_native_kernel(A4)
    assert math.isclose(a.overall,2/math.sqrt(3))
    assert a.yukawas==(('y1',True),('y2',True))
    assert a.scalars==(('lambdaT3',False,1),)
    b=parse_native_kernel(B_MIX)
    assert b.overall==-1
    assert b.yukawas==(('y1',True),('y1',True))
    assert b.scalars==(('lambdaT3ExtraQ020200I1',True,1),
                       ('T3MEU11',True,2),('T3MEU11',False,2))
    e=parse_native_kernel(E_EXTRA)
    assert math.isclose(e.overall,-2*math.sqrt(2/3))
    assert e.yukawas[1][0]=='y2Conj'


def test_formal_pl_pole_cannot_be_parsed_as_finite():
    with pytest.raises(ValueError):parse_native_kernel(A4.replace('hbar*','hbar/(ϵ)*'))
    with pytest.raises(ValueError):parse_native_kernel(A4.replace('Proj[-1]','Proj[1]'))
    with pytest.raises(ValueError):parse_native_kernel(A4.replace('Coupling[lambdaT3, {}, 0]',
                                                            'Coupling[missing, {}, 0]'))


def test_matrix_inverse_mapping_not_zipped_aliases():
    # Source declares r=A m, not m=A r (with non-diagonal mixing).
    A=np.array([[0.,-2/math.sqrt(3)],[1.,1/math.sqrt(3)]])
    true=np.array([.02-.01j,.04+.005j])
    rhs=A@true
    entry={'matchete_symbols':['lambdaT3ExtraQ020020I1','lambdaT3ExtraQ020200I1'],
           'rgbeta_symbols':['lambdaNative1','lambdaNative2'],
           'map_rgbeta_from_matchete':{'real':A.tolist(),'imag':np.zeros_like(A).tolist()}}
    v,meta=matchete_quartics(entry,{'lambdaNative1':rhs[0],'lambdaNative2':rhs[1]})
    np.testing.assert_allclose([v[x] for x in entry['matchete_symbols']],true,atol=1e-14)
    assert meta['mapped_extra_count']==2


def test_native_heavy_sum_and_no_double_count():
    k=parse_native_kernel(A4)
    m=np.array([9.1e4,1e5,1.09e5])
    y1=np.diag([.02,.03,.04]).astype(complex)
    y2=np.diag([.01,.015,.02]).astype(complex)
    q=.015+.005j
    c,per,ledger=evaluate_kernels([k],m,{'y1':y1,'y2':y2},{'lambdaT3':q},None)
    expected=np.array([LOOP*(2/math.sqrt(3))*q*y1[i,i].conjugate()*y2[i,i].conjugate()
            *(1+math.log(MF**2/m[i]**2))/m[i] for i in range(3)])
    np.testing.assert_allclose(np.diag(c),expected,rtol=1e-13,atol=1e-28)
    assert np.linalg.norm(c-per.sum(axis=(0,1)))<1e-28
    assert len(ledger)==1


def test_mixed_U_conjugation_and_hbar_once():
    k=parse_native_kernel(B_MIX)
    U=np.eye(2,dtype=complex)
    U[0,0]=1j
    U[1,1]=-1j
    y=np.eye(3,dtype=complex)*.05
    m=np.full(3,1e5)
    q=.01+.02j
    c,per,_=evaluate_kernels([k],m,{'y1':y},{'lambdaT3ExtraQ020200I1':q},U)
    expected=-LOOP*np.conjugate(q)*(np.conjugate(1j)**2)*(1j**2)*(.05**2)/1e5
    np.testing.assert_allclose(c,np.eye(3)*expected,rtol=1e-13,atol=1e-27)
    assert np.linalg.norm(per)>0


def test_unknown_mixing_and_quartic_refused():
    k=parse_native_kernel(B_MIX)
    with pytest.raises(ValueError):
        evaluate_kernels([k],np.ones(3)*MF,{'y1':np.eye(3)},
                         {'lambdaT3ExtraQ020200I1':.01},None)
    with pytest.raises(ValueError):
        evaluate_kernels([k],np.ones(3)*MF,{'y1':np.eye(3)}, {},np.eye(2))


def test_sqrt_number():
    assert math.isclose(parse_scalar_number('2*Sqrt[2/3]'),2*math.sqrt(2/3))
    with pytest.raises(ValueError):parse_scalar_number('Sin[3]')
