"""Regression checks for the read-only MF/Beta source Ward diagnostic.

These tests do NOT certify native Matchete/RGBeta normalizations or supply model
predictions. Invariant gauge tensors should pass; nonsinglets should be flagged.
"""
from __future__ import annotations
import json
import numpy as np
import pytest

from Numerical.parallel_t3.threshold.rough_ward_gate_v6 import (
    ward_numerator, gauge_ward_report, read_background, check_portal_units,
    _origin_flags, GATE_TOLERANCE)
from Numerical.parallel_t3.wilson.basis import build_basis
from Numerical.parallel_t3.wilson.coordinates import embed_c5


def test_holomorphic_yukawa_u1_charge_cancellation():
    # Two left Weyl fields each -1/4, and H of +1/2.
    tf=np.array([[-.25]],complex)
    ts=np.array([[0,.5j],[-.5j,0]],complex)
    y=np.zeros((1,1,2),complex)
    y[0,0]=[1,1j]
    np.testing.assert_allclose(ward_numerator(y,tf,ts),0,atol=1e-15)
    y[0,0]=[1,0]
    assert np.linalg.norm(ward_numerator(y,tf,ts))>0.1


def test_o2_symmetric_quartic_is_invariant_and_anisotropy_is_detected():
    ts=np.array([[0,1j],[-1j,0]],complex)
    d=np.eye(2)
    q=(np.einsum('ab,cd->abcd',d,d)+np.einsum('ac,bd->abcd',d,d)
       +np.einsum('ad,bc->abcd',d,d))
    np.testing.assert_allclose(ward_numerator(q,None,ts),0,atol=1e-15)
    wrong=q.copy();wrong[0,0,0,0]+=1
    assert np.linalg.norm(ward_numerator(wrong,None,ts))>0


def test_full_weinberg_source_passes_project_gauge_ward():
    basis=build_basis('A_-4')
    c5=np.diag([1+1j,2-3j,3+0.3j])*1e-12
    full=embed_c5(c5,basis)
    quality=gauge_ward_report(full,basis,sector='wilson')
    assert quality['passed'] is True
    assert quality['max_relative']<1e-12


def test_wilson_hypercharge_invalid_component_detected():
    basis=build_basis('A_-4')
    c=np.zeros((basis.nf,basis.nf,basis.ns,basis.ns),complex)
    # LL H H† neutral-hypercharge operator is not Weinberg.
    i=basis.leptons[0][0]
    c[i,i,0,0]=1
    report=gauge_ward_report(c,basis,sector='wilson')
    assert report['passed'] is False
    assert report['ward_relative']['U1']>0.1


def test_flags_are_per_sector_and_per_group():
    ok={'ward_relative':{'SU2':0.,'U1':0.,'SU3':0.}}
    bad={'ward_relative':{'SU2':1.1,'U1':0.,'SU3':0.}}
    found=_origin_flags(ok,ok,bad,ok)
    assert found==['MF_full_quartic_SU2_WardViolation']
    assert _origin_flags(ok,ok,ok,ok)==[]


def test_known_conditional_source_quartic_units_checked_not_silently_replaced():
    a=check_portal_units('A_-4')
    b=check_portal_units('B_-3')
    assert a['lambdaT3_present'] and b['lambdaT3_present']
    # Source-specific CG consistency test, not a claim that any unit's
    # normalization matches the Matchete source.
    assert a['quartic_component_units']['lambdaT3']['real']['ward_relative']['SU2']<1e-9
    assert isinstance(b['noninvariant_unit_names'],list)


def test_wrong_background_checksum_blocks_before_wilson_analysis(tmp_path):
    basis=build_basis('A_-4')
    n=2
    npz=tmp_path/'background.npz'
    np.savez_compressed(npz,mu_gev=np.array([1e5,1e3]),
        Y=np.zeros((n,basis.nf,basis.nf,basis.ns),complex),
        lambda4=np.zeros((n,)+(basis.ns,)*4),
        gY=np.array([.35,.36]),g2=np.array([.65,.64]),g3=np.array([1.0,1.2]))
    meta=tmp_path/'background.json'
    meta.write_text(json.dumps({'model':'A_-4','npz_sha256':'f'*64,
                   'fermion_labels':list(basis.fermion_labels),
                   'scalar_labels':list(basis.scalar_labels)}))
    with pytest.raises(ValueError, match='checksum'):
        read_background(npz,meta,basis)
