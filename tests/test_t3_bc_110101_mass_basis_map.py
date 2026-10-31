"""Regression tests for the native two-invariant B/C source geometry audit.

End-to-end checks deliberately consume committed native Stage22 and B/C mass-
basis evidence. No Wolfram or physical one-loop matching is performed.
"""
from pathlib import Path
import copy
import numpy as np
import pytest
import sympy as sp

from scripts.AuditT3BC110101MassBasisMap import (
    MAT, RGB, NATIVE, RGB_NATIVE, RGB_LIVE, MASS_BASIS, MODELS,
    source_tensor, check_live, exact_map, original_matchete_exact,
    compare_rotations, run_one, load,
)
from scripts.AuditT3BCNativeQuarticMassBasis import load_basis


def inputs():
    return load(NATIVE), load(RGB_NATIVE), load(RGB_LIVE)


def test_two_native_models_exact_same_rank_two_map():
    mat,inv,live=inputs()
    expected=sp.Matrix([[sp.Rational(1,2), 1],[-sp.sqrt(3)/2,0]])
    for model in MODELS:
        m=original_matchete_exact(mat,model)
        r=check_live(live,inv,model)
        a=exact_map(mat,model,r)
        assert a.shape==(2,2)
        assert all(sp.simplify(x)==0 for x in a-expected)
        assert sp.simplify(a.det()-sp.sqrt(3)/2)==0


def test_mutated_live_component_is_rejected():
    _mat,inv,live=inputs()
    fake=copy.deepcopy(live)
    fake['models']['B_-1']['native_covariance_audits'][0]['component_tensor_flat_inputform'][2]='0'
    with pytest.raises(ValueError,match='live RGBeta components'):
        check_live(fake,inv,'B_-1')


def test_mutated_coupling_factor_is_rejected():
    _mat,inv,live=inputs()
    fake=copy.deepcopy(inv)
    records=fake['models']['C_-1']['invariants']
    for r in records:
        if r['symbol']==RGB[0]:
            r['native_coupling_function_inputform']='2*Matrix[bogus][##1] &'
    with pytest.raises(ValueError,match='changed native expression/potential metadata'):
        check_live(live,fake,'C_-1')


def test_rotated_native_geometry_for_both_models():
    mat,inv,live=inputs()
    for model in MODELS:
        native=original_matchete_exact(mat,model)
        rgb=check_live(live,inv,model)
        a=exact_map(mat,model,rgb)
        basis=load_basis(MASS_BASIS/f'{model.replace("-","m")}_final_mass_basis.json',model)
        report=compare_rotations(native,rgb,a,basis.U,trials=5)
        assert report['exact_mass_basis_tensor_map_verified']
        assert report['mass_basis_rotated_component_max_relative_error'] < 1e-12
        assert report['rgbeta_original_vs_rotated_polynomial_max_relative_error'] < 1e-12


def test_results_do_not_claim_physical_certification():
    mat,inv,live=inputs()
    row=run_one('B_-1',mat_source=mat,inv=inv,live=live)
    assert row['status']=='Native110101TensorMapMassBasisChecked'
    assert row['physical_coupling_basis_map_certified'] is False
    assert row['verified_native_rgbeta_coupling_prefactors'] is False
    assert row['verified_native_rgbeta_matchete_lagrangian_sign'] is False
    assert row['full_mass_eigenstate_Matchete_matching_performed'] is False
    assert row['sequential_EFT_matching_performed'] is False


def test_changed_native_matchete_component_rejected_exactly():
    mat,inv,live=inputs()
    fake=copy.deepcopy(mat)
    tensor=next(r for r in fake['models']['B_-1']['tensors']
                if r['symbol']==MAT[0])
    index=next(i for i,x in enumerate(tensor['components']) if x=='1')
    tensor['components'][index]='2'
    rgb=check_live(live,inv,'B_-1')
    with pytest.raises(ValueError,match='not exactly in Matchete tensor span'):
        exact_map(fake,'B_-1',rgb)
