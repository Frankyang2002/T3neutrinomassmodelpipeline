"""Regression for source-derived B/C 110101 quartic potential reconstruction.

Uses saved native evidence and scalar-basis reports; does not run Wolfram or
certify the installed RGBeta runtime's Lam[] or Matchete Match[] semantics.
"""
from __future__ import annotations
import copy
import numpy as np
import pytest
import sympy as sp

from scripts.AuditT3BC110101MassBasisMap import (
    NATIVE, RGB_NATIVE, RGB_LIVE, MASS_BASIS,
    load, check_live, exact_map, original_matchete_exact,
)
from scripts.AuditT3BCNativeQuarticMassBasis import load_basis
from scripts.AuditT3BC110101PotentialMap import (
    MODELS, RGB, NORMALIZATION, EXPECTED_A, MAT_SOURCE,
    metadata_guard, normalization_identity, source_guard,
    check_potential_ensemble, compute_potentials, run_one,
)


def source_data():
    return load(NATIVE),load(RGB_NATIVE),load(RGB_LIVE)


def component_data(model):
    mat, inv, live=source_data()
    rgb=check_live(live,inv,model)
    A=exact_map(mat,model,rgb)
    m=original_matchete_exact(mat,model)
    basis=load_basis(MASS_BASIS/f'{model.replace("-","m")}_final_mass_basis.json',model)
    return m,rgb,A,basis


def test_exact_24_to_6_conversion_and_native_permutations():
    assert normalization_identity()==1
    assert sp.simplify(sp.Integer(6)/sp.factorial(4)*(sp.sqrt(2)**4)-1)==0
    _, inv, _=source_data()
    for model in MODELS:
        meta=metadata_guard(inv,model)
        assert meta['unique_arrangements']==24
        assert meta['rgbeta_internal_coupling_factor']==6
        assert meta['hermitian_conjugate_enabled']


def test_matchete_plus_hc_source_guard():
    source_guard(MAT_SOURCE.read_text(encoding='utf-8-sig'))
    with pytest.raises(ValueError,match='Matchete extra quartic construction changed'):
        source_guard('contribution = coupling[] product If[cg === None, 1, CG[cg, labels]];')


def test_both_models_potential_ensemble_and_mass_basis():
    for model in MODELS:
        mat,rgb,A,basis=component_data(model)
        assert all(sp.simplify(x)==0 for x in A-EXPECTED_A)
        report=check_potential_ensemble(mat,rgb,A,basis.U,trials=20)
        assert report['max_rgbeta_real_to_complex_potential_relative_error']<1e-12
        assert report['max_rgbeta_potential_plus_matchete_lagrangian_relative_error']<1e-12
        assert report['max_mass_basis_potential_relative_error']<1e-12


def test_independent_complex_phases_and_hc():
    mat,rgb,A,basis=component_data('B_-1')
    vec=lambda a,b: np.array([a,b],dtype=complex)
    data=compute_potentials(mat,rgb,A,basis.U,vec(0.21+0.47j,-0.38+0.18j),
                            vec(0.9+0.7j,-0.2+0.4j),
                            vec(0.5-0.3j,0.4+0.1j),
                            vec(-0.6+0.9j,0.7+0.3j))
    assert abs(data['rgbeta_v_from_real_lam']+data['matchete_l_from_native_hc'])<1e-12
    assert abs(data['rgbeta_v_from_real_lam']-data['rgbeta_v_in_mass_basis'])<1e-12
    assert abs(data['complex_rgbeta_polynomial']+data['complex_matchete_polynomial'])<1e-12
    assert abs(data['complex_rgbeta_polynomial'].imag)>1e-5  # h.c. is nontrivial


def test_wrong_native_factor_fails_closed():
    _,inv,_=source_data()
    fake=copy.deepcopy(inv)
    row=next(r for r in fake['models']['C_-1']['invariants']
             if r['symbol']==RGB[1])
    row['native_coupling_function_inputform']=row['native_coupling_function_inputform'].replace('6*','12*')
    with pytest.raises(ValueError,match='factor is not 6'):
        metadata_guard(fake,'C_-1')


def test_missing_arrangement_or_conjugate_fails_closed():
    _,inv,_=source_data()
    bad=copy.deepcopy(inv)
    row=next(r for r in bad['models']['B_-1']['invariants'] if r['symbol']==RGB[0])
    row['native_unique_arrangements_inputform']='{{1, 2, 3, 4}}'
    with pytest.raises(ValueError,match='24 native arrangements'):
        metadata_guard(bad,'B_-1')
    row['native_unique_arrangements_inputform']=next(
        r['native_unique_arrangements_inputform']
        for r in inv['models']['B_-1']['invariants'] if r['symbol']==RGB[0])
    row['self_conjugate']=True
    with pytest.raises(ValueError,match='Hermitian conjugate'):
        metadata_guard(bad,'B_-1')


def test_wrong_sign_fails_potential_comparison():
    mat,rgb,A,basis=component_data('C_-1')
    with pytest.raises(ValueError,match='potential conversion failed'):
        check_potential_ensemble(mat,rgb,-A,basis.U,trials=10)


def test_runner_keeps_uncertified_flags():
    mat,inv,live=source_data()
    for model in MODELS:
        row=run_one(model,mat,inv,live,trials=8)
        assert row['coupling_map_exact_h1']=='-g1/2 - g2'
        assert row['coupling_map_exact_h2']=='sqrt(3)*g1/2'
        assert row['status']=='SourceReconstructedQuarticPotentialMapChecked'
        assert row['physical_coupling_mapping_certified'] is False
        assert row['installed_rgbeta_lam_runtime_reconstruction_performed'] is False
        assert row['installed_matchete_lagrangian_runtime_evaluated'] is False
        assert row['one_loop_matching_performed'] is False
