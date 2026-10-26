"""Test exact Stage 22 diagnostics; production certification is never asserted."""
import copy
import json
from pathlib import Path

import pytest
import sympy as sp

from scripts.Stage22RGBetaAlignedGeometry import (
    CALIBRATION, CG, CG_SCALE, FUNDAMENTAL_REINDEX,
    _e_native_pattern, analyze, reindexed_doublet_tensor, triplet_orthonormality,
    triple_cg, triplet_tensor,
)
from scripts.Stage22TensorAlgebra import exact


DIR = Path('output/general_t3_stage22')


def _actual():
    paths=(DIR/'matchete_tensors.json',DIR/'rgbeta_native_invariants.json',DIR/'rgbeta_native_declarations.json')
    if any(not p.exists() for p in paths):
        pytest.skip('Stage22 native evidence unavailable in local repository; integration skipped')
    return [json.loads(p.read_text(encoding='utf-8-sig')) for p in paths]


def test_triplet_basis_is_orthonormal_and_symmetric():
    assert triplet_orthonormality()
    assert all(CG[a][0,1]==CG[a][1,0] for a in range(3))
    assert sp.simplify(CG_SCALE-sp.sqrt(2)/3**sp.Rational(1,4))==0


def test_triplet_cubic_is_antisymmetric_normalized():
    tensor=triple_cg()
    assert len(tensor)==27
    assert len([v for v in tensor.values() if v!=0])==6
    assert tensor[0,1,2]==sp.I/sp.sqrt(2)
    assert tensor[0,2,1]==-sp.I/sp.sqrt(2)


def test_e_native_structure_is_not_accepted_if_tampered():
    assert _e_native_pattern().count('eps[')==3
    assert _e_native_pattern().count('delS2[')==2
    assert _e_native_pattern().count(')/2')==2


def test_exact_matchete_calibration_against_upload():
    mat, inv, _ = _actual()
    assert CALIBRATION=='A_+0/021010'
    for model, code in [('A_+0','021010'),('A_-2','020101'),('D_+0','021010'),('D_-2','020101')]:
        m=next(t for t in mat['models'][model]['tensors'] if ''.join(map(str,t['counts']))==code)
        r=next(t for t in inv['models'][model]['invariants'] if 'delS2[' in t['native_invariant_expression_inputform'])
        expanded,_=triplet_tensor(r,m)
        assert all(sp.simplify(exact(a)-CG_SCALE*exact(b))==0
                   for a,b in zip(m['components'],expanded['components']))


def test_triple_triplet_e_model_prediction_has_exact_zero_residual():
    mat, inv, _ = _actual()
    for model,code in [('E_+0','021010'),('E_-2','020101')]:
        m=next(t for t in mat['models'][model]['tensors'] if ''.join(map(str,t['counts']))==code)
        r=next(t for t in inv['models'][model]['invariants'] if 'delS2[' in t['native_invariant_expression_inputform'])
        expanded,_=triplet_tensor(r,m)
        assert all(sp.simplify(exact(a)-sp.sqrt(2)*exact(b))==0
                   for a,b in zip(m['components'],expanded['components']))


def test_full_evidence_conditional_status_and_doublet_geometry():
    inputs=_actual()
    out=analyze(*inputs)
    assert out['totals']=={'DirectGeometryCandidate':32,'PseudorealReindexCandidate':8,
                            'CalibratedTripletCGCandidate':6,'IncompatibleNativePair':4,
                            'slots':50,'blocks':46}
    assert out['physical_coupling_map_certified'] is False
    assert out['native_rgbeta_component_basis_verified'] is False
    for model in ('B_-1','C_-1'):
        rotated=[b for b in out['models'][model] if b['status']=='PseudorealReindexCandidate']
        assert len(rotated)==4
        assert all(b['geometric']['r_in_m']==[['1']] and
                   b['component_rotation_applied']==list(FUNDAMENTAL_REINDEX)
                   for b in rotated)


def test_110101_blocks_exhaustive_incompatible_rotation_witness():
    report=analyze(*_actual())
    for model in ('B_-1','C_-1'):
        b=next(b for b in report['models'][model] if b['block']=='110101')
        assert b['status']=='IncompatibleNativePair'
        w=b['exhaustive_witness']
        assert w['rotation_choices_exhausted']==16
        assert w['consistent_I_or_J_alignment_exists'] is False
        assert w['common_full_rank_rotations']==[]
        assert all(len(v)>0 for v in w['native_coupling_rotation_witnesses'].values())


def test_tampered_e_triplet_expression_fails_closed():
    mat, inv, _=_actual()
    m=next(t for t in mat['models']['E_+0']['tensors'] if t['counts']==[0,2,1,0,1,0])
    t=copy.deepcopy(next(t for t in inv['models']['E_+0']['invariants'] if 'delS2[' in t['native_invariant_expression_inputform']))
    t['native_invariant_expression_inputform']=t['native_invariant_expression_inputform'].replace(')/2',')/3',1)
    with pytest.raises(ValueError,match='Unknown native triplet'):
        triplet_tensor(t,m)
