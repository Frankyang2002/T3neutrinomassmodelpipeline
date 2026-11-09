"""Source-CG derived A(-4) triplet map: no Ward projection or coefficient loss."""
from __future__ import annotations

import json
from pathlib import Path
import numpy as np
import pytest

from Numerical.parallel_t3.wilson.basis import build_basis
from Numerical.parallel_t3.wilson.coordinates import embed_c5,embed_formal_llss
from Numerical.parallel_t3.wilson.matchete_a4_map_v8 import (
    SOURCE_REL,SOURCE_SPARSE,SOURCE_REPS,canonical_singlet_cg,
    source_cg,derive_triplet_unitary,realify_complex_unitary,
    source_heavy_from_canonical,convert_a4_boundary)
from Numerical.parallel_t3.wilson.ward_singlet_v7 import _gauge_norm


def source_document():
    return {'IntegratedFields':['F'],'ActiveHeavyFields':['S1','S2'],
            'TreeWilsonTerms':[
                {'Chirality':'PR','TermInputForm':'CG[Bar[T3MixCG],...]'},
                {'Chirality':'PL','TermInputForm':(
                    'Bar[Field[NewScalar1, Scalar, {}, {}]] * '
                    'CG[T3MixCG, {Bar[Index[a]], Bar[Index[b]], Index[c]}] * '
                    'Field[NewScalar2, Scalar, {Index[c]}, {}]')}],
            'CGRegistry':[{'Name':'T3MixCG','RepsInputForm':SOURCE_REPS,
                           'TensorInputForm':SOURCE_SPARSE}]}


def seed(tmp_path):
    path=tmp_path/'eft1_after_F_wilson_seed.json'
    path.write_text(json.dumps(source_document()),encoding='utf-8')
    return path


def source_raw(g):
    """Single-flavour Matchete PL (L L) S1* S2 in native real S2 coordinates."""
    raw=np.zeros((3,3,2,2,8,8),complex)
    singletbar=np.array([1,-1j])/np.sqrt(2)
    triplet=np.array([1,1j])/np.sqrt(2)
    for i in range(2):
        for j in range(2):
            for k in range(3):
                for a in range(2):
                    for b in range(2):
                        v=g[i,j,k]*singletbar[a]*triplet[b]/2
                        raw[0,0,i,j,a,2+2*k+b]+=v
                        raw[0,0,i,j,2+2*k+b,a]+=v
    return raw


def test_source_sparse_exact_and_hash(tmp_path):
    path=seed(tmp_path)
    g,evidence=source_cg(path)
    assert g.shape==(2,2,3)
    assert g[0,0,2]==pytest.approx(1j*3**(-.25))
    assert evidence['source_sha256']
    assert evidence['source_chirality_used']=='PL'


def test_source_sparse_rejects_changed_cg(tmp_path):
    path=seed(tmp_path)
    doc=source_document()
    doc['CGRegistry'][0]['TensorInputForm']=SOURCE_SPARSE.replace('-3^(-1/4)','+3^(-1/4)',1)
    path.write_text(json.dumps(doc))
    with pytest.raises(ValueError,match='T3MixCG changed'):
        source_cg(path)


def test_source_cg_requires_PL_operator(tmp_path):
    path=seed(tmp_path);doc=source_document();doc['TreeWilsonTerms'][1]['Chirality']='PR'
    path.write_text(json.dumps(doc))
    with pytest.raises(ValueError,match='PL LL'):
        source_cg(path)


def test_analytic_triplet_unitary_and_source_cg_reconstruction(tmp_path):
    g,_=source_cg(seed(tmp_path))
    u,scale=derive_triplet_unitary(g)
    expected=np.array([[-1/np.sqrt(2),0,1/np.sqrt(2)],
                       [0,1,0],[-1j/np.sqrt(2),0,-1j/np.sqrt(2)]])
    np.testing.assert_allclose(u,expected,atol=1e-14)
    np.testing.assert_allclose(g.reshape(4,3)@u,
                               scale*canonical_singlet_cg().reshape(4,3),atol=1e-13)
    np.testing.assert_allclose(u.conj().T@u,np.eye(3),atol=1e-13)
    assert scale==pytest.approx(np.sqrt(6)*3**(-.25))


def test_realification_interleaved_is_orthogonal_and_complex_linear(tmp_path):
    g,_=source_cg(seed(tmp_path));u,_=derive_triplet_unitary(g)
    r=realify_complex_unitary(u)
    np.testing.assert_allclose(r.T@r,np.eye(6),atol=1e-13)
    z=np.array([1+2j,-3+4j,2-1j]);v=np.stack((z.real,z.imag),axis=1).ravel()
    out=(r@v).reshape(3,2)
    np.testing.assert_allclose(out[:,0]+1j*out[:,1],u@z,atol=1e-13)


def test_source_exact_map_preserves_all_norm_and_restores_SU2(tmp_path):
    path=seed(tmp_path);g,_=source_cg(path)
    raw=source_raw(g);basis=build_basis('A_-4')
    baseline=embed_formal_llss(raw,basis,np.eye(8),1.)
    r,source=source_heavy_from_canonical(path)
    mapped,report=convert_a4_boundary(baseline,raw,1.,r)
    assert _gauge_norm(baseline,basis)['SU2']==pytest.approx(1.5086889589691435,abs=1e-12)
    assert _gauge_norm(mapped,basis)['SU2']<1e-12
    assert report['LLSS_norm_fraction_retained']==pytest.approx(1.,abs=1e-13)
    assert source['no_coefficient_projection'] is True
    np.testing.assert_allclose(np.linalg.norm(mapped),np.linalg.norm(baseline),atol=1e-20,rtol=1e-13)


def test_hard_weinberg_sector_unchanged(tmp_path):
    path=seed(tmp_path);g,_=source_cg(path);b=build_basis('A_-4');raw=source_raw(g)
    c5=np.diag([1+2j,2,-3j])*1e-10
    original=embed_formal_llss(raw,b,np.eye(8),1)+embed_c5(c5,b)
    r,_=source_heavy_from_canonical(path)
    mapped,report=convert_a4_boundary(original,raw,1.,r)
    np.testing.assert_array_equal(original[:,:,:4,:4],mapped[:,:,:4,:4])
    assert report['outside_LLSS_unchanged'] is True
    assert _gauge_norm(mapped,b)['SU2']<1e-12


def test_wrong_triplet_orientation_fails_closed(tmp_path):
    path=seed(tmp_path);g,_=source_cg(path);b=build_basis('A_-4');raw=source_raw(g)
    before=embed_formal_llss(raw,b,np.eye(8),1)
    r,_=source_heavy_from_canonical(path)
    bad=r.copy();bad[2:,2:]=r[2:,2:].T
    with pytest.raises(ValueError,match='does not restore'):
        convert_a4_boundary(before,raw,1.,bad)


def test_projected_boundary_cannot_masquerade_as_source(tmp_path):
    path=seed(tmp_path);g,_=source_cg(path);b=build_basis('A_-4');raw=source_raw(g)
    original=embed_formal_llss(raw,b,np.eye(8),1)
    r,_=source_heavy_from_canonical(path)
    mapped,_=convert_a4_boundary(original,raw,1.,r)
    with pytest.raises(ValueError,match='does not match saved Matchete formal LLSS'):
        convert_a4_boundary(mapped,raw,1.,r)


def test_current_repo_source_file_if_present():
    root=Path(__file__).resolve().parents[1]
    path=root/SOURCE_REL
    if not path.is_file():
        pytest.skip('A(-4) original Matchete source not in this checkout')
    g,evidence=source_cg(path)
    u,_=derive_triplet_unitary(g)
    np.testing.assert_allclose(u.conj().T@u,np.eye(3),atol=1e-12)
    assert evidence['source_sha256']
