"""Tests for the A(-4) *diagnostic* invariant-subspace boundary projection."""
from __future__ import annotations
import json
import numpy as np
import pytest

from Numerical.parallel_t3.wilson.basis import build_basis
from Numerical.parallel_t3.wilson.coordinates import embed_c5
from Numerical.parallel_t3.wilson.synthetic_invariants import synthetic_llss
from Numerical.parallel_t3.wilson.ward_singlet_v7 import (
    STATUS,project_a4_boundary,_gauge_norm)


def _source():
    b=build_basis('A_-4')
    mat=np.array([[1+2j,.2-.8j,-.1], [.2-.8j,2+.3j,.45j],[-.1,.45j,-1j]])*1e-10
    return b,synthetic_llss(b,mat)


def test_invariant_LLSS_is_fixed_by_projection():
    b,c=_source()
    projected,report=project_a4_boundary(c)
    np.testing.assert_allclose(projected,c,atol=1e-23,rtol=1e-12)
    assert report['status']==STATUS
    assert report['kept_LLSS_frobenius_fraction']==pytest.approx(1.,abs=1e-12)
    assert report['discarded_LLSS_frobenius_fraction']<1e-12
    assert max(report['projected_Ward'].values())<1e-12


def test_c5_and_every_other_non_LLSS_sector_are_untouched():
    b,c=_source()
    c+=embed_c5(np.eye(3)*2e-12,b)
    f=b.leptons[0][0]
    # A very large explicitly non-singlet LLSS component.
    c[f,f,4,4]+=1e-10
    projected,report=project_a4_boundary(c)
    np.testing.assert_array_equal(projected[:,:,0:4,0:4],c[:,:,0:4,0:4])
    np.testing.assert_array_equal(projected[:,:,:4,4:],c[:,:,:4,4:])
    assert report['non_LLSS_tensor_unchanged']
    assert report['discarded_LLSS_frobenius_fraction']>.01
    assert report['original_Ward']['SU2']>1e-3
    assert max(report['projected_Ward'].values())<1e-12
    assert np.linalg.norm(projected-c)>0


def test_orthogonal_and_idempotent():
    b,c=_source();f=b.leptons[0][0]
    c[f,f,4,4]+=2e-10
    out,diag=project_a4_boundary(c)
    twice,again=project_a4_boundary(out)
    np.testing.assert_allclose(twice,out,atol=1e-23,rtol=1e-12)
    assert diag['orthogonal_projection_residual']<1e-11
    assert again['discarded_LLSS_frobenius_fraction']<1e-12


def test_reject_zero_boundary():
    b,_=_source()
    with pytest.raises(ValueError,match='Zero LLSS'):
        project_a4_boundary(np.zeros((b.nf,b.nf,b.ns,b.ns),complex))


def test_reject_unsupported_scalar_family():
    _,c=_source()
    with pytest.raises(ValueError,match='only for A_-4'):
        project_a4_boundary(c,model='D_-2')


def test_reject_bad_geometry_and_pair_symmetry():
    b,c=_source()
    with pytest.raises(ValueError,match='geometry'):
        project_a4_boundary(c[:10])
    c[b.leptons[0][0],b.leptons[1][0],4,4]+=1e-8
    with pytest.raises(ValueError,match='exchange symmetry'):
        project_a4_boundary(c)


def test_reject_unrepaired_other_sectors():
    b,c=_source()
    # Wrong hypercharge in the HH block cannot be repaired by LLSS-only work.
    c[b.leptons[0][0],b.leptons[0][0],0,0]+=1e-10
    with pytest.raises(ValueError,match='still violates'):
        project_a4_boundary(c)


def test_json_report_and_no_conjectured_native_conversion():
    _,c=_source()
    _,report=project_a4_boundary(c)
    json.dumps(report,allow_nan=False)
    assert report['not_physical_matching'] is True
    assert report['unverified_Matchete_RGBeta_native_scalar_basis_map'] is True
    assert report['all_fermions_integrated_at_MF'] is True
    assert report['both_scalars_remain_active_to_MS'] is True
    assert len(report['projection_components'])==6
