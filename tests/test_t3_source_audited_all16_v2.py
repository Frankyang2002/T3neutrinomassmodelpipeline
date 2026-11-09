"""Regression for the six blocked v1 source-basis cases; no ODE required.

These tests read the original tracked Workstream B tensors, not synthetic
physics output. Running the v2 batch remains the end-to-end numerical test.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from Numerical.parallel_t3.background.contracts import tag
from Numerical.parallel_t3.wilson.basis import build_basis
from Numerical.parallel_t3.wilson.coordinates import embed_formal_llss
from Numerical.parallel_t3.wilson.ward_singlet_v7 import _gauge_norm
from Numerical.parallel_t3.wilson.source_audited_all16 import (
    _self_named_from_interleaved, covariance_corrected_mass,
    source_scalar_map,
)
from Numerical.parallel_t3.threshold.models import (
    SELF_S1, SELF_S2, _mixed_rotation, from_parameters,
)

ROOT=Path(__file__).resolve().parents[1]
BLOCKED=('A_-2','B_-1','C_-1','D_+0','E_-2','E_+0')


def _saved(model:str):
    base=ROOT/'output/parallel_t3/background/B_all16_qualitative_v1'/tag(model)
    manifest=json.loads((base/'fermion_boundary_manifest.json').read_text(encoding='utf-8-sig'))
    native=ROOT/manifest['source_matchete_seed'].replace('\\','/')
    assert native.is_file()
    return base,manifest,native


@pytest.mark.parametrize('model',BLOCKED)
def test_blocked_model_source_map_is_orthogonal_and_exact(model):
    _,meta,seed=_saved(model)
    r,report=source_scalar_map(ROOT,model,seed)
    b=build_basis(model)
    np.testing.assert_allclose(r.T@r,np.eye(b.ns-4),atol=2e-12)
    assert report['no_Wilson_projection'] is True
    if model in ('B_-1','C_-1'):
        diag=meta['mixed_native_diagnostics']['mixing_U_muS']
        u=np.asarray(diag['real'])+1j*np.asarray(diag['imag'])
        np.testing.assert_allclose(r,_mixed_rotation(u).T,atol=2e-12)
    else:
        assert report['source_self_bilinear_field_packing']=='Minus_all_then_Plus_all'
        assert report['source_CG_orientation_verified']


@pytest.mark.parametrize('model',BLOCKED)
def test_saved_formal_wilson_becomes_full_gauge_singlet(model):
    base,_,seed=_saved(model)
    r,_=source_scalar_map(ROOT,model,seed)
    with np.load(base/'fermion_boundary_formal.npz',allow_pickle=False) as data:
        formal=np.asarray(data['C_LLSS_formal'],complex)
    b=build_basis(model)
    c=embed_formal_llss(formal,b,r,1.)
    assert np.linalg.norm(c)>0
    residual=_gauge_norm(c,b)
    assert max(residual.values())<1e-7,(model,residual)
    # The mapping is orthogonal, with no deleted Wilson tensor components.
    before=embed_formal_llss(formal,b,np.eye(len(r)),1.)
    assert np.linalg.norm(c)==pytest.approx(np.linalg.norm(before),rel=1e-12)


@pytest.mark.parametrize('model',('A_-2','D_+0','E_-2','E_+0'))
def test_self_triplet_scalar_mass_is_covariant_and_isospectral(model):
    old=from_parameters(model,1.8e6,2.1e6,mu_sq=1.8e4+1.4e4j)
    mass,report=covariance_corrected_mass(ROOT,model,old.mass2_canonical_gev2)
    assert report['corrected'] is True
    assert max(report['corrected_mass_Ward'].values())<1e-10
    np.testing.assert_allclose(np.linalg.eigvalsh(mass),
                               np.linalg.eigvalsh(old.mass2_canonical_gev2),
                               atol=1e-6)


def test_named_triplet_source_packing_exact():
    P=_self_named_from_interleaved(3)
    canonical=np.array([11,12,21,22,31,32])
    np.testing.assert_array_equal(P@canonical,[11,21,31,12,22,32])
    np.testing.assert_array_equal(P@P.T,np.eye(6))
