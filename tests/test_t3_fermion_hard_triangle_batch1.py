"""Checks for Batch 1 F-hard source resolution and C5 normalization.

The actual 16-model source reproduction and end-to-end numerical calculation
are performed by scripts/RunT3FHardTriangleBatch1.py, not by these tests.
"""
from __future__ import annotations

from types import SimpleNamespace
import numpy as np

from Numerical.parallel_t3.wilson.fermion_hard_triangle_batch1 import (
    _basis_transform, _uv_at_mf, convert_hard_projector_normalisation, MF,
)
from Numerical.parallel_t3.threshold.matching import master_weinberg_template, project_master
from Numerical.lo_certification.completeness.hard_matching import weinberg_tensor


def test_native_hard_to_master_Weinberg_convention():
    hard=weinberg_tensor()
    master=master_weinberg_template()
    np.testing.assert_allclose(hard,2*master,rtol=2e-15,atol=2e-15)
    sample=np.array([[0.1+0.02j,0.003+0.001j,0.02j],
                     [0.003+0.001j,0.4,0.05],
                     [0.02j,0.05,0.9-0.2j]])*1e-12
    full=np.einsum('pq,ijAB->pqijAB',sample,hard,optimize=True)
    c5,err=project_master(full)
    np.testing.assert_allclose(c5,2*sample,rtol=1e-13,atol=1e-27)
    assert err<1e-13
    np.testing.assert_allclose(convert_hard_projector_normalisation(sample),2*sample)


def test_heavy_rotation_is_full_tensor_not_a_projection():
    rng=np.random.default_rng(25)
    k,n=3,8
    raw=rng.normal(size=(k,3,3,2,2,n,n))+1j*rng.normal(size=(k,3,3,2,2,n,n))
    X=rng.normal(size=(n,n))
    R,_=np.linalg.qr(X)
    transformed=_basis_transform(raw,R,1.0)
    restored=np.einsum('kfgijAB,aA,bB->kfgijab',transformed,R,R,optimize=True)
    np.testing.assert_allclose(restored,raw,rtol=2e-14,atol=2e-14)


def test_native_uv_mf_endpoint_without_changing_trajectory():
    raw={'mu_uv_gev':np.array([1e6,1e5]),
         'uv__MF':np.stack([np.eye(3)*90e3,np.eye(3)*MF]).astype(complex),
         'uv__y1':np.ones((2,3,3),complex)*.012,
         'uv__y2':np.ones((2,3,3),complex)*.022,
         'uv__y2Conj':np.ones((2,3,3),complex)*.007}
    bg=SimpleNamespace(arrays=raw)
    uv=_uv_at_mf(bg,'E_-2')
    assert uv.mu_gev==MF and not uv.representation.majorana_fermion
    assert set(uv.general_couplings)=={'y2Conj'}
    np.testing.assert_array_equal(uv.MF,np.eye(3)*MF)
    np.testing.assert_array_equal(raw['uv__y1'][0],np.ones((3,3),complex)*.012)
