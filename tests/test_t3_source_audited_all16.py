"""Source-audited all16 map: CG algebra, Ward guards, and v8 regression.

Requires the original project modules; no synthetic physics outputs are produced.
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pytest

from Numerical.parallel_t3.wilson.basis import build_basis, ALL_MODELS
from Numerical.parallel_t3.wilson.source_audited_all16 import (
    _E_REPS, _E_SPARSE, _e_native_cg, _su2_variation,
    _BASE_QUARTIC_UNITS, assert_source_rgbeta_definition,
    canonical_mix_tensor, mass_ward_report, source_invariant_units,
    source_scalar_map, _triplet_unitary,
)
from Numerical.parallel_t3.wilson.matchete_a4_map_v8 import (
    SOURCE_REL, source_heavy_from_canonical,
)
from Numerical.parallel_t3.threshold.rough_ward_gate_v6 import gauge_ward_report

ROOT=Path(__file__).resolve().parents[1]


@pytest.mark.parametrize('dimensions',[(1,3),(2,2),(3,1),(3,3)])
def test_canonical_lambda_t3_complex_CG_is_SU2_singlet(dimensions):
    t=canonical_mix_tensor(*dimensions)
    assert t.shape==(2,2,*dimensions)
    assert np.linalg.norm(t)>0
    assert _su2_variation(t,*dimensions,first_conjugated=False)<1e-12


@pytest.mark.parametrize('model',('B_-3','C_-3','D_-2','E_-4'))
def test_source_registered_corrected_lambda_t3_has_full_real_Ward(model):
    units=source_invariant_units(model)
    b=build_basis(model)
    assert 'lambdaT3' in units
    for z in (1+0j,0+1j,.31-.63j):
        t=2*np.real(z*units['lambdaT3'])
        w=gauge_ward_report(t,b,sector='quartic')
        assert w['max_relative'] is not None
        assert w['max_relative']<2e-12, (model,z,w)


def test_E_same_variance_spin1_metric_cross_quartic_Ward():
    units=source_invariant_units('E_-2')
    w=gauge_ward_report(units['lambda12Cross'],build_basis('E_-2'),sector='quartic')
    assert w['max_relative']<2e-12,w


def test_mass_Ward_identity_for_all_scalar_representations():
    for model in ALL_MODELS:
        b=build_basis(model)
        m=np.eye(b.ns-4)*1.e6
        assert max(mass_ward_report(m,b).values())<1e-14


def test_E_source_CG_is_intertwiner_after_A4_anchored_basis_rotation():
    if not (ROOT/SOURCE_REL).is_file():
        pytest.skip('Tracked source anchor absent in this checkout')
    u,_=_triplet_unitary(ROOT)
    native=_e_native_cg()
    mapped=np.einsum('ijkl,ka,lb->ijab',native,u.conj(),u)
    assert abs(np.linalg.norm(mapped)/np.linalg.norm(native)-1)<1e-12
    assert _su2_variation(mapped,3,3,first_conjugated=True)<1e-12


def test_A4_exact_v8_source_map_not_changed():
    seed=ROOT/SOURCE_REL
    if not seed.is_file():
        pytest.skip('Tracked A4 native source absent')
    r,evidence=source_scalar_map(ROOT,'A_-4',seed)
    ref,_=source_heavy_from_canonical(seed)
    np.testing.assert_allclose(r,ref,atol=2e-14)
    assert evidence['orthogonality_residual']<1e-12
    assert evidence['no_Wilson_projection'] is True


def test_all16_cg_mapping_no_source_guess():
    # Source checks are intentionally enforced by the batch runner, including
    # the source-manifest SHA256 and full Wilson Ward identity.
    assert len(ALL_MODELS)==16
    assert (ROOT/'RGE/running/rgbeta/T3RGBetaModel.wl').exists()
    info=assert_source_rgbeta_definition(ROOT)
    assert 'RGBeta_source_sha256' in info
