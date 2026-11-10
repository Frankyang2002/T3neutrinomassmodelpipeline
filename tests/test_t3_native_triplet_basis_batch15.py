"""Matchete CG-driven SU2 source geometry: independent weak Ward regression."""
from __future__ import annotations
import importlib.util
from pathlib import Path
import math
import numpy as np
import pytest

P=Path(__file__).resolve().parents[1]/'Numerical/parallel_t3/complete/native_triplet_basis_batch15.py'
spec=importlib.util.spec_from_file_location('t3_cg15',P)
b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)


def native_cg():
    # Exact six components from source CGRegistry T3MixCG; source SHA is
    # verified separately in the on-repository integration test.
    x=np.zeros((2,2,3),complex);a=3**(-.25)
    x[0,0,0]=a;x[0,0,2]=1j*a;x[0,1,1]=-a
    x[1,0,1]=-a;x[1,1,0]=-a;x[1,1,2]=1j*a
    return x


def real_generator(t):
    d=t.shape[0];out=np.zeros((2*d,2*d),complex)
    for p in range(d):
        for q in range(d):
            a,c=t[p,q].real,t[p,q].imag
            out[2*p:2*p+2,2*q:2*q+2]=1j*np.array([[c,a],[-a,c]])
    return out


def source_weak_polynomial(model='A_-4'):
    x=native_cg();C=np.zeros((2,2,8,8),complex)
    for i in range(2):
        for j in range(2):
            for trip in range(3):
                for ising in range(2):
                    for itrip in range(2):
                        v=x[i,j,trip]*(1,-1j)[ising]*(1,1j)[itrip]/2
                        if model=='A_-4':
                            a,c=ising,2+2*trip+itrip
                        else:
                            a,c=2*trip+itrip,6+ising
                        C[i,j,a,c]+=v;C[i,j,c,a]+=v
    return (C+C.swapaxes(0,1)+C.swapaxes(2,3)+C.transpose(1,0,3,2))/4


def test_matchete_cartesian_generator_is_derived_not_assumed():
    gs=b.generators_from_cg(native_cg())
    assert np.allclose(gs[2],np.array([[0,0,-1j],[0,0,0],[1j,0,0]]))
    assert all(np.allclose(g,g.conj().T) for g in gs)
    assert np.allclose(gs[0]@gs[1]-gs[1]@gs[0],1j*gs[2])
    assert not np.allclose(gs[2],b.spherical_spin_one()[2])


def test_exact_source_CG_change_to_spherical_basis():
    U,meta=b.source_from_spherical(native_cg())
    assert meta['cg_triplet_spherical_intertwining_residual_max']<1e-12
    assert np.allclose(U.conj().T@U,np.eye(3),atol=1e-12)
    for source,ref in zip(b.generators_from_cg(native_cg()),b.spherical_spin_one()):
        assert np.allclose(U.conj().T@source@U,ref,atol=1e-12)
    R=b.complex_to_interleaved_real(U)
    assert np.allclose(R.T@R,np.eye(6),atol=1e-12)


@pytest.mark.parametrize('model,offset',[('A_-4',2),('D_-2',0)])
def test_weak_ward_failures_reproduced_exactly_then_closed(model,offset):
    C=source_weak_polynomial(model)
    spherical=b.spherical_spin_one()
    cart=b.generators_from_cg(native_cg())
    wrong=[];right=[]
    for t,ts,tc in zip(b.fundamental(),spherical,cart):
        old=np.zeros((8,8),complex);old[offset:offset+6,offset:offset+6]=real_generator(ts)
        src=np.zeros((8,8),complex);src[offset:offset+6,offset:offset+6]=real_generator(tc)
        wrong.append(b.relative_operator_ward(C,t,old))
        right.append(b.relative_operator_ward(C,t,src))
    assert wrong==pytest.approx([math.sqrt(2+math.sqrt(2)),math.sqrt(2+math.sqrt(2)),2],abs=1e-12)
    assert max(right)<1e-12


@pytest.mark.parametrize('model',('A_-4','D_-2'))
def test_source_to_reference_orthogonal_heavy_rotations(model):
    # Test local representation map independently of the repository parser.
    U,_=b.source_from_spherical(native_cg())
    R=b.complex_to_interleaved_real(U)
    out=np.eye(8);offset=2 if model=='A_-4' else 0
    out[offset:offset+6,offset:offset+6]=R
    assert np.linalg.norm(out.T@out-np.eye(8))<1e-12
    assert np.allclose(out[2:8,0:2],0) if model=='A_-4' else np.allclose(out[6:8,:6],0)


def test_malformed_or_rank_deficient_CG_blocks_instead_of_approximating():
    with pytest.raises(ValueError,match='triplet'):
        b.generators_from_cg(np.zeros((2,2,3)))
    with pytest.raises(ValueError,match='2x2x3'):
        b.generators_from_cg(np.zeros((2,2,2)))


def test_optional_project_archive_CG_agrees_with_exact_six_values():
    root=Path(__file__).resolve().parents[1]
    f=(root/'output/general_observable/fermion_first_all16_v1/A_m4/data/eft1_after_F_wilson_seed.json')
    if not f.is_file():
        pytest.skip('Actual Matchete CG archived only in local project checkout')
    import json
    tree=json.loads(f.read_text(encoding='utf-8-sig'))
    cg=b.triplet_cg_from_seed(tree)
    assert np.allclose(cg,native_cg(),rtol=1e-11,atol=1e-14)
