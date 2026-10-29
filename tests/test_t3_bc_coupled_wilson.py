"""Numerical Wilson master and strict interface regression tests."""
from __future__ import annotations
from dataclasses import dataclass
from types import SimpleNamespace
import numpy as np
import pytest
import sympy as sp

from Numerical.experimental.T3BCCoupledWilsonTensor import (
 PAIRS, flavor_pairs, independent_tensor, rhs_master, interpolate_kernels,
 integrate_coupled, flavor_to_matrix, _symmetric_lookup, _sympy_expression_lambdas,
 weinberg_projection, matchete_from_rgbeta,
)


def symseed():
    a,b,m=sp.symbols('y1 y2 MF')
    items={(1,1,5,5):sp.conjugate(a)**2/m,
           (2,2,7,7):sp.conjugate(b)**2/m,
           (1,2,5,7):2*sp.conjugate(a)*sp.conjugate(b)/m}
    for i,j,a,b in list(items):
        for ii,jj,aa,bb in ((j,i,a,b),(i,j,b,a),(j,i,b,a)):
            items[(ii,jj,aa,bb)]=items[(i,j,a,b)]
    return SimpleNamespace(nonzero_items=lambda:items.items())


def test_physical_flavor_matrices_symmetric():
    y1=np.arange(9).reshape(3,3)*(.03+.01j)
    y2=np.arange(9).reshape(3,3)*(.02-.005j)
    f=flavor_pairs(np.array([1500.,2000.,2700.]),y1,y2)
    assert f.shape==(3,6)
    assert np.allclose(flavor_to_matrix(f[1]),flavor_to_matrix(f[1]).T)


@pytest.mark.parametrize('bad',[[1,2],[1,-2,3],[1,2,float('nan')]])
def test_invalid_masses_rejected(bad):
    with pytest.raises(ValueError):flavor_pairs(np.asarray(bad),np.ones((3,3)),np.ones((3,3)))


def test_mixed_majorana_flavor_lift_uses_ordered_sum_not_average():
    """Regression: the previous cross-sector average halved the saved direct beta."""
    m=np.array([1600.0,2350.0,3100.0])
    y1=np.array([[.11+.04j,.02-.03j,.06+.01j],
                 [.03+.02j,.09+.01j,-.01+.03j],
                 [.02-.01j,.06+.03j,.07+.02j]])
    y2=np.array([[.07-.02j,.05+.02j,.03-.01j],
                 [-.01+.04j,.04+.05j,.08-.02j],
                 [.02+.02j,.01-.03j,.09+.01j]])
    f11,f12,f22=flavor_pairs(m,y1,y2)
    a=np.conjugate(y1);b=np.conjugate(y2)
    # This is exactly the two-orientation formula in the existing production
    # WeinbergFlavorMatching/RunT3DirectWilsonSequential full-flavor lift.
    expected_cross=(a/m[None,:])@b.T+(b/m[None,:])@a.T
    np.testing.assert_allclose(flavor_to_matrix(f12),expected_cross,rtol=1e-14,atol=1e-20)
    np.testing.assert_allclose(flavor_to_matrix(f11),(a/m[None,:])@a.T,rtol=1e-14,atol=1e-20)
    np.testing.assert_allclose(flavor_to_matrix(f22),(b/m[None,:])@b.T,rtol=1e-14,atol=1e-20)
    assert np.linalg.norm(f12)>0
    # Avoid a vacuous symmetrical-only check that would also pass an average.
    assert np.linalg.norm(flavor_to_matrix(f12)-expected_cross/2)>0.4*np.linalg.norm(expected_cross)


def test_cross_sector_tree_lift_matches_majorana_two_orientation_boundary():
    """The tree mixed symbol is lifted to F12, not F12/2."""
    m=np.array([1400.0,1900.0,2300.0])
    y1=np.array([[.01,.02,.03],[.04,.05,.06],[.07,.08,.09]],dtype=complex)
    y2=np.array([[.03,.02,.01],[.06,.05,.04],[.09,.08,.07]],dtype=complex)
    factors=flavor_pairs(m,y1,y2)
    y1s,y2s=sp.symbols('y1 y2')
    mf=sp.Symbol('MF',nonzero=True)
    mixed=sp.conjugate(y1s)*sp.conjugate(y2s)/mf
    terms={(1,2,5,7):mixed,(2,1,5,7):mixed,
           (1,2,7,5):mixed,(2,1,7,5):mixed}
    seed=SimpleNamespace(nonzero_items=lambda:terms.items())
    c=independent_tensor(seed,factors)
    expected=(np.conjugate(y1)/m[None,:])@np.conjugate(y2).T+ (np.conjugate(y2)/m[None,:])@np.conjugate(y1).T
    np.testing.assert_allclose(flavor_to_matrix(c[0,1,4,6,:]),expected,rtol=1e-14,atol=1e-20)


def test_seed_lift_and_mixed_symmetrization():
    f=np.array([[1,2,3,4,5,6],[2,3,4,5,6,7],[4,5,6,7,8,9]],dtype=complex)
    C=independent_tensor(symseed(),f)
    assert C.shape==(3,3,12,12,6)
    assert np.allclose(C[0,0,4,4,:],f[0])
    assert np.allclose(C[1,1,6,6,:],f[2])
    assert np.allclose(C[0,1,4,6,:],2*f[1])
    assert np.allclose(C,C.swapaxes(0,1))


def test_seed_lift_handles_matchete_nonzero_mass_symbol():
    # Production ScalarOnlyTensorAdapters._term_prefactor constructs
    # sp.Symbol("MF", nonzero=True), not the unassumed Symbol("MF").
    y1,y2=sp.symbols('y1 y2')
    mf=sp.Symbol('MF', nonzero=True)
    terms={(1,1,5,5):sp.I*sp.conjugate(y1)**2/(8*mf),
           (1,2,5,7):sp.conjugate(y1)*sp.conjugate(y2)/(3*mf)}
    for (i,j,a,b),x in list(terms.items()):
        for key in ((j,i,a,b),(i,j,b,a),(j,i,b,a)):
            terms[key]=x
    seed=SimpleNamespace(nonzero_items=lambda:terms.items())
    factors=np.asarray([[2,3,4,5,6,7],[8,9,10,11,12,13],[14,15,16,17,18,19]],dtype=complex)
    result=independent_tensor(seed,factors)
    np.testing.assert_allclose(result[0,0,4,4,:],1j*factors[0]/8)
    np.testing.assert_allclose(result[0,1,4,6,:],factors[1]/3)
    assert np.isfinite(result).all()


def test_seed_lift_rejects_ambiguous_mass_assumptions():
    y=sp.Symbol('y1')
    unassumed=sp.Symbol('MF')
    nonzero=sp.Symbol('MF',nonzero=True)
    expr=y*y*(1/unassumed+1/nonzero)
    seed=SimpleNamespace(nonzero_items=lambda:[((1,1,5,5),expr)])
    with pytest.raises(ValueError,match='ambiguous symbolic assumptions'):
        independent_tensor(seed,np.ones((3,6)))


def test_seed_lift_rejects_missing_symbols():
    seed=SimpleNamespace(nonzero_items=lambda:[((1,1,1,1),sp.Symbol('unknown')+1)])
    with pytest.raises(ValueError,match='unresolved'):independent_tensor(seed,np.ones((3,6)))


def fake_kernels(quartic=0.,gauge=0.,yukawa=0.,gamma=0.):
    q=np.zeros((12,12,12,12),dtype=complex)
    for a in range(12):q[a,a,a,a]=quartic
    sg=np.zeros_like(q)
    sg[0,0,0,0]=gauge
    y=np.zeros((3,3,12),dtype=complex)
    y[0,2,0]=yukawa;y[2,0,0]=yukawa
    gs=np.eye(12,dtype=complex)*gamma
    gf=np.eye(3,dtype=complex)*gamma
    return {'y':y,'quartic':q,'scalarg':sg,
            'mixed':np.zeros((3,3,12,12),dtype=complex),
            'gammas':gs,'gammaf':gf}


def test_zero_rhs_exact():
    c=np.ones((3,3,12,12,6),dtype=complex)*1e-8
    assert np.max(abs(rhs_master(c,fake_kernels())))==0


def test_scalar_pair_index_and_normalization():
    c=np.zeros((3,3,12,12,6),dtype=complex)
    c[0,0,0,0,0]=2.
    r=rhs_master(c,fake_kernels(quartic=3.))
    assert r[0,0,0,0,0]==6.
    rg=rhs_master(c,fake_kernels(gauge=4.))
    assert rg[0,0,0,0,0]==-8.


def test_gamma_scalar_and_fermion_permutations():
    c=np.zeros((3,3,12,12,6),dtype=complex)
    c[0,0,0,0,0]=2.
    r=rhs_master(c,fake_kernels(gamma=3.))
    assert r[0,0,0,0,0]==24.  # 2 scalar + 2 fermion gamma factors


def test_conjugate_coefficient_term_is_present():
    c=np.zeros((3,3,12,12,6),dtype=complex)
    c[2,2,0,0,0]=1+2j
    r=rhs_master(c,fake_kernels(yukawa=.2))
    # Nonzero generation into LL/eC channels from C* and Yukawa mixing.
    assert np.max(abs(r))>0


def test_source_basis_stays_symmetric_under_master():
    rng=np.random.default_rng(3)
    c=rng.normal(size=(3,3,12,12,2))+1j*rng.normal(size=(3,3,12,12,2))
    c=(c+c.swapaxes(0,1)+c.swapaxes(2,3)+c.swapaxes(0,1).swapaxes(2,3))/4
    r=rhs_master(c,fake_kernels(yukawa=.2,gamma=.3))
    assert np.allclose(r,r.swapaxes(0,1))
    assert np.allclose(r,r.swapaxes(2,3))


def test_quartic_lookup_is_fully_permutation_symmetric():
    a=_symmetric_lookup({(1,1,2,3):3+0j})
    assert a[0,0,1,2]==3 and a[2,0,1,0]==3
    assert np.count_nonzero(a)==12


def test_quartic_symbol_guard():
    with pytest.raises(ValueError,match='unmapped'):
        _sympy_expression_lambdas({(1,1,1,1):sp.Symbol('unknownNewCoupling')})


def test_analytic_projection_full_components():
    c=np.zeros((3,3,12,12,6),dtype=complex)
    template={(1,1,3,3):.5,(1,1,3,4):.25j,(1,1,4,4):-.5,
              (2,2,1,1):.5,(2,2,1,2):.25j,(2,2,2,2):-.5,
              (1,2,1,3):-.25,(1,2,1,4):-.25j,(1,2,2,3):-.25j}
    vals=np.arange(1,7,dtype=float)*(.005+.003j)
    for (i,j,a,b),factor in template.items():
        for ii,jj in ((i,j),(j,i)):
            for aa,bb in ((a,b),(b,a)):
                c[ii-1,jj-1,aa-1,bb-1,:]=factor*vals
    projected,res=weinberg_projection(c,template,np.array([0,1,2]))
    assert res<1e-12
    assert np.allclose(projected,vals)


def test_projection_refuses_non_weinberg_support():
    template={(1,1,3,3):.5,(1,1,3,4):.25j,(1,1,4,4):-.5,
              (2,2,1,1):.5,(2,2,1,2):.25j,(2,2,2,2):-.5,
              (1,2,1,3):-.25,(1,2,1,4):-.25j,(1,2,2,3):-.25j}
    c=np.zeros((3,3,12,12,6),dtype=complex)
    c[0,0,2,2,:]=1
    c[0,1,0,2,:]=2
    c[0,1,0,3,:]=3
    with pytest.raises(ValueError,match='not a pure Weinberg'):
        weinberg_projection(c,template,np.array([0,1,2]))


def test_interpolation_and_grid_range():
    t=np.exp(np.linspace(np.log(1e5),np.log(1e3),17))
    k=[fake_kernels(quartic=float(i)) for i in range(len(t))]
    middle=interpolate_kernels(np.log(t),k,(np.log(t[0])+np.log(t[1]))/2)
    assert np.isclose(middle['quartic'][0,0,0,0],.5)
    with pytest.raises(ValueError,match='outside EFT1'):
        interpolate_kernels(np.log(t),k,np.log(t[0])+1)


def test_closed_system_zero_coupling_identity():
    mu=np.geomspace(1e5,1e3,17)
    initial=np.zeros((3,3,12,12,1),dtype=complex);initial[0,0,0,0,0]=1e-11
    trajectory,sol=integrate_coupled(initial,mu,[fake_kernels() for _ in mu])
    assert sol.success
    assert np.allclose(trajectory[...,-1],initial,rtol=1e-9,atol=1e-20)


def test_stage29_matrix_inverse_is_not_assumed_diagonal():
    h={'format':'stage29-conditional-uv-quartic-handoff-v1',
       'physical_coupling_mapping_certified':False,
       'models':{'B_-1':{'blocks':[{'rgbeta_symbols':['lambdaT3GeneralQ110101I1','lambdaHHdagS1barS2barCross'],
               'matchete_symbols':['mat1','mat2'],
               'signed_map_exact':[['0','2*sqrt(3)/3'],['-1','-sqrt(3)/3']]}]}}}
    state=SimpleNamespace(extra={'lambdaT3GeneralQ110101I1':2*np.sqrt(3)/3},
                          baseline=SimpleNamespace(lambdaHHdagS1barS2barCross=-1-np.sqrt(3)/3))
    mapped=matchete_from_rgbeta('B_-1',h,state)
    assert np.isclose(mapped['mat1'],1) and np.isclose(mapped['mat2'],1)


def test_stage29_requires_unverified_status():
    h={'format':'stage29-conditional-uv-quartic-handoff-v1',
       'physical_coupling_mapping_certified':True,'models':{}}
    with pytest.raises(ValueError,match='certification changed'):
        matchete_from_rgbeta('B_-1',h,None)


def test_context_source_paths_if_repo_present():
    from pathlib import Path
    root=Path(__file__).resolve().parents[1]
    f=root/'output/general_t3_stage22/matchete_tensors.json'
    if not f.is_file():pytest.skip('full public repository not available in isolated container')
    import json
    from Numerical.experimental.T3BCCoupledWilsonTensor import extra_matchete_templates
    data=json.loads(f.read_text(encoding='utf-8'))
    for model in ('B_-1','C_-1'):
        templates=extra_matchete_templates(data,model)
        assert len(templates)==7


def brute_component(C,k,i,j,a,b,f):
    """Independent index-loop transcription of the project's seven equations."""
    y=k['y'];quartic=k['quartic'];sg=k['scalarg'];mixed=k['mixed']
    gs=k['gammas'];gf=k['gammaf'];n=3;s=12
    v=0j
    for kk in range(n):
        for ll in range(n):
            for d in range(s):
                v+=y[i,j,d]*np.conj(y[kk,ll,d])*C[kk,ll,a,b,f]
                v+=(2*y[i,ll,d]*y[j,kk,d]+y[i,j,d]*y[kk,ll,d])*np.conj(C[kk,ll,a,b,f])
    for aa,bb in ((a,b),(b,a)):
        for c in range(s):
            for d in range(s):
                v+=.5*(quartic[aa,bb,c,d]-sg[aa,bb,c,d])*C[i,j,c,d,f]
    for aa,bb in ((a,b),(b,a)):
        for ii,jj in ((i,j),(j,i)):
            for ll in range(n):
                for c in range(s):
                    v+=mixed[jj,ll,bb,c]*C[ii,ll,aa,c,f]
    for aa,bb,ii,jj in ((a,b,i,j),(b,a,j,i)):
        for kk in range(n):
            for ll in range(n):
                for c in range(s):
                    v+=y[jj,ll,aa]*y[ii,kk,c]*C[ll,kk,bb,c,f]
    for aa,bb in ((a,b),(b,a)):
        for d in range(s):v+=gs[bb,d]*C[i,j,aa,d,f]
    for ii,jj in ((i,j),(j,i)):
        for kk in range(n):v+=gf[jj,kk]*C[ii,kk,a,b,f]
    return v


@pytest.mark.parametrize('a,b,i,j',[(0,0,0,0),(0,1,0,1),(1,2,1,2)])
def test_seven_term_numeric_contractions_against_explicit_index_loops(a,b,i,j):
    rng=np.random.default_rng(34)
    C=(rng.normal(size=(3,3,12,12,1))+1j*rng.normal(size=(3,3,12,12,1)))*1e-5
    C=(C+C.swapaxes(0,1)+C.swapaxes(2,3)+C.swapaxes(0,1).swapaxes(2,3))/4
    y=(rng.normal(size=(3,3,12))+1j*rng.normal(size=(3,3,12)))*.05
    y=(y+y.swapaxes(0,1))/2
    s=(rng.normal(size=(12,12,12,12))*.001)
    s=(s+s.swapaxes(0,1)+s.swapaxes(2,3)+s.swapaxes(0,1).swapaxes(2,3))/4
    gamma_s=np.diag(rng.normal(size=12)*.1)
    gamma_f=np.diag(rng.normal(size=3)*.1)
    mixed=(rng.normal(size=(3,3,12,12))+1j*rng.normal(size=(3,3,12,12)))*.02
    k={'y':y,'quartic':s,'scalarg':s*.3,'mixed':mixed,
       'gammas':gamma_s,'gammaf':gamma_f}
    expected=brute_component(C,k,i,j,a,b,0)
    got=rhs_master(C,k)[i,j,a,b,0]
    assert np.isclose(expected,got,rtol=1e-11,atol=2e-16)
