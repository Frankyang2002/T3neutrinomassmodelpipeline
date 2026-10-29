"""Self-bilinear exact mass basis and topology-normal propagator tests."""
import numpy as np
import pytest

from Numerical.experimental.T3SelfBilinearMassBasis import (
    SelfBilinearMassBasis, invariant_metric, real_intertwiner,
    rotate_component_tensor, evaluate_tensor,
)


@pytest.mark.parametrize('dim',[1,3])
def test_integer_isospin_real_intertwiner(dim):
    C=invariant_metric(dim);W=real_intertwiner(dim)
    assert np.allclose(W.conj().T@W,np.eye(dim))
    assert np.allclose(W.T@C@W,np.eye(dim))
    assert np.allclose(C@W.conjugate(),W)


@pytest.mark.parametrize('dim',[1,3])
@pytest.mark.parametrize('mu',[500+600j,-700+20j,0j])
def test_potential_and_hessian_spectra(dim,mu):
    basis=SelfBilinearMassBasis.build('S1',dim,2.3e6,mu)
    rng=np.random.default_rng(199+dim)
    for _ in range(3):
        real=rng.normal(size=dim)
        imag=rng.normal(size=dim)
        assert np.isclose(basis.potential_old(basis.fields(real,imag)),
                          basis.potential_mass_basis(real,imag),rtol=2e-13)
    assert np.allclose(np.sort(basis.eigen_masses_squared),
                       [2.3e6-2*abs(mu),2.3e6+2*abs(mu)])


@pytest.mark.parametrize('dim',[1,3])
@pytest.mark.parametrize('mu',[complex(70,18),complex(-60,200),0j])
def test_full_propagator_against_real_hessian_inverse(dim,mu):
    basis=SelfBilinearMassBasis.build('S2',dim,500,mu)
    C=invariant_metric(dim)
    B=2*mu.conjugate()*C
    H=np.block([[500*np.eye(dim)+B.real,-B.imag],
                [-B.imag,500*np.eye(dim)-B.real]])
    p2=31
    G=np.linalg.inv(p2*np.eye(2*dim)+H)
    xx=G[:dim,:dim];yy=G[dim:,dim:]
    xy=G[:dim,dim:];yx=G[dim:,:dim]
    normal_exact=.5*(xx+yy+1j*(yx-xy))
    abnormal_exact=.5*(xx-yy+1j*(yx+xy))
    gn,ga=basis.propagators(p2)
    assert np.allclose(gn,normal_exact,rtol=2e-13,atol=2e-17)
    assert np.allclose(ga,abnormal_exact,rtol=2e-13,atol=2e-17)
    assert np.allclose(gn,np.eye(dim)*np.trace(gn)/dim)
    assert (np.linalg.norm(ga)>0)==(mu!=0)


@pytest.mark.parametrize('dim',[1,3])
def test_tensor_transport_for_supplied_arbitrary_vertex(dim):
    rng=np.random.default_rng(dim+17)
    basis=SelfBilinearMassBasis.build('S1',dim,2000,150+32j)
    fields=('S1','S1dag','H','S1')
    t=rng.normal(size=(dim,dim,2,dim))+1j*rng.normal(size=(dim,dim,2,dim))
    rotated=rotate_component_tensor(t,fields,basis)
    assert len(rotated)==8
    for _ in range(5):
        R=rng.normal(size=dim);I=rng.normal(size=dim)
        H=rng.normal(size=2)+1j*rng.normal(size=2)
        S=basis.fields(R,I)
        original=evaluate_tensor(t,(S,S.conjugate(),H,S))
        vals={'S1_R':R,'S1_I':I,'H':H}
        moved=sum(evaluate_tensor(v,[vals[name] for name in names])
                  for names,v in rotated.items())
        assert np.allclose(original,moved,rtol=3e-12,atol=1e-10)


def test_excludes_tachyonic_mass_and_unsupported_dimensions():
    assert SelfBilinearMassBasis.build('S1',1,10,6).tachyonic
    with pytest.raises(ValueError):real_intertwiner(2)
    with pytest.raises(ValueError):SelfBilinearMassBasis.build('S1',2,12,1)
    with pytest.raises(ValueError):SelfBilinearMassBasis.build('H',1,12,1)
    with pytest.raises(ValueError):SelfBilinearMassBasis.build('S1',1,12,complex(float('nan'),0))


def test_zero_self_bilinear_propagator_is_unsplit():
    basis=SelfBilinearMassBasis.build('S2',3,1200,0)
    g,a=basis.propagators(4)
    assert np.allclose(g,np.eye(3)/1204)
    assert np.allclose(a,0)


def test_vertex_transport_rejects_missing_selected_fields():
    basis=SelfBilinearMassBasis.build('S1',3,100,4)
    with pytest.raises(ValueError):rotate_component_tensor(np.ones((2,2)),('H','H'),basis)
    with pytest.raises(ValueError):rotate_component_tensor(np.ones((2,2)),('S1','H'),basis)


@pytest.mark.parametrize('text,expected',[
    ('Sqrt[3/11]',np.sqrt(3/11)),
    ('2*Sqrt[3/11]',2*np.sqrt(3/11)),
    ('3^(-1/4)',3**(-.25)),
    ('(-I)/3^(1/4)',-1j/3**.25),
    ('(-I)/Sqrt[2]',-1j/np.sqrt(2)),
    ('-(1/Sqrt[2])',-1/np.sqrt(2)),
])
def test_saved_stage22_matchete_numeric_inputform(text,expected):
    from Numerical.experimental.T3SelfBilinearMassBasis import parse_matchete_constant
    assert np.isclose(parse_matchete_constant(text),expected,rtol=3e-14)


def test_native_extra_tensor_catalogue_transport(tmp_path):
    import json
    from Numerical.experimental.T3SelfBilinearMassBasis import native_tensors_for_model
    # Synthetic A(+0) catalogue with the same field dimensions and number
    # of independent quartic slots as the saved Stage22 source.
    tensor={'symbol':'lambdaT3ExtraQ021010I1',
            'species':['Hdag','Hdag','S1','S2'],
            'dimensions':[2,2,1,3],
            'components':['1']*12,'potential_prefactor':'1'}
    rows=[dict(tensor,symbol=f'lambdaT3ExtraQ021010I{i}') for i in range(1,6)]
    obj={'format':'stage22-exact-tensors-v1',
         'source':'MatcheteInvariantTensors','model_count':16,
         'quartic_slot_count':50,'models':{'A_+0':{'tensors':rows}}}
    path=tmp_path/'stage22.json'
    path.write_text(json.dumps(obj),encoding='utf-8')
    basis=SelfBilinearMassBasis.build('S1',1,2500,40+25j)
    results=native_tensors_for_model(path,'A_+0',basis)
    assert len(results)==5
    assert set(results[0]['rotated'])=={
        ('Hdag','Hdag','S1_R','S2'),('Hdag','Hdag','S1_I','S2')}
    for item in results:
        for value in item['rotated'].values():
            assert value.shape==(2,2,1,3)
    obj['models']['A_+0']['tensors'][0]['components'][0]='unsafe_symbol'
    path.write_text(json.dumps(obj),encoding='utf-8')
    with pytest.raises(ValueError,match='Unsupported'):
        native_tensors_for_model(path,'A_+0',basis)
