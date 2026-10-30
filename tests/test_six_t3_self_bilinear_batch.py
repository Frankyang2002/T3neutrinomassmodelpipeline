"""Standalone algebra/regression checks for six-self-bilinear batch.
Full Matchete/RGBeta and source-CG validation occurs in live batch, not here.
"""
import numpy as np
import pytest
import mpmath as mp
from scripts.VerifySixT3SelfBilinearBatch import (
    MODELS, MASS_NAMES, EXTRA_NAMES, triangle_integral, candidate_channels,
    lift_flavour, validate_native_result, tag,
)


def test_six_models_and_expected_mass_symbols():
    assert set(MODELS)=={'A_-2','A_+0','D_-2','D_+0','E_-2','E_+0'}
    assert MASS_NAMES == {-2:'muT3GeneralS2Sq',0:'muT3GeneralS1Sq'}
    assert EXTRA_NAMES[-2][1]=='lambdaT3GeneralQ020101I1'
    assert EXTRA_NAMES[0][1]=='lambdaT3GeneralQ021010I1'
    assert tag('A_-2')=='A_m2'


@pytest.mark.parametrize('x,y,z', [(1e10,1.306122448979592e6,1.02e6),
                                    (1.21e10,1.14e6,9.8e5),
                                    (1.44e10,2e6,1.2e6)])
def test_triangle_matches_independent_log_identity(x,y,z):
    actual=triangle_integral(x**.5,y,z,dps=55)
    with mp.workdps(55):
        masses=[mp.mpf(str(t)) for t in (x,y,z)]
        identity=sum(q*mp.log(q)/mp.fprod(q-v for v in masses if v!=q)
                     for q in masses)
        assert abs(actual-identity)/max(abs(actual),abs(identity))<mp.mpf('1e-13')


@pytest.mark.parametrize('alpha',[-2,0])
def test_complex_phases_and_relative_sign(alpha):
    p=np.diag([2,3,4]).astype(complex)
    m=np.diag([3,4,6]).astype(complex)
    mu=1+2j
    extra=.02+.01j
    got=candidate_channels(p,m,lam_t3=.4+.2j,extra_rgb=extra,mu_rgb=mu,alpha=alpha)
    phase=mu/abs(mu) if alpha==-2 else mu.conjugate()/abs(mu)
    expected=-extra.conjugate()*phase*.5*(p-m)
    assert np.allclose(got['anomalous'],expected)
    assert np.allclose(got['normal'],(.4+.2j)*.5*(p+m))
    assert np.allclose(got['total'],got['normal']+got['anomalous'])
    opposite=candidate_channels(p,m,lam_t3=.4+.2j,extra_rgb=extra,mu_rgb=mu,
                                alpha=alpha,matchete_sign=+1)
    assert np.allclose(got['anomalous'],-opposite['anomalous'])


def test_zero_coupling_and_zero_mu():
    p=np.eye(3,dtype=complex)*2
    m=np.eye(3,dtype=complex)*4
    for a in (0,-2):
        assert np.linalg.norm(candidate_channels(p,m,lam_t3=1,extra_rgb=0,
                                  mu_rgb=2,alpha=a)['anomalous'])==0
        assert np.linalg.norm(candidate_channels(p,m,lam_t3=1,extra_rgb=1,
                                  mu_rgb=0,alpha=a)['anomalous'])==0
    with pytest.raises(ValueError,match='invalid spectral sign'):
        candidate_channels(p,m,lam_t3=1,extra_rgb=1,mu_rgb=1,alpha=2)


def test_flavour_lift_complex_symmetric():
    y1=np.array([[1,2j,3],[-.2,1,3j],[1j,2,1.5]],dtype=complex)
    y2=np.array([[.5,2,1j],[1,2j,.3],[.4j,1,-.2]],dtype=complex)
    weights=np.array([.1,.2,.3])
    got=lift_flavour(weights,y1,y2)
    ref=sum(weights[r]*(np.outer(y1[:,r].conj(),y2[:,r].conj())+
                        np.outer(y2[:,r].conj(),y1[:,r].conj())) for r in range(3))
    assert np.allclose(got,ref)
    assert np.allclose(got,got.T)


@pytest.mark.parametrize('model,dim,species',[
    ('A_-2',3,'S2dag'),('A_+0',1,'S1dag'),
    ('D_-2',1,'S2dag'),('D_+0',3,'S1dag'),
    ('E_-2',3,'S2dag'),('E_+0',3,'S1dag')])
def test_native_tensor_certificate_schema(model,dim,species):
    o={'format':'t3-six-self-native-bilinear-tensor-v1',
       'status':'NativeMassTensorExtracted','model':model,'dimension':dim,
       'species':[species,species],'native_metric_proportional_identity':True,
       'mass_operator_in_full_lagrangian':True,'operator_names':['T3ExtraM_X'],
       'rgbeta_matchete_vertex_map_certified':False,
       'native_metric_factor_inputform':'1',
       'native_metric_inputform':[['1' if i==j else '0' for j in range(dim)]
                                  for i in range(dim)]}
    assert validate_native_result(o,model)==1
    o['native_metric_inputform'][0][0]='2'
    with pytest.raises(ValueError,match='not proportional'):
        validate_native_result(o,model)


def test_reject_nondiagonal_false_symmetry():
    bad=np.array([[1,2,0],[0,1,0],[0,0,1]],dtype=complex)
    with pytest.raises(ValueError,match='bad symmetric C5 channels'):
        candidate_channels(bad,np.eye(3),lam_t3=1,extra_rgb=1,mu_rgb=1,alpha=0)

@pytest.mark.parametrize('model,dim',[
    ('A_-2',3),('A_+0',1),('D_-2',1),('D_+0',3),('E_-2',3),('E_+0',3)])
def test_independent_complex_phase_mass_hessian_and_propagators(model,dim):
    # Independent real/imaginary Hessian with the declared spherical metric.
    m=1e6
    mu=1e4*np.exp(.35j)
    C=(np.eye(1) if dim==1 else np.array([[0,0,1],[0,-1,0],[1,0,0]]))
    B=2*mu.conjugate()*C
    hessian=np.block([[m*np.eye(dim)+B.real,-B.imag],
                      [-B.imag,m*np.eye(dim)-B.real]])
    eigen=np.linalg.eigvalsh(hessian)
    want=np.sort(np.repeat([m-2*abs(mu),m+2*abs(mu)],dim))
    assert np.allclose(eigen,want,atol=1e-7)
    p2=8e5
    G=np.linalg.inv(p2*np.eye(2*dim)+hessian)
    xx,yy=G[:dim,:dim],G[dim:,dim:]
    xy,yx=G[:dim,dim:],G[dim:,:dim]
    normal=.5*(xx+yy+1j*(yx-xy))
    anomalous=.5*(xx-yy+1j*(yx+xy))
    gplus,gminus=[1/(p2+v) for v in (m+2*abs(mu),m-2*abs(mu))]
    assert np.allclose(normal,.5*(gplus+gminus)*np.eye(dim),rtol=3e-13)
    assert np.allclose(anomalous,.5*(gplus-gminus)*mu/abs(mu)*C,rtol=3e-13,atol=1e-18)
