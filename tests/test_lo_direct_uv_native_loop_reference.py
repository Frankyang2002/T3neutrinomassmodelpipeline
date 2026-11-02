"""Independent reference physics checks (synthetic, not physical T3 certification)."""
import numpy as np
import pytest

from Numerical.lo_certification.direct_uv.native_amplitude_contract import ContractError
from Numerical.lo_certification.direct_uv.native_uv_loop_reference import (
    b0_exact, cartesian_hessian, decode, encode, independent_reference,
    potential, ratio_function, rephase_input,
)


def ma_toy(y=0.2, majorana=100., a=42000., b=1000.):
    # Complex scalar mass splitting B and Majorana single Weyl field.
    return {'fermion_mass_gev':encode([[majorana]]),
            'scalar_A_gev2':encode([[a]]),
            'scalar_B_gev2':encode([[b]]),
            'renormalization_scale_gev':1.e7,
            'vertices':[{'name':'y_HC','field':encode(np.array([[[y]]]*3))}]}


def dirac_toy():
    M=np.array([[0.,150.],[150.,0.]],complex)
    A=np.array([[40000.,2000.],[2000.,52000.]],complex)
    y1=np.zeros((3,2,2),complex)
    y2=np.zeros((3,2,2),complex)
    y1[:,0,0]=[.1,.2,.3]
    y2[:,1,1]=[.2,.1,-.4]
    return {'fermion_mass_gev':encode(M),'scalar_A_gev2':encode(A),
            'scalar_B_gev2':encode(np.zeros((2,2))),
            'renormalization_scale_gev':1.e7,
            'vertices':[{'name':'y1_HC','conjugate':encode(y1)},
                        {'name':'y2_HC','field':encode(y2)}]}


def test_b0_exact_log_primitive_against_numeric_integration():
    from scipy.integrate import quad
    cases=[(1e10,1e6,1e7),(4e4,4e4*(1+1e-9),1000.),
           (4e4,4e4*(1+0.1),1000.),(1e6,1e10,2e5)]
    for a,b,mu in cases:
        numeric=quad(lambda t:-np.log((t*a+(1-t)*b)/(mu**2)),0,1,
                     epsabs=2e-11,epsrel=2e-11)[0]
        assert abs(b0_exact(a,b,mu)-numeric)<5e-10


def test_hessian_reconstructs_quadratic_potential_with_complex_mixing():
    A=np.array([[5.,1+1.3j],[1-1.3j,9.]],complex)
    B=np.array([[.1+.5j,.3-.2j],[.3-.2j,-.4+.2j]],complex)
    H=cartesian_hessian(A,B)
    x=np.array([.6,-.9,.3,.1])
    phi=(x[:2]+1j*x[2:])/np.sqrt(2)
    assert np.allclose(x@H@x/2,potential(A,B,phi),atol=1e-14)


def test_single_complex_scalar_unbroken_limit_exact_zero():
    d=ma_toy(b=0.)
    r=independent_reference(d)
    m=decode(r['complex_symmetric_matrix_gev'],'m',(3,3))
    assert np.linalg.norm(m)<1e-14
    assert r['diagnostics']['uv_pole_relative_to_natural_amplitude']<1e-13


def test_majorana_split_scalar_closed_form_no_spurious_factor_two():
    d=ma_toy(y=.13,b=700.)
    r=independent_reference(d)
    nu=decode(r['complex_symmetric_matrix_gev'],'m',(3,3))
    y=.13; mf=100.;a=42000.;b=700.;mu=1.e7
    # For nu psi y phi, phi=(R+iI)/sqrt2, relative I vertex is i:
    predicted=(1/(16*np.pi**2))*(y*y*mf/2)*(
        b0_exact(mf*mf,a+b,mu)-b0_exact(mf*mf,a-b,mu))
    assert np.allclose(nu, predicted*np.ones((3,3)),rtol=3e-10,atol=1e-16)
    assert abs(predicted)>0


def test_dirac_mixed_scalar_full_mass_nonzero_and_pole_free():
    r=independent_reference(dirac_toy())
    z=decode(r['complex_symmetric_matrix_gev'],'nu',(3,3))
    assert np.linalg.norm(z)>1e-9
    assert np.allclose(z,z.T)
    assert r['diagnostics']['uv_pole_relative_to_natural_amplitude']<1e-12
    # Same-source terms vanish by Dirac chirality and complex scalar flow.
    assert np.linalg.norm(decode(r['source_vertex_pair_matrices_gev']['y1_HC x y1_HC'],'y1'))<1e-14
    assert np.linalg.norm(decode(r['source_vertex_pair_matrices_gev']['y2_HC x y2_HC'],'y2'))<1e-14


def test_full_complex_field_rephasing_covariance():
    d=dirac_toy()
    baseline=independent_reference(d)
    p=np.exp(1j*np.array([.34,-.93]));s=np.exp(1j*np.array([.47,.27]))
    after=independent_reference(rephase_input(d,p,s))
    old=decode(baseline['complex_symmetric_matrix_gev'],'a')
    new=decode(after['complex_symmetric_matrix_gev'],'b')
    assert np.linalg.norm(old-new)/np.linalg.norm(old)<1e-8


def test_wrong_handedness_or_yukawa_orientation_is_detected():
    d=dirac_toy()
    wrong=dict(d)
    wrong['vertices']=[dict(v) for v in d['vertices']]
    v=wrong['vertices'][0]
    v['field']=v.pop('conjugate')
    baseline=independent_reference(d)
    altered=independent_reference(wrong)
    a=decode(baseline['complex_symmetric_matrix_gev'],'a')
    b=decode(altered['complex_symmetric_matrix_gev'],'b')
    assert np.linalg.norm(a-b)>1e-8*max(np.linalg.norm(a),1e-100)


def test_nonzero_pole_rejected_instead_of_subtracted():
    # Deliberately use both scalar field and conjugate with a Majorana line.
    d=ma_toy(b=0.)
    v=d['vertices'][0]
    v['conjugate']=v['field']
    with pytest.raises(ContractError,match='UV pole'):
        independent_reference(d)


def test_non_symmetric_mass_rejected():
    d=dirac_toy()
    d['fermion_mass_gev']=encode([[1,2],[0,3]])
    with pytest.raises(ContractError,match='not complex symmetric'):
        independent_reference(d)


def test_ratio_function_near_equal_mass():
    z=np.array([1,1+1e-9,1-1e-9,1+1e-5,1e-4,1e5])
    f=ratio_function(z)
    assert np.all(np.isfinite(f))
    assert abs(f[0]-1)<1e-14
    assert abs(f[1]-(1+0.5e-9))<1e-14


def test_all16_synthetic_full_mass_runner_and_signed_labels(tmp_path, monkeypatch):
    import json
    from Numerical.lo_certification.direct_uv import native_uv_loop_reference as ref
    from Numerical.lo_certification.direct_uv.native_amplitude_contract import model_dictionary
    monkeypatch.setattr(ref,'source_audit',lambda root:{'status':'source_forms_verified','errors':[]})
    monkeypatch.setattr(ref,'native_portal_symbol_trace',lambda root:{'source_level_bare_to_indexed_portal_map_established':False})
    monkeypatch.setattr(ref,'optional_native_evidence',lambda root:{})
    directory=tmp_path/'inputs';directory.mkdir()
    for model in model_dictionary():
        name=model['model']
        ds1=model['SU2_dimensions']['S1'];ds2=model['SU2_dimensions']['S2'];df=model['SU2_dimensions']['F']
        maj=model['fermion_type']=='Majorana'
        nc=3*df;nf=nc if maj else 2*nc;ns=ds1+ds2
        M=np.eye(nc,dtype=complex)*900
        if not maj:
            z=np.zeros_like(M)
            M=np.block([[z,M],[M,z]])
        A=np.eye(ns,dtype=complex)*1e6
        A[0,ds1]=A[ds1,0]=800.
        Y=np.zeros((3,nf,ns),complex)
        Z=np.zeros((3,nf,ns),complex)
        Z[:,0 if maj else nc,0]=[.001,.002,.003]
        Y[:,0,ds1]=[.002,.001,-.001]
        obj={'model':name,'fermion_mass_gev':encode(M),
             'scalar_A_gev2':encode(A),'scalar_B_gev2':encode(np.zeros((ns,ns))),
             'renormalization_scale_gev':1e7,
             'vertices':[{'name':'y1_HC','conjugate':encode(Z)},
                         {'name':'y2_HC','field':encode(Y)}]}
        (directory/(name+'.json')).write_text(json.dumps(obj))
    result=ref.run_all(tmp_path,input_dir=directory,compare_production=False)
    assert result['model_count']==16
    assert len(result['reports_by_model'])==16
    assert result['physical_source_certified_models']==0
    assert all(r['uv_pole_relative']<1e-10 for r in result['models'])
    assert all(len(result['reports_by_model'][m['model']]['takagi_masses_ev_from_independent_singular_values'])==3
               for m in result['models'])
