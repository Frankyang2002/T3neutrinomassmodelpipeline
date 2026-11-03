"""Batch-3 scalar spectral identity and 16-case synthetic report; not physical certification."""
import json
from hashlib import sha256
from pathlib import Path

import numpy as np
import pytest

from Numerical.lo_certification.direct_uv.native_amplitude_contract import model_dictionary
from Numerical.lo_certification.direct_uv.native_uv_loop_reference import encode, decode
from Numerical.lo_certification.direct_uv.certification import charge_bases
from Numerical.lo_patterns_16 import analytic_explanation as a


def pack(z):return encode(np.asarray(z,dtype=complex))

def demo(kind='A',alpha=0, *, complex_phase=True):
    """No physical particle content: a minimal, explicitly provisional Weyl toy."""
    cls={x['model']:x for x in model_dictionary()}
    label=f'{kind}_{alpha:+d}'
    info=cls[label]
    d1=info['SU2_dimensions']['S1'];d2=info['SU2_dimensions']['S2'];df=info['SU2_dimensions']['F']
    maj=info['fermion_type']=='Majorana'
    nf=3*df*(1 if maj else 2);ns=d1+d2
    m0=230.
    if maj:M=np.eye(nf,dtype=complex)*m0
    else:
        half=nf//2
        M=np.block([[np.zeros((half,half)),np.eye(half)*m0],
                    [np.eye(half)*m0,np.zeros((half,half))]]).astype(complex)
    A=np.diag([4.2e4]*d1+[5.3e4]*d2).astype(complex)
    qs=charge_bases(kind,alpha)['scalar']
    pairs=[(s,t) for s in range(d1) for t in range(d1,ns) if abs(qs[s]-qs[t])<1e-10]
    pairs=pairs[:2 if kind=='E' else 1]
    Z=np.zeros((3,nf,ns),complex);Y=np.zeros_like(Z)
    for idx,(s,t) in enumerate(pairs):
        d=(1200+170*idx)*(np.exp(.4j) if complex_phase else 1.)
        A[t,s]=d;A[s,t]=d.conjugate()
        Z[:,(0 if maj else nf//2),s]=np.array([.013+.004j,.009-.002j,.003+.001j])*(idx+1)
        Y[:,0,t]=np.array([.008-.003j,.011+.005j,-.004+.002j])
    return {'model':label,'fermion_mass_gev':pack(M),'scalar_A_gev2':pack(A),
            'scalar_B_gev2':pack(np.zeros_like(A)),
            'vertices':[{'name':'y1_HC','conjugate':pack(Z)},
                        {'name':'y2_HC','field':pack(Y)}],
            'renormalization_scale_gev':1e7},pairs


def exact_eigenbasis_in_toy(data):
    """Independent scalar matrix-function check for degenerate heavy M²=mF² I.

    Uses direct eigensums, not the tested divided-difference expression.
    """
    M=decode(data['fermion_mass_gev'],'M');A=decode(data['scalar_A_gev2'],'A')
    z=decode(data['vertices'][0]['conjugate'],'z')
    y=decode(data['vertices'][1]['field'],'y')
    val,U=np.linalg.eigh(A)
    m2=float((M@M.conj().T)[0,0].real)
    mm=np.zeros((3,3),complex)
    for r,s2 in enumerate(val):
        zr=np.einsum('ais,s->ai',z,U[:,r].conj())
        yr=np.einsum('ais,s->ai',y,U[:,r])
        t=zr@M.conj().T@yr.T
        mm+=a.LOOP*a.b0(m2,float(s2))*(t+t.T)
    return mm


def test_scalar_divided_difference_near_equal_scalar_and_heavy_masses():
    m2=10000.
    assert a.divdiff(m2,m2,m2)==pytest.approx(-.5/m2,rel=1e-12)
    assert a.divdiff(m2,1.5e4,1.5e4*(1+1e-11))==pytest.approx(a.b0_derivative(m2,1.5e4),rel=1e-8)


def test_exact_kernel_reconstructs_independently_diagonalised_complex_scalar():
    for label,alpha in [('A',0),('B',-1),('C',-1),('D',0),('E',0)]:
        data,_=demo(label,alpha)
        d1=next(x for x in model_dictionary() if x['model']==data['model'])['SU2_dimensions']['S1']
        result=a.reconstruct(data,label,alpha,d1)
        other=exact_eigenbasis_in_toy(data)
        assert a.relative(result['exact'],other)<1e-9, data['model']
        assert np.linalg.norm(result['exact'])>0


def test_rephased_complex_delta_not_conjugated():
    data,pairs=demo('A',0)
    A=decode(data['scalar_A_gev2'],'A')
    val,U=np.linalg.eigh(A)
    d1=1
    got=a.reconstruct(data,'A',0,d1)['exact']
    assert a.relative(got,exact_eigenbasis_in_toy(data))<1e-9
    # Swapping delta with its conjugate without transforming Yukawas is a
    # genuine phase mistake and changes the neutrino-mass matrix.
    s,t=pairs[0]
    A[t,s]=A[t,s].conjugate(); A[s,t]=A[t,s].conjugate()
    changed={**data,'scalar_A_gev2':pack(A)}
    bad=a.reconstruct(changed,'A',0,d1)['exact']
    assert a.relative(got,bad)>0.1


def test_lambda_sign_flip_is_exactly_odd_and_leading_only_approximate():
    data,_=demo('A',0)
    A=decode(data['scalar_A_gev2'],'A')
    original=a.reconstruct(data,'A',0,1)
    flipped=A.copy();flipped[:1,1:]*=-1;flipped[1:,:1]*=-1
    new=a.reconstruct({**data,'scalar_A_gev2':pack(flipped)},'A',0,1)
    assert a.relative(original['exact'],-new['exact'])<1e-12
    assert original['linear_relative_error']>0
    zero=A.copy();zero[:1,1:]=0;zero[1:,:1]=0
    z=a.reconstruct({**data,'scalar_A_gev2':pack(zero)},'A',0,1)
    assert np.linalg.norm(z['exact'])<1e-30


def test_rejects_unsupported_scalar_B_or_extra_yukawa():
    data,_=demo('A',0)
    B=np.eye(4)*50
    with pytest.raises(a.ExplanationError,match='Scalar B is nonzero'):
        a.reconstruct({**data,'scalar_B_gev2':pack(B)},'A',0,1)
    extra={**data,'vertices':data['vertices']+[{'name':'extra','field':pack(np.zeros((3,12,4)))}]}
    with pytest.raises(a.ExplanationError,match='New vertex sources'):
        a.reconstruct(extra,'A',0,1)


def test_rejects_wrong_charge_cross_mixing():
    data,_=demo('A',0)
    mat=decode(data['scalar_A_gev2'],'A')
    # S1 has Q=0, S2 has Q=+2,+1,0; index 1 cannot mix with S1.
    mat[1,0]=100;mat[0,1]=100
    with pytest.raises(a.ExplanationError,match='charge conservation'):
        a.reconstruct({**data,'scalar_A_gev2':pack(mat)},'A',0,1)


def put(path,record):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(record),encoding='utf-8')


def fingerprint(path):return sha256(path.read_bytes()).hexdigest()


def synthetic_project(tmp_path,monkeypatch):
    root=tmp_path/'project';root.mkdir()
    study=root/a.STUDY
    prev=root/a.BATCH02
    base=a.INPUTS
    audit_file=root/a.CG_AUDIT
    names=[m['model'] for m in model_dictionary()]
    frozen_rows=[];report2=[];audit_rows=[];weights={}
    for model in names:
        family,alpha=model[0],int(model[2:]);item=next(x for x in model_dictionary() if x['model']==model)
        dims=item['SU2_dimensions']
        data,pairs=demo(family,alpha)
        inp=root/base/f'{model}.json';put(inp,data)
        repath=(f'Numerical/benchmarks/fitted_reference/T3_dS1_{dims["S1"]}_dS2_{dims["S2"]}'
                f'_dF_{dims["F"]}_alpha_{"m" if alpha<0 else "p"}{abs(alpha)}.json')
        bench=root/repath
        y=np.eye(3,dtype=complex)*1e-3
        put(bench,{'base_state':{'ordinary':{'y1':pack(y),'y2':pack(y),
                  'lambdaT3':{'real':.1,'imag':0}}}})
        candidate=a.reconstruct(data,family,alpha,dims['S1'])
        masses=np.sort(np.linalg.svd(candidate['exact'],compute_uv=False)*1e9).tolist()
        record={'model':model,'status':'conditional_uv_evaluated',
                'physical_normalisation_certified':False,
                'source':{'benchmark':repath,'benchmark_sha256':fingerprint(bench),
                          'input':str((base/f'{model}.json').as_posix()),'input_sha256':fingerprint(inp)},
                'mnu_gev':pack(candidate['exact']),'takagi_masses_ev_sorted':masses}
        put(study/'models'/f'{model}.json',record)
        frozen_rows.append({'model':model,'status':'conditional_uv_evaluated',
                            'benchmark_sha256':fingerprint(bench),
                            'component_input_sha256':fingerprint(inp)})
        projections={}
        for factor,key in ((2.,'lambda_t2_linear_matrix_residual'),
                           (.5,'lambda_t05_linear_matrix_residual'),
                           (-1.,'lambda_sign_flip_matrix_residual')):
            changed=decode(data['scalar_A_gev2'],'A').copy()
            changed[:dims['S1'],dims['S1']:]*=factor
            changed[dims['S1']:,:dims['S1']]*=factor
            ans=a.reconstruct({**data,'scalar_A_gev2':pack(changed)},family,alpha,dims['S1'])['exact']
            projections[key]=np.linalg.norm(ans-factor*candidate['exact']) / np.linalg.norm(candidate['exact'])
        report2.append({'model':model,'physical_normalisation_certified':False,
                        'benchmark_sha256':fingerprint(bench),
                        'component_input_sha256':fingerprint(inp),
                        **dict(zip(['m1_ev','m2_ev','m3_ev'],masses)),**projections})
        ch=[]
        for p in candidate['channels']:
            if abs(p['delta_gev2'])==0:continue
            q=p['charge'];ch.append({'charge_numeric':q,'nonzero':True,
                                    'complex_weight_numeric':{'real':1.,'imag':0.}})
        weights[model]={'channels':ch,'relative_physical_normalization_certified':False}
    put(audit_file,{'status':'charge_resolved_group_coefficients','count':16})
    summary={'format':'t3-conditional-pattern-study-all16-v1',
             'physical_normalisation_certified':False,
             'counts':{'conditional_evaluated':16,'blocked':0,'physical_source_certified':0},
             'models':frozen_rows,
             'provenance':{'revision':{'local_head':'synthetic'},
               'source_audit':{'source_files':{}},
               'source_audit_sha256':{'charge_resolved_loop_weights.json':fingerprint(audit_file)}}}
    put(study/'all16_summary.json',summary)
    put(study/'channel_weights.json',{'status':'native_group_CG_candidate_only_not_physical_normalization',
                                       'models':weights})
    put(prev/'batch02_analysis.json',{'format':'t3-chatA-readonly-batch02-numerical-patterns-v1',
      'physical_normalisation_certified':False,
      'counts':{'models':16,'source_scan_rows':288,'physical_source_certified':0},
      'input_run_sha256':fingerprint(study/'all16_summary.json'),
      'models':report2})
    # This is a synthetic fixture; source SHA check of a real checkout is
    # separately verified by the existing Batch-2 tests.
    monkeypatch.setattr(a,'_source_hash_guard',lambda *v,**k:None)
    return root


def test_all16_synthetic_exact_reports_are_read_only_and_uncertified(tmp_path,monkeypatch):
    root=synthetic_project(tmp_path,monkeypatch)
    initial={str(p.relative_to(root)):fingerprint(p) for p in root.rglob('*') if p.is_file()}
    answer=a.explain(root)
    assert answer['model_count']==16
    assert answer['physical_normalisation_certified'] is False
    assert max(r['exact_reconstruction_relative_error'] for r in answer['models'])<1e-10
    assert all(r['m3_ev']>0 for r in answer['models'])
    assert all(len(r['channels'])==2 for r in answer['models'] if r['family']=='E')
    assert (root/a.BATCH03/'BATCH03_ANALYTICAL_REPORT.md').is_file()
    assert (root/a.BATCH03/'HANDOFF.json').is_file()
    after={str(p.relative_to(root)):fingerprint(p) for p in root.rglob('*') if p.is_file()}
    for path,h in initial.items():assert after[path]==h


def test_blocks_output_overwrite(tmp_path,monkeypatch):
    root=synthetic_project(tmp_path,monkeypatch)
    a.explain(root)
    with pytest.raises(a.ExplanationError,match='already populated'):
        a.explain(root)


def test_blocks_changed_batch02_provenance(tmp_path,monkeypatch):
    root=synthetic_project(tmp_path,monkeypatch)
    p=root/a.BATCH02/'batch02_analysis.json'
    b=json.loads(p.read_text());b['input_run_sha256']='0'*64;put(p,b)
    with pytest.raises(a.ExplanationError,match='SHA256 does not match'):
        a.explain(root)
    assert not (root/a.BATCH03).exists()


def test_blocks_false_certificate(tmp_path,monkeypatch):
    root=synthetic_project(tmp_path,monkeypatch)
    p=root/a.BATCH02/'batch02_analysis.json'
    b=json.loads(p.read_text());b['physical_normalisation_certified']=True;put(p,b)
    with pytest.raises(a.ExplanationError,match='False physical certification'):
        a.explain(root)
    assert not (root/a.BATCH03).exists()


def test_blocks_wrong_native_weight_support(tmp_path,monkeypatch):
    root=synthetic_project(tmp_path,monkeypatch)
    p=root/a.STUDY/'channel_weights.json';w=json.loads(p.read_text())
    w['models']['A_+0']['channels'][0]['charge_numeric']=17
    put(p,w)
    with pytest.raises(a.ExplanationError,match='native CG charge support'):
        a.explain(root)
    assert not (root/a.BATCH03).exists()


def test_uses_direct_analytic_source_propagator_not_production_takagi(monkeypatch):
    from Numerical.direct_one_loop import core
    data,_=demo('B',-1)
    def forbidden(*a,**k):
        raise AssertionError('Production full loop/Takagi must not be invoked')
    monkeypatch.setattr(core,'evaluate',forbidden)
    monkeypatch.setattr(core,'takagi',forbidden)
    assert np.linalg.norm(a.reconstruct(data,'B',-1,2)['exact'])>0


def test_complex_scalar_and_heavy_field_rephasing_keeps_neutrino_matrix():
    from copy import deepcopy
    for kind,alpha in [('A',0),('C',-1),('E',-2)]:
        data,_=demo(kind,alpha)
        dim=next(x for x in model_dictionary() if x['model']==data['model'])['SU2_dimensions']['S1']
        base=a.reconstruct(data,kind,alpha,dim)['exact']
        m=decode(data['fermion_mass_gev'],'M')
        aa=decode(data['scalar_A_gev2'],'A')
        nf,ns=len(m),len(aa)
        df=np.exp(1j*np.arange(nf)*.171)
        ds=np.exp(1j*np.arange(ns)*.293)
        out=deepcopy(data)
        out['fermion_mass_gev']=pack(df[:,None]*m*df[None,:])
        out['scalar_A_gev2']=pack(ds.conj()[:,None]*aa*ds[None,:])
        for index,row in enumerate(out['vertices']):
            if 'field' in row:
                orig=decode(data['vertices'][index]['field'],'Y')
                row['field']=pack(orig*df[None,:,None]*ds[None,None,:])
            if 'conjugate' in row:
                orig=decode(data['vertices'][index]['conjugate'],'Z')
                row['conjugate']=pack(orig*df[None,:,None]*ds.conj()[None,None,:])
        got=a.reconstruct(out,kind,alpha,dim)['exact']
        assert a.relative(base,got)<1e-9, data['model']


def test_wrong_loop_factor_and_omitted_pair_transpose_detected():
    data,_=demo('A',0)
    result=a.reconstruct(data,'A',0,1)['exact']
    assert a.relative(result,2*result)>0.4  # catches incorrect loop multiplicity
    # Symmetric mass cannot equal just the untransposed source cross term
    Z=decode(data['vertices'][0]['conjugate'],'Z')
    Y=decode(data['vertices'][1]['field'],'Y')
    z=Z[:,6,0]
    y=Y[:,0,3]
    one=np.outer(z,y)
    assert np.linalg.norm(one-one.T)>1e-5*np.linalg.norm(one)


def test_divided_difference_against_feynman_parameter_integral():
    from scipy.integrate import quad
    for m2,a2,b2 in [(2.3e10,1e6,1.3e6),(4.5e4,5.0e4,8.1e4),
                      (4.5e4,4.5e4,4.5e4)]:
        # B0 difference quotient = -int dx [ln(D_b)-ln(D_a)]/(b-a).
        if a2==b2:
            ref=-quad(lambda x: (1-x)/(x*m2+(1-x)*a2),0,1,
                      epsabs=1e-14,epsrel=1e-12)[0]
        else:
            ref=-quad(lambda x: np.log1p((1-x)*(b2-a2)/
                     (x*m2+(1-x)*a2))/(b2-a2),0,1,
                     epsabs=1e-14,epsrel=1e-12)[0]
        assert a.divdiff(m2,a2,b2)==pytest.approx(ref,rel=3e-9,abs=1e-17)
