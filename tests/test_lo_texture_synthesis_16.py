"""Batch-5: complex Weyl-texture factorisation, source guards, and 16-model toy integration.

The synthetic fixtures verify code logic only; they are not physical T3 results.
Set RUN_LO_TEXTURE_INTEGRATION=1 to run against the current frozen local checkout.
"""
from __future__ import annotations

from hashlib import sha256
import json
import os
from pathlib import Path

import numpy as np
import pytest

from Numerical.lo_certification.direct_uv.native_amplitude_contract import model_dictionary
from Numerical.lo_certification.direct_uv.certification import charge_bases
from Numerical.lo_patterns_16 import texture_synthesis as s


def pack(z): return s.encode(np.asarray(z,complex))

def put(p, obj):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def sha(p):return sha256(p.read_bytes()).hexdigest()


def toy_yukawas(seed=0):
    rng=np.random.default_rng(seed)
    y1=.002*(np.eye(3)+rng.normal(size=(3,3))*.45+1j*rng.normal(size=(3,3))*.1)
    y2=.002*(np.eye(3)+rng.normal(size=(3,3))*.45+1j*rng.normal(size=(3,3))*.1)
    return y1,y2


def demo_tensors(seed=0,nf=5,ns=4):
    rng=np.random.default_rng(seed)
    c1=rng.normal(size=(3,nf,ns))+1j*rng.normal(size=(3,nf,ns))
    c2=rng.normal(size=(3,nf,ns))+1j*rng.normal(size=(3,nf,ns))
    y1,y2=toy_yukawas(seed)
    z=np.einsum('ir,rfs->ifs',np.conj(y1),c1)
    y=np.einsum('ir,rfs->ifs',np.conj(y2),c2)
    data={'model':'A_+0','vertices':[{'name':'y1_HC','conjugate':pack(z)},
                                    {'name':'y2_HC','field':pack(y)}]}
    frozen={'base_state':{'ordinary':{'y1':pack(y1),'y2':pack(y2)}}}
    return data,frozen,c1,c2


def test_complex_yukawa_component_tensor_extracted_exactly():
    inp,frozen,c1,c2=demo_tensors(18)
    t=s.extract_source_tensors(inp,frozen,'A_+0')
    assert s.relative(t['C1'],c1)<1e-13
    assert s.relative(t['C2'],c2)<1e-13
    assert max(x['relative_reconstruction_error'] for x in t['tensor_reconstruction'].values())<1e-13


def test_complex_yukawa_replacement_bilinear_and_conjugation():
    inp,frozen,c1,c2=demo_tensors(29)
    source=s.extract_source_tensors(inp,frozen,'A_+0')
    a,b=toy_yukawas(30)
    got=s.substitute_textures(inp,source,a,b)
    assert s.relative(s.decode(got['vertices'][0]['conjugate'],'z'),
                      np.einsum('ir,rfs->ifs',np.conj(a),c1))<1e-12
    assert s.relative(s.decode(got['vertices'][1]['field'],'y'),
                      np.einsum('ir,rfs->ifs',np.conj(b),c2))<1e-12
    # Confusing transpose/conjugation changes the source for complex Yukawas.
    assert s.relative(s.decode(got['vertices'][0]['conjugate'],'z'),
                      np.einsum('ir,rfs->ifs',a,c1))>1e-2
    assert s.decode(inp['vertices'][0]['conjugate'],'z').shape==(3,5,4)


def test_single_yukawa_swap_is_not_both_yukawa_swap():
    inp,frozen,_,_=demo_tensors(7)
    origin=s.extract_source_tensors(inp,frozen,'A_+0')
    alt1,alt2=toy_yukawas(8)
    only1=s.substitute_textures(inp,origin,alt1,origin['y2'])
    only2=s.substitute_textures(inp,origin,origin['y1'],alt2)
    assert s.relative(s.decode(only1['vertices'][1]['field'],'z'),
                      s.decode(inp['vertices'][1]['field'],'z'))<1e-12
    assert s.relative(s.decode(only2['vertices'][0]['conjugate'],'z'),
                      s.decode(inp['vertices'][0]['conjugate'],'z'))<1e-12


def test_extract_source_vertex_name_or_orientation_wrong():
    inp,frozen,_,_=demo_tensors(2)
    inp['vertices'][0]['name']='y1'
    with pytest.raises(s.SynthesisError,match='orientation'):
        s.extract_source_tensors(inp,frozen,'A_+0')


def test_extract_source_wrong_model():
    inp,frozen,_,_=demo_tensors(2)
    with pytest.raises(s.SynthesisError,match='model label'):
        s.extract_source_tensors(inp,frozen,'B_-1')


def test_rank_deficient_yukawa_refused():
    inp,frozen,_,_=demo_tensors(3)
    y=s.decode(frozen['base_state']['ordinary']['y1'],'y1')
    y[2]=y[0]
    frozen['base_state']['ordinary']['y1']=pack(y)
    with pytest.raises(s.SynthesisError,match='rank/condition'):
        s.extract_source_tensors(inp,frozen,'A_+0')


def test_extremely_ill_conditioned_yukawa_refused():
    y=np.diag([1e-3,1e-9,1e-15]).astype(complex)
    inp=np.ones((3,3,2),complex)
    with pytest.raises(s.SynthesisError,match='rank/condition'):
        s._extract_channel_tensor(y,inp,'toy')


def test_nonfinite_yukawa_refused():
    inp,frozen,_,_=demo_tensors(3)
    y=s.decode(frozen['base_state']['ordinary']['y1'],'y1');y[0,0]=np.inf
    frozen['base_state']['ordinary']['y1']=pack(y)
    with pytest.raises(s.SynthesisError,match='nonfinite'):
        s.extract_source_tensors(inp,frozen,'A_+0')


def test_nonfinite_replacement_refused():
    inp,frozen,_,_=demo_tensors(3)
    t=s.extract_source_tensors(inp,frozen,'A_+0')
    wrong=np.zeros((3,3),complex);wrong[0,0]=np.nan
    with pytest.raises(s.SynthesisError,match='Nonfinite counterfactual'):
        s.substitute_textures(inp,t,wrong,t['y2'])


def test_original_vertices_not_modified():
    inp,frozen,_,_=demo_tensors(3)
    old=json.dumps(inp,sort_keys=True)
    t=s.extract_source_tensors(inp,frozen,'A_+0')
    s.substitute_textures(inp,t,*toy_yukawas(55))
    assert json.dumps(inp,sort_keys=True)==old


def test_scalar_mixing_rescaled_not_diagonal():
    A=np.diag([8000.,9000.,10000.,11000.]).astype(complex)
    A[2,0]=11+2j;A[0,2]=11-2j
    A[3,0]=-3j;A[0,3]=3j
    data={'scalar_A_gev2':pack(A)}
    new=s.scalar_mixing_scale(data,1,2.2)
    B=s.decode(new['scalar_A_gev2'],'A')
    assert np.diag(B)==pytest.approx(np.diag(A))
    assert B[2,0]==pytest.approx(A[2,0]*2.2)
    assert B[0,2]==pytest.approx(A[0,2]*2.2)
    assert s.relative(A,s.decode(data['scalar_A_gev2'],'A'))==0


def test_scalar_negative_mixing_factor_refused():
    A={'scalar_A_gev2':pack(np.eye(4))}
    with pytest.raises(s.SynthesisError,match='positive'):
        s.scalar_mixing_scale(A,1,-1)


def test_scalar_nonhermitian_refused():
    A=np.eye(4,dtype=complex);A[1,0]=2j
    with pytest.raises(s.SynthesisError,match='Hermitian'):
        s.scalar_mixing_scale({'scalar_A_gev2':pack(A)},1,2)


def test_complex_phasors_exact_coherence_and_opposing_pair():
    from Numerical.lo_patterns_16.flavour_geometry import geometry
    A=np.diag([1.1+.4j,2.,3.]);B=np.diag([-1.07-.4j,3.,4.]);C=np.diag([.004,2.,3.]);M=A+B+C
    groups=[{'matrix_ev':x,'mass_gev':10000.*(k+1),'multiplicity':1}
            for k,x in enumerate([A,B,C])]
    got=geometry(M,groups,'toy')
    obj={'m1_ev':got['masses_ev'][0],
         'm1_coherence_C':got['modes'][0]['coherence_fraction'],
         'group_strength_ev':got['fermi_group_strength_ev'],
         'phasors':got['projections'][:3]}
    exp={'lightest_mode':got['modes'][0]}
    info=s.phasor_explanation(obj,exp,'toy')
    assert info['cancelled_fraction']>0.98
    assert info['strongest_destructive_pair']['phase_cosine']<-0.9
    assert abs(info['phasor_sum_ev']['real']-got['masses_ev'][0])<1e-10


def test_phasor_conservation_negative_control():
    from Numerical.lo_patterns_16.flavour_geometry import geometry
    H=np.diag([1.,2.,3.]);J=np.diag([.3,1.,1.])
    x=geometry(H+J,[{'matrix_ev':H,'mass_gev':100,'multiplicity':1},
                    {'matrix_ev':J,'mass_gev':200,'multiplicity':1}], 'toy')
    obj={'m1_ev':x['masses_ev'][0],'m1_coherence_C':x['modes'][0]['coherence_fraction'],
         'group_strength_ev':x['fermi_group_strength_ev'],
         'phasors':x['projections'][:2]}
    obj['phasors'][0]['amplitude_ev']['real']*=1.1
    with pytest.raises(s.SynthesisError,match='do not sum'):
        s.phasor_explanation(obj,{'lightest_mode':x['modes'][0]},'toy')


def test_mode_sum_is_not_singular_value_sum():
    from Numerical.lo_patterns_16.flavour_geometry import geometry
    A=np.diag([1.,2.,3.]);B=np.diag([-1.,1.,1.]);H=A+B
    answer=geometry(H,[{'matrix_ev':A,'mass_gev':100,'multiplicity':1},
                       {'matrix_ev':B,'mass_gev':200,'multiplicity':1}], 'toy')
    assert answer['masses_ev'][0]==0
    assert not np.allclose(np.sort(np.linalg.svd(H,compute_uv=False)),
                           np.sort(np.linalg.svd(A,compute_uv=False))+
                           np.sort(np.linalg.svd(B,compute_uv=False)))


def test_path_traversal_rejected(tmp_path):
    with pytest.raises(s.SynthesisError,match='escapes checkout'):
        s._resolve_local(tmp_path.resolve(),'../../steal.json')


def test_existing_output_dir_refusal_early(tmp_path):
    root=tmp_path.resolve()
    target=root/s.TARGET;target.mkdir(parents=True)
    put(target/'saved.json',{'important':True})
    # A source runner should fail before reading missing source data.
    with pytest.raises(s.SynthesisError,match='already populated'):
        s.analyse(root,make_plots=False)
    assert json.loads((target/'saved.json').read_text())=={'important':True}


def test_disallow_outside_output_directory(tmp_path):
    with pytest.raises(s.SynthesisError,match='not a contained'):
        s.analyse(tmp_path,output_dir=tmp_path/'outside',make_plots=False)


def test_read_missing_file_refused(tmp_path):
    with pytest.raises(s.SynthesisError,match='missing'):
        s.read(tmp_path/'missing.json')


def test_decode_shape_refused():
    with pytest.raises(s.SynthesisError,match='bad array shape'):
        s.decode(pack(np.eye(3)),'toy',(2,2))


def _synthetic_model_data(item,seed):
    """Build a 3-heavy-generation source-conditional Weyl toy with fixed CG maps."""
    label=item['model'];kind=label[0];alpha=int(label[2:])
    d1=item['SU2_dimensions']['S1'];d2=item['SU2_dimensions']['S2'];df=item['SU2_dimensions']['F']
    maj=item['fermion_type']=='Majorana'
    nf=(3*df if maj else 6*df);ns=d1+d2
    mass=np.array([120.,175.,260.])
    D=np.diag(np.repeat(mass,df))
    if maj:M=D.astype(complex)
    else:
        zero=np.zeros_like(D)
        M=np.block([[zero,D],[D,zero]]).astype(complex)
    A=np.diag([1.8e5]*d1+[2.3e5]*d2).astype(complex)
    charges=charge_bases(kind,alpha)['scalar']
    pairs=[(i,j) for i in range(d1) for j in range(d1,ns) if abs(charges[i]-charges[j])<1e-10]
    assert len(pairs)>=1
    C1=np.zeros((3,nf,ns),complex);C2=np.zeros_like(C1)
    for ch,(i,j) in enumerate(pairs):
        A[j,i]=(650.+ch*10.);A[i,j]=A[j,i]
        for r in range(3):
            index=r*df
            C1[r,index if maj else 3*df+index,i]=(.8+.04*r)*(.7+.1*ch)
            C2[r,index,j]=(.9-.04*r)*(.7+.1*ch)
    y1,y2=toy_yukawas(seed)
    z=np.einsum('ir,rfs->ifs',np.conj(y1),C1)
    y=np.einsum('ir,rfs->ifs',np.conj(y2),C2)
    inp={'model':label,'renormalization_scale_gev':1e7,
         'fermion_mass_gev':pack(M),'scalar_A_gev2':pack(A),
         'scalar_B_gev2':pack(np.zeros_like(A)),
         'vertices':[{'name':'y1_HC','conjugate':pack(z)},
                     {'name':'y2_HC','field':pack(y)}]}
    frozen={'base_state':{'ordinary':{'y1':pack(y1),'y2':pack(y2),
                   'lambdaT3':{'real':.05 if kind=='E' else .1,'imag':0}}}}
    return inp,frozen


def fixture_16(root):
    """Consistent synthetic B1/B3/B4 source chain, NO production source mutation."""
    from Numerical.lo_patterns_16.analytic_explanation import reconstruct
    names=[x['model'] for x in model_dictionary()]
    basis=[];d3=[];d4=[];manifest3={};manifest4={}
    for i,item in enumerate(model_dictionary()):
        name=item['model'];inp,frozen=_synthetic_model_data(item,i+1)
        input_file=root/'output/general_observable/lo_direct_one_loop/provisional_inputs'/f'{name}.json'
        bench_file=root/'Numerical/benchmarks/fitted_reference'/f'{name}.json'
        put(input_file,inp);put(bench_file,frozen)
        src={'input':str(input_file.relative_to(root)).replace('\\','/'),
             'benchmark':str(bench_file.relative_to(root)).replace('\\','/'),
             'input_sha256':sha(input_file),'benchmark_sha256':sha(bench_file)}
        rec=reconstruct(inp,name[0],int(name[2:]),item['SU2_dimensions']['S1'])
        model_file=root/s.STUDY/'models'/f'{name}.json'
        put(model_file,{'model':name,'status':'conditional_uv_evaluated',
                        'physical_normalisation_certified':False,
                        'mnu_gev':pack(rec['exact']),'source':src})
        basis.append({'model':name})
        groups=[{'fermion_mass_gev':float(np.sqrt(h['m2'])),
                 'eigenspace_multiplicity':h['dimension'],
                 'matrix_ev':pack(h['matrix']*1e9),
                 'matrix_norm_ev':float(np.linalg.norm(h['matrix'])*1e9)}
                for h in rec['heavy_groups']]
        d3.append({'model':name,'physical_normalisation_certified':False,
                   'heavy_mass_groups':groups})
        info=s.mode_summary(rec['exact'],rec['heavy_groups'],name)
        d4.append({'model':name,'m1_over_m3':info['m1_over_m3'],
                   'lightest_mode':{'coherence_fraction':info['m1_coherence_C']}})
        put(root/s.BATCH03/'models'/f'{name}.json',
            {'exact_closed_kernel_mnu_gev':pack(rec['exact'])})
        manifest3[name]={'model_detail_sha256':sha(model_file),
                         'benchmark_sha256':src['benchmark_sha256']}
        manifest4[name]={'batch1_model_sha256':sha(model_file),
                         'candidate_weyl_input_sha256':src['input_sha256']}
    summary={'format':'t3-conditional-pattern-study-all16-v1',
             'physical_normalisation_certified':False,
             'counts':{'conditional_evaluated':16,'blocked':0,'physical_source_certified':0},
             'provenance':{'revision':{'local_head':'synthetic'},'source_audit':{'source_files':{}},
                           'source_audit_sha256':{}},'models':basis}
    put(root/s.STUDY/'all16_summary.json',summary)
    shasum=sha(root/s.STUDY/'all16_summary.json')
    put(root/s.BATCH03/'batch03_analytic_explanation.json',
        {'format':'t3-chatA-batch03-analytic-scalar-kernel-v1','model_count':16,
         'physical_normalisation_certified':False,'source_summary_sha256':shasum,'models':d3})
    sha3=sha(root/s.BATCH03/'batch03_analytic_explanation.json')
    put(root/s.BATCH03/'SOURCE_MANIFEST.json',
        {'source_summary_sha256':shasum,'physical_source_certified':False,'models':manifest3})
    put(root/s.BATCH04/'batch04_flavour_analysis.json',
        {'format':'t3-chatA-batch04-flavour-geometry-v1','model_count':16,
         'physical_normalisation_certified':False,
         'source_sha256':{'source_run_sha256':shasum,'source_batch03_sha256':sha3},
         'models':d4})
    put(root/s.BATCH04/'SOURCE_MANIFEST.json',
        {'source_batch03_sha256':sha3,'physical_normalisation_certified':False,
         'models':manifest4})
    return root


def test_full_synthetic_all16_control_matrix_and_preserve_sources(tmp_path):
    root=fixture_16(tmp_path.resolve())
    before={str(p):sha(p) for p in root.rglob('*') if p.is_file()}
    result=s.analyse(root,make_plots=False)
    assert result['model_count']==16
    assert result['physical_normalisation_certified'] is False
    assert all(len(m['common_mixing_controlled'])==4 for m in result['models'])
    assert sum(len(m['C_family_swaps']) for m in result['models'])==21
    assert len(result['C_texture_interpolations'])==22
    out=root/s.TARGET
    assert (out/'BATCH05_ANALYTICAL_SYNTHESIS.md').is_file()
    assert len(json.loads((out/'batch05_synthesis.json').read_text())['models'])==16
    assert all(sha(Path(k))==v for k,v in before.items())
    assert (out/'batch05_C_texture_interpolation.csv').is_file()


def test_synthetic_source_hash_tamper_refused_before_write(tmp_path):
    root=fixture_16(tmp_path.resolve())
    infile=root/'output/general_observable/lo_direct_one_loop/provisional_inputs/C_+1.json'
    x=read_json(infile);x['model']='C_-1';put(infile,x)
    with pytest.raises(Exception,match='SHA|hash'):
        s.analyse(root,make_plots=False)
    assert not (root/s.TARGET).exists()


def test_synthetic_false_physical_certification_refused(tmp_path):
    root=fixture_16(tmp_path.resolve())
    m=root/s.BATCH04/'batch04_flavour_analysis.json'
    x=read_json(m);x['physical_normalisation_certified']=True;put(m,x)
    with pytest.raises(s.SynthesisError,match='certification'):
        s.analyse(root,make_plots=False)
    assert not (root/s.TARGET).exists()


def test_synthetic_missing_reference_model_refused(tmp_path):
    root=fixture_16(tmp_path.resolve())
    m=root/s.BATCH04/'batch04_flavour_analysis.json'
    x=read_json(m);x['models']=x['models'][:-1];put(m,x)
    with pytest.raises(s.SynthesisError,match='model catalogue'):
        s.analyse(root,make_plots=False)
    assert not (root/s.TARGET).exists()


def test_synthetic_model_group_amplitude_tamper_refused(tmp_path):
    root=fixture_16(tmp_path.resolve())
    m=root/s.BATCH03/'batch03_analytic_explanation.json'
    x=read_json(m);x['models'][0]['heavy_mass_groups'][0]['matrix_ev']['real'][0][0]*=1.2;put(m,x)
    # SHA lineage is checked before any per-model tensor recalculation.
    with pytest.raises(s.SynthesisError,match='lineage SHA'):
        s.analyse(root,make_plots=False)
    assert not (root/s.TARGET).exists()


def read_json(p):return json.loads(p.read_text())


@pytest.mark.skipif(os.environ.get('RUN_LO_TEXTURE_INTEGRATION')!='1',
                    reason='Requires actual frozen 16-model checkout and all Batch 1/3/4 outputs')
def test_opt_in_real_local_frozen_full_integration():
    root=Path(__file__).resolve().parents[1]
    result=s.analyse(root,output_dir=root/s.BASE/'analysis_batch05_pytest_integration_01',make_plots=False)
    assert result['model_count']==16
    assert all(m['physical_normalisation_certified'] is False for m in result['models'])
    assert all(m['tensor_reconstruction']['y1_HC']['relative_reconstruction_error']<5e-9
               for m in result['models'])
