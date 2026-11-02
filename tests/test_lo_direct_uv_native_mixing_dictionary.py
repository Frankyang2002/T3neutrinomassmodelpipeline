"""Synthetic regression tests; NOT claims of a live Matchete execution."""
from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import tempfile

import numpy as np
import pytest

from Numerical.lo_certification.direct_uv.native_mixing_dictionary import (
    CLASSES, REPRESENTATIVES, SOURCE_FILES, AUDIT_DIR, MixingAuditError, all_sixteen,
    audit, charge, complex_number, finite_hessian, generator_rotations,
    mixing_block, models, native_models, native_tensor, read,
    rotated_tensor, source_guard, source_tensor,
)

FLAGS = ('T3_EXPERIMENTAL_GENERAL_FERMIONS',
         'T3_EXPERIMENTAL_GENERAL_QUADRATICS',
         'T3_EXPERIMENTAL_GENERAL_QUARTICS',
         'T3_EXPERIMENTAL_FLAVOR_MASSES')


def seed_native():
    rows = []
    tensors = {}
    for cls,(d1,d2,df,_alphas) in CLASSES.items():
        alpha=REPRESENTATIVES[cls]
        t=np.zeros((2,2,d1,d2),complex)
        # H0 H0 neutrality condition independent of alpha.
        s=0
        tt=int(1-(d1-d2)/2+s)
        assert 0<=tt<d2
        t[1,1,s,tt]=1.5
        if d1==d2==3:
            t[1,1,1,2]=-0.5
        tensors[cls]=t
        reduced=[d for d in (2,2,d1,d2) if d>1]
        flat=t.reshape(-1)
        rows.append({'class':cls,'model':f'{cls}_{alpha:+d}','alpha':alpha,
            'dS1':d1,'dS2':d2,'dF':df,
            'native_build_status':'Success','native_full_lagrangian_validation':True,
            'mixing_operator_name':'T3ScalarMix_GenericCG','mixing_coupling':'lambdaT3',
            'native_nonhc_operator_present':True,'native_hc_operator_present':True,
            'Higgs_fields_in_nonhc_operator':'H H',
            'scalar_fields_in_nonhc_operator':'S1 Bar[S2]',
            'native_nonhc_LInt_coefficient':'+lambdaT3 * T3MixCG (no 1/2)',
            'physical_source_to_Weyl_certified':False,
            'native_cg':{'symbol':'T3MixCG','tensor_accessor':'Matchete`GetCGTensor',
                'dimensions':reduced,
                'full_dimensions_including_singlets':[2,2,d1,d2],
                'registered_in_assembled_interaction':True,
                'components_inputform':[str(float(z.real)) for z in flat],
                'components_real':flat.real.tolist(),
                'components_imag':flat.imag.tolist()}})
    report={'format':'t3-ws4-native-mixing-tensor-five-classes-v1',
      'status':'five_classes_native_registered_mixing_tensors_extracted',
      'restricted_environment':{'isolated_to_wolframscript_process':True,
         'effective':{k:'0' for k in FLAGS}},
      'native_rgbeta_mixing_coupling_dictionary_certified':False,
      'native_spinor_and_mass_dictionary_certified':False,
      'matching_or_running_performed':False,
      'models':rows,'native_source_sha256':{}}
    return report,tensors


def write_json(path:Path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value),encoding='utf-8')


def seed_root(root:Path):
    report,tensors=seed_native()
    contents={
      'Lagrangian/interactions/T3Topology.wl':
        '"Expression" -> PlusHc[ lambdaT3[] H[i] H[j] s1 s2 CG[T3MixCG, labels] ]',
      'Lagrangian/model/SU2Invariants.wl':
        '{2, 2, d1, d2} {False, False, False, True} {1, 2}',
      'scripts/BuildT3DirectOneLoopInputs.py':
        "cg['T3MixCG'][1,1,s,t] Vcoeff=-lam*v*v*cg['T3MixCG'][1,1,s,t]/2 "
        'A[ds1+t,s]+=Vcoeff A[s,ds1+t]+=Vcoeff.conjugate()',
    }
    from Numerical.lo_certification.direct_uv.native_mixing_dictionary import digest
    for rel in SOURCE_FILES:
        path=root/rel
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(contents.get(rel,'verified synthetic source'),encoding='utf-8')
        report['native_source_sha256'][rel]=digest(path)
    gen={'status':'generator_algebra_verified',
         'representations':[{'status':'generator_algebra_verified','dimension':d,
            'U_columns_are_weight_eigenvectors':np.eye(d,dtype=int).tolist()}
            for d in (2,3)]}
    cases=[]; cov=[]
    for name,key,dims in models():
        d1,d2,df,a=dims
        native=tensors[name[0]]
        components=[{'indices':[i+1,j+1,s+1,t+1],'cg':str(z.real)}
            for (i,j,s,t),z in np.ndenumerate(native) if abs(z)>1e-14]
        record={'model':{'dS1':d1,'dS2':d2,'dF':df,'alpha':a},
                'vertices':{'T3MixCG':{'dimensions':[2,2,d1,d2],
                                      'components':components}}}
        cases.append(record)
        cov.append({'model':record['model'],'all_covariant':True})
        lam=complex(0.31,-0.47)
        vev=246.
        A=np.eye(d1+d2,dtype=complex)*1e6
        expected=mixing_block(native,lam,vev)
        A[d1:,:d1]=expected
        A[:d1,d1:]=expected.conj().T
        encoded=lambda v: {'real':v.real.tolist(),'imag':v.imag.tolist()}
        write_json(root/'output/general_observable/lo_direct_one_loop/provisional_inputs'/f'{name}.json',
                   {'model':name,'scalar_A_gev2':encoded(A),
                    'scalar_B_gev2':encoded(np.zeros_like(A))})
        write_json(root/'Numerical/benchmarks/fitted_reference'/f'{key}.json',
                   {'representation':{'d_s1':d1,'d_s2':d2,'d_f':df,'alpha':a},
                    'numerical':{'vev_gev':vev},
                    'base_state':{'ordinary':{'lambdaT3':{'real':lam.real,'imag':lam.imag}}}})
    write_json(root/'output/general_observable/lo_direct_one_loop/source_audits/component_support_audit.json',
               {'cases':cases,'errors':[]})
    write_json(root/'output/general_observable/lo_direct_one_loop/source_audits/matchete_generator_analysis.json',gen)
    write_json(root/'output/general_observable/lo_direct_one_loop/source_audits/cg_covariance_audit.json',
               {'status':'cg_generator_covariance_verified','models':cov})
    # Simulate the real frozen inputs' file-provenance contract, including
    # Windows-style relative paths. The production audit refuses mismatches.
    from Numerical.lo_certification.direct_uv.native_mixing_dictionary import digest
    audit_hash = digest(root/AUDIT_DIR/'component_support_audit.json')
    generator_hash = digest(root/AUDIT_DIR/'matchete_generator_analysis.json')
    for name,key,_ in models():
        path = root/'output/general_observable/lo_direct_one_loop/provisional_inputs'/f'{name}.json'
        record = read(path)
        record['provenance'] = {
            'cg_support_path': str(AUDIT_DIR/'component_support_audit.json').replace('/', '\\'),
            'generator_analysis_path': str(AUDIT_DIR/'matchete_generator_analysis.json').replace('/', '\\'),
            'cg_support_sha256': audit_hash,
            'generator_analysis_sha256': generator_hash,
            'benchmark_sha256': digest(root/'Numerical/benchmarks/fitted_reference'/f'{key}.json'),
        }
        write_json(path, record)
    return report,tensors


@pytest.mark.parametrize('cls',CLASSES)
def test_each_native_class_tensor(cls):
    report,tensors=seed_native()
    np.testing.assert_allclose(native_tensor(next(x for x in report['models'] if x['class']==cls),cls),tensors[cls])


def test_all_five_native_classes():
    report,_=seed_native()
    assert len(native_models(report))==5


def test_native_doublet_vev_orientation():
    report,ts=seed_native()
    assert all(abs(ts[cls][1,1]).max()>0 for cls in CLASSES)
    assert charge(2,0.5,1)==0


def test_nonsymmetric_higgs_rejected():
    r,_=seed_native()
    row=r['models'][0]
    flat=np.array(row['native_cg']['components_real']).reshape((2,2,1,3))
    flat[0,1,0,0]=3.
    row['native_cg']['components_real']=flat.reshape(-1).tolist()
    with pytest.raises(MixingAuditError,match='Higgs Bose'):
        native_models(r)


def test_nonrestricted_native_rejected():
    r,_=seed_native()
    r['restricted_environment']['effective'][FLAGS[0]]='1'
    with pytest.raises(MixingAuditError,match='nonrestricted'):
        native_models(r)


def test_wrong_orientation_rejected():
    r,_=seed_native()
    r['models'][0]['scalar_fields_in_nonhc_operator']='Bar[S1] S2'
    with pytest.raises(MixingAuditError,match='operator orientation'):
        native_models(r)


def test_wrong_native_tensor_dimensions_rejected():
    r,_=seed_native()
    r['models'][0]['native_cg']['dimensions']=[2,2,3,1]
    with pytest.raises(MixingAuditError,match='expected registered'):
        native_models(r)


def test_complex_conjugate_lower_left_action():
    t=np.array([[[0,0],[0,0]],[[0,0],[0,3]]],complex).reshape((2,2,1,2))
    lam=0.7+0.3j
    x=mixing_block(t,lam,2.)
    assert abs(x[1,0]+6*lam)<1e-13
    assert abs(x[1,0]+6*lam.conjugate())>1e-3


def test_higgs_vev_factor_two_negative_control():
    t=np.ones((2,2,1,1),complex)
    x=mixing_block(t,2+0j,5.)
    assert x[0,0]==-25
    assert x[0,0]!=-50


def test_complex_hessian_exact_cartesian():
    A=np.array([[3,1+2j],[1-2j,4]],complex)
    _h,r=finite_hessian(A)
    assert r<1e-13


def test_nonunitary_generator_rejected():
    data={'representations':[{'status':'generator_algebra_verified','dimension':2,
       'U_columns_are_weight_eigenvectors':[[1,1],[0,1]]}]}
    with pytest.raises(MixingAuditError,match='Nonunitary'):
        generator_rotations(data)


def test_weight_rotation_matches_explicit_sum():
    rs=np.random.default_rng(13)
    t=rs.normal(size=(2,2,1,3)) + 1j*rs.normal(size=(2,2,1,3))
    u2=np.array([[1,1],[-1,1]])/np.sqrt(2)
    u3=np.eye(3)
    got=rotated_tensor(t,{1:np.eye(1),2:u2,3:u3},1,3)
    expect=np.einsum('abcd,ai,bj,ck,dl->ijkl',t,u2,u2,np.eye(1),u3,optimize=False)
    np.testing.assert_allclose(got,expect)


def test_all16_end_to_end_synthetic():
    with tempfile.TemporaryDirectory() as temp:
        root=Path(temp)
        r,ts=seed_root(root)
        native_path=root/'export.json'
        write_json(native_path,r)
        answer=audit(root,native_path)
        assert answer['native_registered_mixing_classes_checked']==5
        assert answer['actual_frozen_inputs']['models_checked']==16
        assert answer['actual_frozen_inputs']['max_scalar_A_mixing_relative_residual']<1e-12
        assert answer['full_physical_source_to_Weyl_certified_models']==0


def test_cg_negative_control_all16():
    with tempfile.TemporaryDirectory() as temp:
        root=Path(temp)
        r,_=seed_root(root)
        cases_path=root/'output/general_observable/lo_direct_one_loop/source_audits/component_support_audit.json'
        x=read(cases_path)
        x['cases'][0]['vertices']['T3MixCG']['components'][0]['cg']='4'
        write_json(cases_path,x)
        with pytest.raises(MixingAuditError,match='newly registered Matchete mixing CG differs'):
            all_sixteen(root,native_models(r))


def test_frozen_mixing_negative_control():
    with tempfile.TemporaryDirectory() as temp:
        root=Path(temp)
        r,_=seed_root(root)
        p=root/'output/general_observable/lo_direct_one_loop/provisional_inputs/A_-4.json'
        x=read(p)
        x['scalar_A_gev2']['real'][3][0]+=10.
        write_json(p,x)
        with pytest.raises(MixingAuditError,match='scalar-A not Hermitian|mixing-to-scalar-A mismatch'):
            all_sixteen(root,native_models(r))


def test_stale_native_export_blocked():
    with tempfile.TemporaryDirectory() as temp:
        root=Path(temp)
        r,_=seed_root(root)
        path=root/'Lagrangian/interactions/T3Topology.wl'
        path.write_text(path.read_text()+' changed',encoding='utf-8')
        with pytest.raises(MixingAuditError,match='Source changed'):
            source_guard(root,r)


def test_absent_native_export_blocked():
    with tempfile.TemporaryDirectory() as temp:
        with pytest.raises(MixingAuditError,match='Missing required file'):
            audit(Path(temp),Path(temp)/'absent.json',check_frozen=False)


def test_frozen_scalar_B_negative_control():
    with tempfile.TemporaryDirectory() as temp:
        root=Path(temp)
        r,_=seed_root(root)
        p=root/'output/general_observable/lo_direct_one_loop/provisional_inputs/B_-1.json'
        x=read(p)
        x['scalar_B_gev2']['real'][0][0]=1e4
        write_json(p,x)
        with pytest.raises(MixingAuditError,match='nonzero scalar B'):
            all_sixteen(root,native_models(r))


def test_nonhc_half_factor_negative_control():
    r,_=seed_native()
    r['models'][0]['native_nonhc_LInt_coefficient']='+lambdaT3 /2 * T3MixCG'
    with pytest.raises(MixingAuditError,match='operator orientation'):
        native_models(r)


def test_current_committed_audit_path_used():
    with tempfile.TemporaryDirectory() as temp:
        root=Path(temp)
        r,_=seed_root(root)
        assert (root/AUDIT_DIR/'component_support_audit.json').is_file()
        answer=all_sixteen(root,native_models(r))
        assert answer['models_checked']==16
        assert answer['frozen_input_provenance_sha256_verified'] is True
        assert answer['component_support_audit_path'] == (AUDIT_DIR/'component_support_audit.json').as_posix()


def test_missing_canonical_audit_rejected_even_if_old_path_available():
    with tempfile.TemporaryDirectory() as temp:
        root=Path(temp)
        r,_=seed_root(root)
        p=root/AUDIT_DIR/'component_support_audit.json'
        old=root/'Reports/output/full/multiplicity/component_support_audit.json'
        old.parent.mkdir(parents=True, exist_ok=True)
        old.write_bytes(p.read_bytes())
        p.unlink()
        with pytest.raises(MixingAuditError,match='Missing required file'):
            all_sixteen(root,native_models(r))


def test_frozen_component_provenance_sha_mismatch_rejected():
    with tempfile.TemporaryDirectory() as temp:
        root=Path(temp)
        r,_=seed_root(root)
        p=root/'output/general_observable/lo_direct_one_loop/provisional_inputs/A_-4.json'
        obj=read(p)
        obj['provenance']['cg_support_sha256']='0'*64
        write_json(p,obj)
        with pytest.raises(MixingAuditError,match='SHA256 differs from frozen input provenance'):
            all_sixteen(root,native_models(r))


def test_frozen_component_provenance_wrong_directory_rejected():
    with tempfile.TemporaryDirectory() as temp:
        root=Path(temp)
        r,_=seed_root(root)
        p=root/'output/general_observable/lo_direct_one_loop/provisional_inputs/A_-4.json'
        obj=read(p)
        obj['provenance']['cg_support_path']='Reports\\output\\full\\multiplicity\\component_support_audit.json'
        write_json(p,obj)
        with pytest.raises(MixingAuditError,match='different source-audit directory'):
            all_sixteen(root,native_models(r))
