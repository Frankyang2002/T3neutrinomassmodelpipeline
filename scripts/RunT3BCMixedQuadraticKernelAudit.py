"""Read-only B/C mixed quadratic resummation + source-leading matching audit.

Intended repository path: scripts/RunT3BCMixedQuadraticKernelAudit.py
Invoked after RunT3BCRunningOneScaleTriangle. Writes new reports only.

1. Independent old-basis resolvent versus exact eigenbasis three-mass loop.
2. Actual seven native Matchete extra-quartic tensors rotated at running U.
3. Full-mass vs leading-1/MF C5 with the *existing* v22 F/Y/X weights.
4. Hard+soft MSbar loop identity across multiple scales, no Wilson matching.

Does NOT claim native full UV amplitude or complete EFT counterterms.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import re
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Numerical.experimental.T3BCMixedQuadraticKernel import (
    require, positive_mass_matrix, euclidean_matrix_kernel,
    signed_exact_master, strict_heavy_mass_master, hard_soft_master,
    invariant_quartic_contraction,
)

INPUT_DIR = ROOT/'output/studies/t3_bc_running_one_scale_triangle_v1'
OUTPUT_DIR = ROOT/'output/studies/t3_bc_mixed_quadratic_kernel_v2'
NATIVE_SOURCE = ROOT/'output/general_t3_stage22/matchete_tensors.json'
MODELS=('B_-1','C_-1')
SOURCE_V22=ROOT/'scripts/DirectVsSequentialOneLoop.experimental.py'
SOURCE_V18=ROOT/'scripts/EvaluateT3BCThreeFermionTriangle.experimental.py'
SOURCE_V16=ROOT/'scripts/AssembleT3BCFormalFiniteScaleEnvelope.experimental.py'
SOURCE_V17=ROOT/'scripts/CompareT3BCFiniteEnvelopeToTriangle.experimental.py'
FORMAT='t3-bc-running-mixed-quadratic-original-basis-loop-audit-v2'


def read(path:Path) -> dict:
    require(path.is_file(),f'missing source {path}')
    result=json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(result,dict),'source JSON must be object')
    return result


def sha(path:Path)->str:
    require(path.is_file(),f'missing hashed source {path}')
    return sha256(path.read_bytes()).hexdigest()


def complex_array(obj:dict, shape:tuple[int,...]) -> np.ndarray:
    require(isinstance(obj,dict) and set(obj)=={'real','imag'},'bad complex array schema')
    a=np.asarray(obj['real'],dtype=float)+1j*np.asarray(obj['imag'],dtype=float)
    require(a.shape==shape and np.isfinite(a).all(),f'invalid complex shape {shape}')
    return a


def load_result(model:str,folder:Path)->tuple[Path,dict,np.ndarray,np.ndarray,np.ndarray]:
    require(model in MODELS,'unsupported B/C model')
    tag=model.replace('-','m')
    path=folder/f'{tag}_running_one_scale.json'
    data=read(path)
    require(data.get('format')=='t3-bc-common-uv-scale-restricted-direct-matchete-triangle-v1' and
            data.get('status')=='RunningCommonScaleRestrictedTriangleEvaluatedNotPhysicalC5' and
            data.get('model')==model,'not the checked running-input B/C source')
    for k in ('physical_C5_certified','full_fermion_X_Y_sequential_RG_transport',
              'general_extra_mixed_mass_source_matching_certified',
              'native_matchete_rgbeta_action_dictionary_certified'):
        require(data.get(k) is False,f'{k} changed; review rather than assume')
    require(data.get('common_running_parameters_verified_same_UV_scale') is True,
            'not a common UV endpoint')
    check=data.get('running_input_direct_vs_sequential_parity',{})
    require(check.get('status')=='RunningInputRestrictedFullMassTriangleParityPassed'
            and float(check.get('max_relative_residual',1))<2e-8,
            'v1 independent restricted triangle parity missing')
    s=data['scalar_mass_basis']; A=complex_array(s['mass_squared_flavour_matrix_gev2'],(2,2))
    U=complex_array(s['unitary_flavour_rotation_U'],(2,2))
    eig,_=positive_mass_matrix(A)
    require(np.linalg.norm(U.conj().T@A@U-np.diag(eig))<1e-9,
            'source U does not diagonalize scalar A')
    require(np.allclose(np.sqrt(eig),[data['muX_gev'],data['muY_gev']],rtol=1e-10),
            'running mass/eigenvalue mismatch')
    mf=np.asarray(data['MF_takagi_masses_gev'],float)
    require(mf.shape==(3,) and np.all(np.diff(mf)>0), 'bad heavy mass inventory')
    return path,data,A,U,mf


def native_quartic_transport(model:str,U:np.ndarray,A:np.ndarray,native_file:Path)->dict:
    """Actual Matchete-native seven invariant tensors; imports audited source methods."""
    from Numerical.experimental.GeneralT3MixedScalarMassBasis import MixedDoubletMassBasis
    from scripts import AuditT3BCNativeQuarticMassBasis as native
    require(native_file.is_file(),'missing actual native Matchete quartic tensor inventory')
    basis=MixedDoubletMassBasis.build(model,float(A[0,0].real),float(A[1,1].real),complex(A[0,1]))
    require(np.linalg.norm(basis.U-U)<1e-8,'running eigenbasis convention unexpectedly changed')
    vertices=native.native_tensors(native.read_json(native_file),model)
    require(len(vertices)==7,'not seven source quartics')
    results={}; rotated={}
    for i,(name,(fields,tensor)) in enumerate(vertices.items()):
        pieces=native.rotate_tensor(tensor,fields,U)
        error=native.verify_polynomial(fields,tensor,pieces,basis,trials=16,seed=4100+i)
        results[name]={'fields':list(fields),'rotated_field_signatures':len(pieces),
                        'native_component_norm':float(np.linalg.norm(tensor)),
                        'max_relative_polynomial_and_hc_error':error}
        rotated[name]=pieces
    pair=native.shared_independent_block(rotated['lambdaT3ExtraQ110101I1'],
                                          rotated['lambdaT3ExtraQ110101I2'])
    # Zero-HH source quartics cannot create two H external fields via a single
    # quartic insertion in the unbroken theory; scalar mass insertions only
    # alter scalar propagators, not the number of external Higgs legs.
    holomorphic=['lambdaT3ExtraQ020020I1','lambdaT3ExtraQ020200I1']
    require(set(holomorphic).issubset(results),'two source LLHH extra quartics missing')
    return {'status':'NativeSevenQuarticsRotatedAtRunningEigenbasis',
            'native_matchete_tensors_sha256':sha(native_file),
            'quartic_count':len(vertices),
            'source_vertices':results,
            'two_independent_110101_contractions':pair,
            'two_extra_holomorphic_LLHH_quartics':holomorphic,
            'baseline_topology_LLHH_quartic':'lambdaT3',
            'no_new_HH_vertex_from_mixed_quadratic_insertions':True,
            'all_native_baseline_topology_yukawa_vertices_rotated':False,
            'physical_full_interaction_normalization_certified':False}


def source_module(path:Path,label:str):
    require(path.is_file(),f'missing source code {path}')
    spec=importlib.util.spec_from_file_location(label,path)
    require(spec is not None and spec.loader is not None,'source import unavailable')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def _assemble(masses:np.ndarray,y1:np.ndarray,y2:np.ndarray,
              coefficient_for_mass)->np.ndarray:
    result=np.zeros((3,3),complex)
    for i,m in enumerate(masses):
        a11,a12,a22=map(complex,coefficient_for_mass(float(m)))
        K=np.array([[2*a11,a12],[a12,2*a22]],complex)
        Y=np.column_stack((y1[:,i].conj(),y2[:,i].conj()))
        result+=Y@K@Y.T
    return result


def strict_full_mass_response(model:str,data:dict,U:np.ndarray,
                              masses:np.ndarray)->dict:
    """Rebuild source WXX,WXY,WYY to quantify the strict-1/MF truncation.

    Does not infer new operator coefficients: both full and strict amplitudes
    use v22's existing matched F/Y/X source polynomials and same running U.
    """
    import math
    tag=model.replace('-','m')
    v22=source_module(SOURCE_V22,'t3bc_v22_mixed_v2_'+tag)
    v18=source_module(SOURCE_V18,'t3bc_v18_mixed_v2_'+tag)
    v16=source_module(SOURCE_V16,'t3bc_v16_mixed_v2_'+tag)
    v17=source_module(SOURCE_V17,'t3bc_v17_mixed_v2_'+tag)
    weights=v22.sequential_weights(ROOT,tag,v16,v17,v18)
    quartics={k:complex(float(v['real']),float(v['imag']))
              for k,v in data['conditional_matchete_quartics'].items()}
    require(set(quartics)==set(v22.Q_NAMES),'unexpected restricted matched quartic axes')
    numeric=tuple(v18.polynomial_values(w,v22.substitution(U,quartics))
                  for w in weights)
    mx=float(data['muX_gev']);my=float(data['muY_gev'])
    y1=complex_array(data['running_heavy_eigenbasis_y1_3x3'],(3,3))
    y2=complex_array(data['running_heavy_eigenbasis_y2_3x3'],(3,3))
    def strict_coeff(mf:float):
        phi=np.array([strict_heavy_mass_master(mf*mf,mx*mx,mx*mx),
                      strict_heavy_mass_master(mf*mf,my*my,my*my),
                      strict_heavy_mass_master(mf*mf,mx*mx,my*my)])
        ks=sum(numeric[i]*phi[i] for i in range(3))/(16*math.pi**2*mf)
        return tuple(complex(z) for z in ks)
    def full_coeff(mf:float):
        return v22.expected_ordered(numeric,mf,mx,my,v17)
    full=_assemble(masses,y1,y2,full_coeff)
    strict=_assemble(masses,y1,y2,strict_coeff)
    source=complex_array(data['candidate_ordered_source_C5_3x3_gev_inverse_NOT_PHYSICAL'],(3,3))
    parity=float(np.linalg.norm(full-source)/max(np.linalg.norm(source),1e-40))
    require(parity<2e-8,f'v22 source weight parity failed at running U: {parity}')
    relative=float(np.linalg.norm(strict-full)/max(np.linalg.norm(full),1e-40))
    v=float(data['reference_EW_vev_gev'])
    mass_full=np.sort(np.linalg.svd(-v*v/2*full*1e9,compute_uv=False))[::-1]
    mass_strict=np.sort(np.linalg.svd(-v*v/2*strict*1e9,compute_uv=False))[::-1]
    return {'source_matchete_full_mass_parity_relative':parity,
            'leading_1_over_MF_restricted_C5_relative_error':relative,
            'full_mass_reference_singular_ev_NOT_PHYSICAL':mass_full.tolist(),
            'leading_1_over_MF_reference_singular_ev_NOT_PHYSICAL':mass_strict.tolist(),
            'full_mass_and_strict_use_same_yukawas_and_quartics':True,
            'v22_fixed_order_matchete_source_polynomials_reused':True,
            'full_renormalized_EFT_counterterms_computed':False}


def complex_majorana_spinor_stress(model:str,data:dict,source_path:Path) -> dict:
    """Synthetic *complex* MF demonstrates whether the archived all-real
    benchmark has any power to resolve the heavy-spinor conjugation choice.

    Synthetic phases change the UV input; this is NOT a physical parameter
    scan or an alternative prediction. The two conditional prescriptions use
    the identical source Matchete mass function and identical three masses.
    """
    from Numerical.experimental.T3BCRunningOneScaleTriangle import (
        takagi_majorana, assemble_ordered_c5,
    )
    tag=model.replace('-','m')
    v22=source_module(SOURCE_V22,'t3bc_v22_spinor_stress_'+tag)
    direct_path,_,_,_=v22.source_gates(ROOT,tag)
    require(sha(direct_path)==data['frozen_source_parity']['native_direct_source_sha256'],
            'native 480-term source changed since v1')
    compiled,denom,counts,_=v22.compile_source(
        direct_path.read_text(encoding='utf-8-sig'),model)
    require(len(compiled)==480, 'stress matching source missing 480 terms')
    mf=np.asarray(data['MF_takagi_masses_gev'],float)
    # Complex diagonal phases ensure the Takagi unitary is genuinely complex;
    # the real benchmark's y1@V vs y1@V* agreement is otherwise automatic.
    phases=np.array([.19,-.13,.27])
    synthetic_M=np.diag(mf*np.exp(1j*phases))
    masses,V,checks=takagi_majorana(synthetic_M)
    require(np.linalg.norm(V.imag)>0.05,'synthetic heavy rotation accidentally real')
    y1=complex_array(data['running_heavy_eigenbasis_y1_3x3'],(3,3))
    y2=complex_array(data['running_heavy_eigenbasis_y2_3x3'],(3,3))
    U=complex_array(data['scalar_mass_basis']['unitary_flavour_rotation_U'],(2,2))
    mx=float(data['muX_gev']);my=float(data['muY_gev'])
    quartics={k:complex(v['real'],v['imag'])
              for k,v in data['conditional_matchete_quartics'].items()}
    def ev(mf0,mx0,my0,U0,q0):
        return v22.direct_ordered(compiled,denom,mf0,mx0,my0,U0,q0)
    first,_=assemble_ordered_c5(masses,y1@V.conj(),y2@V,mx,my,U,quartics,ev)
    second,_=assemble_ordered_c5(masses,y1@V,y2@V,mx,my,U,quartics,ev)
    delta=float(np.linalg.norm(first-second)/max(np.linalg.norm(first),1e-40))
    require(np.isfinite(delta),'synthetic complex-MF response nonfinite')
    return {'status':'SyntheticComplexHeavyMassPhaseSensitivityCalculated',
            'synthetic_phase_radians':phases.tolist(),
            'same_three_takagi_masses_gev':masses.tolist(),
            'Takagi_checks':checks,
            'complex_heavy_rotation_imaginary_norm':float(np.linalg.norm(V.imag)),
            'alternative_y1_V_vs_Vconj_relative_C5_difference':delta,
            'real_benchmark_difference_original':float(data['alternative_same_V_spinor_rotation_relative_C5_difference']),
            'synthetic_only_no_UV_running_recomputed':True,
            'correct_physical_heavy_spinor_convention_certified':False}


def run_model(model:str,source_dir:Path,native_file:Path,output_dir:Path,
              *,skip_native:bool=False,skip_source_weights:bool=False,
              skip_spinor_stress:bool=False)->dict:
    input_path,data,A,U,mf=load_result(model,source_dir)
    eigen=np.linalg.eigvalsh(A)
    kernels=[];master=[]
    for idx,mass in enumerate(mf):
        test=euclidean_matrix_kernel(A,float(mass)**2)
        covariance=invariant_quartic_contraction(A,float(mass)**2,seed=320+idx)
        require(covariance<3e-12,'arbitrary quartic kernel not basis covariant')
        test['source_basis_covariant_arbitrary_two_scalar_vertex_residual']=covariance
        kernels.append(test)
        for i,j in ((0,0),(0,1),(1,1)):
            full=signed_exact_master(float(mass)**2,eigen[i],eigen[j])
            strict=strict_heavy_mass_master(float(mass)**2,eigen[i],eigen[j])
            scales={label:hard_soft_master(float(mass)**2,eigen[i],eigen[j],s*s)
                    for label,s in {'MF':float(mass),'Y':float(np.sqrt(eigen[1])),
                                    'X':float(np.sqrt(eigen[0]))}.items()}
            err=abs(full-strict)/max(abs(full),1e-30)
            master.append({'heavy_column':idx,'scalar_pair':{(0,0):'XX',(0,1):'XY',(1,1):'YY'}[(i,j)],
                           'full_mass_signed_phi':full,'strict_1_over_MF_phi':strict,
                           'relative_mass_function_correction':err,'hard_soft_MSbar':scales})
    native_status=(native_quartic_transport(model,U,A,native_file)
                   if not skip_native else {'status':'NativeSourceAuditSkippedByExplicitOption'})
    matching=(strict_full_mass_response(model,data,U,mf)
              if not skip_source_weights else {'status':'SourceWeightAuditSkippedByExplicitOption'})
    spinor=(complex_majorana_spinor_stress(model,data,input_path)
            if not skip_spinor_stress else {'status':'ComplexMFStressSkippedByExplicitOption'})
    result={'format':FORMAT,'status':('FullRestrictedSourceKernelChecksPassed' if not (skip_native or skip_source_weights or skip_spinor_stress)
                      else 'KernelOnlyPartialExplicitSkips'),
            'model':model,'input_running_one_scale_sha256':sha(input_path),
            'native_quartic_source':native_status,
            'running_scalar_A_gev2':{'real':A.real.tolist(),'imag':A.imag.tolist()},
            'heavy_mass_momentum_kernels':kernels,
            'hard_soft_mass_function_checks':master,
            'full_mass_vs_strict_source_coefficient':matching,
            'complex_majorana_spinor_convention_stress':spinor,
            'max_original_basis_eigenbasis_relative_error':max(x['original_basis_vs_spectral_relative'] for x in kernels),
            'largest_zero_insertion_kernel_error':max(x['dyson_relative_errors']['0'] for x in kernels),
            'max_1_over_MF_mass_function_relative_error':max(x['relative_mass_function_correction'] for x in master),
            'physical_C5_certified':False,'full_sequential_EFT_matching_performed':False,
            'native_action_interpackage_dictionary_certified':False,
            'new_matchete_match_with_independent_mixed_quadratic_performed':False,
            'scope':'Mathematical mixed scalar propagator resummation and running-U native extra quartic rotation; restricted archived y^2 lambda source weights; not a new physical Wilson matching'}
    output_dir.mkdir(parents=True,exist_ok=True)
    path=output_dir/(model.replace('-','m')+'_mixed_kernel.json')
    require(not path.exists(),'refusing to overwrite existing report: '+str(path))
    tmp=path.with_suffix('.json.tmp')
    tmp.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    tmp.replace(path)
    return result


def main()->int:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input-dir',type=Path,default=INPUT_DIR)
    ap.add_argument('--output-id',default='t3_bc_mixed_quadratic_kernel_v2')
    ap.add_argument('--native-file',type=Path,default=NATIVE_SOURCE)
    ap.add_argument('--model',choices=MODELS)
    ap.add_argument('--skip-native',action='store_true')
    ap.add_argument('--skip-source-weights',action='store_true')
    ap.add_argument('--skip-spinor-stress',action='store_true')
    ns=ap.parse_args()
    require(bool(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,60}',ns.output_id)),
            'invalid output ID')
    out=ROOT/'output/studies'/ns.output_id
    require(not out.exists(),'output already exists; use new --output-id')
    rows=[];failures=[]
    for m in ((ns.model,) if ns.model else MODELS):
        try:
            v=run_model(m,ns.input_dir,ns.native_file,out,skip_native=ns.skip_native,
                        skip_source_weights=ns.skip_source_weights,
                        skip_spinor_stress=ns.skip_spinor_stress)
            rows.append({'model':m,'status':v['status'],
                         'mass_resummation_error_max':v['max_original_basis_eigenbasis_relative_error'],
                         'omitted_mixing_kernel_error_max':v['largest_zero_insertion_kernel_error'],
                         'strict_mass_function_difference_max':v['max_1_over_MF_mass_function_relative_error'],
                         'native_quartics':v['native_quartic_source']['status'],
                         'C5_strict_vs_full_source':v['full_mass_vs_strict_source_coefficient'].get('leading_1_over_MF_restricted_C5_relative_error'),
                         'complex_MF_spinor_choice_relative_C5_change':v['complex_majorana_spinor_convention_stress'].get('alternative_y1_V_vs_Vconj_relative_C5_difference')})
            print(m, 'PASS', 'kernel_residual',v['max_original_basis_eigenbasis_relative_error'],
                  'zero_insertion',v['largest_zero_insertion_kernel_error'],flush=True)
        except Exception as exc:
            failures.append({'model':m,'type':type(exc).__name__,'error':str(exc)})
            print(m,'BLOCKED',exc,flush=True)
    out.mkdir(parents=True,exist_ok=True)
    summary={'format':FORMAT,'status':('CompleteRestrictedKernelAudits' if not failures and all(
                 x['status']=='FullRestrictedSourceKernelChecksPassed' for x in rows) else
                 'CompletedSubsetByExplicitSkip' if not failures else 'PartialBlocked'),
             'models':rows,'errors':failures,'physical_C5_certified':False,
             'complete_sequential_matching_certified':False}
    (out/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('SUMMARY',summary['status'],out/'summary.json',flush=True)
    return 0 if not failures else 3

if __name__=='__main__':raise SystemExit(main())
