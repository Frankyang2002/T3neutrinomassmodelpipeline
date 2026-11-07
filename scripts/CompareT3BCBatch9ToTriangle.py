#!/usr/bin/env python3
"""Batch 10 v3: compare B/C Batch 9 with archived independent direct triangle.

Local, read-only physics inputs. No Matchete, matching, or RGE calls.
Requires current GitHub-committed B/C direct all-heavy Matchete 480-term
ordered A5 sources, not the missing historical v14--v17 source audits.

The native direct all-heavy amplitude is a previously established (v22)
three-propagator full-mass triangle reference, but its native source alone
is NOT an independent Feynman-rule normalization certificate.

Compare full complex symmetric C5 matrices for identical running muF inputs.
Report both the *full mass* direct coefficient (different MF power counting)
and an independently extrapolated leading-1/MF coefficient, retaining exact
MX/MY dependence. The latter is extracted from the direct source's asymptotic
logarithm while explicitly removing the known O(log(MF)/MF^3) correction.
No coefficient is fitted to Batch 9. Separately compare its log slope with
Batch 9's independently projected F-first pole.

These are conditional source/kernel checks, not a renormalized physical C5.
The mandatory EFT output order remains UV -> all F -> EFT1(SM+X+Y) -> all
X,Y -> SMEFT. The direct all-heavy expression is a *cross-check* only.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import importlib.util
import json
import math
from pathlib import Path
import re
import sys
from typing import Any, Callable

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
STEMS = {'B_-1':'B_m1','C_-1':'C_m1'}
Q_NAMES = {'lambdaT3','lambdaT3ExtraQ020020I1','lambdaT3ExtraQ020200I1'}
BATCH9_FORMAT = 't3-mixed-f-inherited-scalar-common-pole-v1'
BATCH9_STATUS = 'CommonCouplingFInheritedAndScalarRawPolesCancelled_Conditional'
DIRECT_FORMAT = 't3-matchete-native-heavy-mass-eigenbasis-v1'
DIRECT_STATUS = 'OneLoopDirectSymbolicWeinbergProbePassedNumericalC5Pending'
OUTPUT_BASE = Path('output/general_sequential_pilot/mixed_triangle_compare_batch10')

class ComparisonError(ValueError):
    pass

def require(ok: bool, message: str) -> None:
    if not ok:
        raise ComparisonError(message)

def digest(path: Path) -> str:
    require(path.is_file(), 'Required source missing: '+str(path))
    return sha256(path.read_bytes()).hexdigest()

def read_json(path: Path) -> dict:
    require(path.is_file(), 'Required JSON missing: '+str(path))
    data=json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(data,dict), 'Expected JSON object: '+str(path))
    return data

def encode_matrix(x: np.ndarray, shape=(3,3)) -> dict:
    x=np.asarray(x,dtype=complex)
    require(x.shape==shape and np.isfinite(x).all(), 'Expected finite complex matrix of shape '+str(shape))
    return {'real':x.real.tolist(),'imag':x.imag.tolist()}

def decode_matrix(x: Any, name: str, shape=(3,3)) -> np.ndarray:
    require(isinstance(x,dict) and 'real' in x and 'imag' in x,
            name+': missing real and imag')
    z=np.asarray(x['real'],float)+1j*np.asarray(x['imag'],float)
    require(z.shape==shape and np.isfinite(z).all(), name+': invalid complex matrix')
    return z

def relative_frobenius(a: np.ndarray, b: np.ndarray) -> float:
    a=np.asarray(a,complex);b=np.asarray(b,complex)
    return float(np.linalg.norm(a-b)/max(np.linalg.norm(a),np.linalg.norm(b),1e-40))

def symmetry_gate(z: np.ndarray, label: str, tolerance: float=1e-11) -> float:
    require(z.shape==(3,3),label+': expected 3x3')
    v=relative_frobenius(z,z.T)
    require(v<tolerance,label+': not symmetric under ordinary transpose')
    return v

def source_record(root: Path,path: Path) -> dict:
    try:
        p=path.relative_to(root)
    except ValueError:
        p=path
    return {'path':str(p).replace('\\','/'), 'sha256':digest(path)}

def import_source(path: Path, name: str):
    require(path.is_file(),'Missing existing source script: '+str(path))
    spec=importlib.util.spec_from_file_location(name,path)
    require(spec is not None and spec.loader is not None,'Cannot import '+str(path))
    mod=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

def identical(a: np.ndarray,b: np.ndarray,label: str,tolerance: float=5e-13):
    require(relative_frobenius(np.asarray(a),np.asarray(b))<=tolerance,
            'Saved numerical input mismatch: '+label)

def validate_direct_report(report: dict,model: str,source: Path) -> None:
    stem=STEMS[model]
    require(report.get('format')==DIRECT_FORMAT and
            report.get('model')==model and
            report.get('status')==DIRECT_STATUS and
            report.get('one_loop_match_stage_completed') is True and
            report.get('matching_approach') is not None and
            report.get('weinberg_extraction_status')=='Success' and
            report.get('weinberg_term_count')==236 and
            report.get('weinberg_holomorphic_term_count')==118 and
            report.get('weinberg_conjugate_term_count')==118 and
            report.get('ordered_coefficient_is_canonical_flavour_symmetric_c5') is False and
            report.get('all_17_quartics_present') is True and
            report.get('both_yukawa_interactions_present') is True and
            report.get('heavy_mass_basis_detailed_check_passed') is True and
            report.get('production_sequential_matching_performed') is False and
            report.get('numerical_c5_computed') is False and
            report.get('weinberg_ordered_coefficient_inputform_file')==source.name,
            model+': direct Matchete report schema/source certificate changed')

def assemble_direct(masses: np.ndarray,inp: dict,
                    coefficients: Callable[...,tuple[complex,complex,complex]]) -> np.ndarray:
    """Use the archived v22 ordered A5 -> C5 map without extra factors.

    Each ordered source coefficient already includes 1/(16 pi^2 MF).
    C5 = sum_r L_r [[2a11,a12],[a12,2a22]] L_r^T, L_r=(y1*,y2*).
    """
    mx,my=map(float,inp['ms'])
    require(0<mx<my and np.min(masses)>my,'Need MF > MY > MX > 0')
    out=np.zeros((3,3),complex)
    for r,m in enumerate(masses):
        a=np.asarray(coefficients(float(m),mx,my,inp['U'],inp['couplings']),complex)
        require(a.shape==(3,) and np.isfinite(a).all(),
                'Invalid ordered native direct coefficients')
        k=np.array([[2*a[0],a[1]],[a[1],2*a[2]]],complex)
        L=np.column_stack((inp['y1'][:,r].conj(),inp['y2'][:,r].conj()))
        out+=L@k@L.T
    symmetry_gate(out,'Direct all-heavy C5')
    return out

def fit_log_plus_power(factors: tuple[float, float, float, float],
                       evaluations: dict[float, np.ndarray]) -> tuple[np.ndarray,np.ndarray]:
    """Four-point extraction of leading dimension-five constant and log.

    For t=base*u, write t C5(t*MF) = A + B log(u)
         + [D + E log(u)]/u**2 + O(u**-4 log(u)/base**4).
    A = a + b*log(base), B = b. The fitted D,E remove the entire
    one-loop relative MF^-2 correction, including its mass logarithm.

    This is a linear solve fixed by heavy-mass power counting, not a fit
    to the Batch 9 target matrix. A and B have dimensions GeV^-1.
    """
    require(len(factors)==4 and all(t>0 and math.isfinite(t) for t in factors),
            'Four positive, finite rescaling factors required')
    base=float(factors[0])
    u=np.asarray(factors,float)/base
    x=np.log(u)
    design=np.column_stack((np.ones(4),x,1/u**2,x/u**2))
    require(np.linalg.cond(design)<1000,
            'Asymptotic fit design ill-conditioned')
    matrices=[np.asarray(evaluations[t],complex) for t in factors]
    require(all(a.shape==(3,3) and np.isfinite(a).all() for a in matrices),
            'Asymptotic samples must be finite 3x3 complex matrices')
    v=np.stack(matrices).reshape((4,9))
    A,B,D,E=np.linalg.solve(design,v).reshape((4,3,3))
    # Undo reference to log(u) to recover the t-independent strict coefficient.
    a=A-B*math.log(base)
    b=B
    symmetry_gate(a,'Asymptotic leading coefficient')
    symmetry_gate(b,'Asymptotic logarithmic slope')
    return a,b


def leading_from_direct(inp: dict, coefficients: Callable,base:float=128.) -> dict:
    """Extract strict dimension-five C5 without knowing Batch 9 C5.

    Fixed 4-point linear extraction analytically removes the universal
    O(log(t)/t**2) correction of a one-loop three-mass scalar triangle.
    Compare two independent overlapping windows, requiring convergence.
    The next omitted terms scale as O(log(t)/t**4), not O(log(t)/t**2).
    """
    require(32<=base<=4096 and math.isfinite(base),'Bad asymptotic scale')
    mf=np.asarray(inp['mf'],float)
    factors=tuple(float(base*2**j) for j in range(5))
    evaluations={}
    for factor in factors:
        evaluations[factor]=factor*assemble_direct(factor*mf,inp,coefficients)
    a1,b1=fit_log_plus_power(factors[:4],evaluations)
    a2,b2=fit_log_plus_power(factors[1:],evaluations)
    return {'leading':a2,'slope':b2,
            'asymptotic_matrix_convergence_relative':relative_frobenius(a1,a2),
            'asymptotic_slope_convergence_relative':relative_frobenius(b1,b2),
            'samples':list(factors),
            'early_intercept':a1,'early_slope':b1,
            'fit_ansatz':'t*C5(t*MF)=a+b*ln(t)+(c+d*ln(t))/t^2+O(ln(t)/t^4)',
            'subleading_log_over_t2_removed':True}

def input_gates(root: Path,model: str):
    stem=STEMS[model]
    p9=root/'output/general_sequential_pilot/mixed_inherited_pole_batch9/mixed_inherited_11'/stem/'result.json'
    p8=root/'output/general_sequential_pilot/mixed_fixed_order_batch6/mixed_fixed_08'/stem/'result.json'
    pr=root/'output/general_observable/mixed_adapter_v3'/stem/'running.json'
    pp=(root/'output/general_observable/mixed_f_synchronised_v1'/stem/'matchete'/
        f'{stem}_final_parameters.json')
    ph=(root/'output/general_sequential_pilot/mixed_scalar_hard_batch6/mixed_hard_08'/stem/
        'hard_c5_candidate.txt')
    pf=(root/'output/general_observable/lo_certification_spinor'/
        'native_hard_e515d82d118284ba_native01'/stem/
        f'{stem}_native_one_loop_raw_hard_llhh_sector.inputform.txt')
    pw=root/'Lagrangian/RunMatching.wl'
    pd=root/'output/general_t3_mass_eigen_matchete/one_loop_probe'/(
        f'{stem}_final_one_loop_probe_weinberg_ordered_coefficient.inputform.txt')
    pj=root/'output/general_t3_mass_eigen_matchete/one_loop_probe'/(
        f'{stem}_final_one_loop_probe.json')
    pathmap={'batch9':p9,'batch8':p8,'running':pr,'parameters':pp,
             'scalar_hard':ph,'F_first_native_raw':pf,'weinberg_definition':pw,
             'direct_480_terms':pd,'native_direct_report':pj}
    for p in pathmap.values():digest(p)
    b9,b8=read_json(p9),read_json(p8)
    require(b9.get('format')==BATCH9_FORMAT and b9.get('model')==model and
            b9.get('status')==BATCH9_STATUS and
            b9.get('physical_C5_certified') is False,
            'Missing passed Batch9 v11 evidence for '+model)
    require(b8.get('model')==model and b8.get('quartic_map_is_conditional') is True and
            b8.get('physical_C5_certified') is False,
            model+': unexpected Batch 8 numerical evidence')
    require(b9['source_batch8']['fixed_result_sha256']==digest(p8) and
            b9['source_batch8']['hard_scalar_sha256']==digest(ph)==b8['hard_scalar_source_sha256'] and
            b9['source_archived_F_first_one_loop']['sha256']==digest(pf) and
            b9['scalar_ordered_to_symmetric_C5_convention']['source_sha256']==digest(pw) and
            b9['scalar_ordered_to_symmetric_C5_convention']
              ['factor_from_saved_batch8_Asym_to_C5']==2,
            model+': source fingerprint or fixed A5/C5 normalization changed')
    require(b8['source_fingerprints']['source_fingerprints']['running']==digest(pr) and
            b8['source_fingerprints']['source_fingerprints']['parameters']==digest(pp),
            model+': Batch 8 frozen running or parameters file changed')
    require(b9.get('source_pole_check',{}).get(
                'common_pole_cancellation_verified_in_raw_projected_sectors') is True and
            b9.get('common_scale_finite_matrix_relative_change',1)>-1 and
            b9.get('common_scale_finite_matrix_relative_change',1)<2e-8,
            model+': inherited raw-pole/common-scale checks not passed')
    validate_direct_report(read_json(pj),model,pd)
    # Reuse the project's numerical input decoder; never re-evaluate Match.
    from Numerical.running.T3MixedBatch6Numerics import numerical_inputs
    inp=numerical_inputs(model,read_json(pr),read_json(pp),conditional_quartic_map=True)
    require(set(inp['couplings'])==Q_NAMES,'Native quartic set changed')
    identical(inp['mf'],np.asarray(b8['fermion_mass_eigenvalues_gev'],float),'MF')
    identical(inp['ms'],np.asarray(b8['scalar_mass_eigenvalues_gev'],float),'MX/MY')
    identical(inp['U'],decode_matrix(b8['scalar_U'],'U',(2,2)),'U')
    require(b8['common_couplings_at']=='muF' and b9['source_couplings_frozen_at']=='muF' and
            float(b8['muS_gev'])==float(b9['comparison_scale_gev']),
            model+': frozen couplings/comparison-scale mismatch')
    c8=decode_matrix(b8['hard_finite_c5_gev_inverse'],'Batch 8 finite A5')
    p8pole=decode_matrix(b8['hard_pole_residue_c5_gev_inverse'],'Batch 8 pole A5')
    identical(2*c8,decode_matrix(b9['conditional_raw_hard_finite_c5_gev_inverse'],
             'Batch 9 scalar finite C5'),'A5->C5 finite')
    identical(2*p8pole,decode_matrix(b9['hard_scalar_pole_residue_c5_gev_inverse'],
             'Batch 9 scalar pole C5'),'A5->C5 pole')
    c5=decode_matrix(b9['conditional_common_scale_sum_c5_gev_inverse'],'Batch9 combined')
    c9f=decode_matrix(b9['conditional_raw_F_inherited_finite_plus_log_c5_gev_inverse'],
                      'Batch 9 inherited F finite')
    identical(c5,2*c8+c9f,'Batch 9 matrix decomposition')
    symmetry_gate(c5,'Batch 9 C5')
    f_pole=decode_matrix(b9['archived_F_inherited_pole_residue_c5_gev_inverse'],'F-first pole')
    symmetry_gate(f_pole,'F-first pole')
    return inp,b9,c5,f_pole,pathmap

def compare_model(root:Path, model:str, v22, tolerance:float, sample_base:float) -> dict:
    inp,b9,batch9,f_pole,paths=input_gates(root,model)
    native=paths['direct_480_terms'].read_text(encoding='utf-8-sig')
    compiled,denom,shapes,quartics=v22.compile_source(native,model)
    require(len(compiled)==480 and set(quartics)==Q_NAMES and
            shapes=={'y1_i*y1_j':144,'y1_i*y2_j':192,'y2_i*y2_j':144},
            'Invalid native direct all-heavy source inventory')
    def coefficient(mf,mx,my,U,Qs):
        return v22.direct_ordered(compiled,denom,mf,mx,my,U,Qs)
    full=assemble_direct(np.asarray(inp['mf'],float),inp,coefficient)
    asym=leading_from_direct(inp,coefficient,sample_base)
    leading=asym['leading']
    # The triangle large-MF scale-log coefficient is -2 P_F, with
    # PF from independently archived F-first native one-loop LLHH.
    # Check this identity; it is *not* used in the asymptotic extraction.
    slope_expected=-2*f_pole
    slope_residual=relative_frobenius(asym['slope'],slope_expected)
    full_residual=relative_frobenius(batch9,full)
    leading_residual=relative_frobenius(batch9,leading)
    truncation_residual=relative_frobenius(full,leading)
    convergence=asym['asymptotic_matrix_convergence_relative']
    slope_convergence=asym['asymptotic_slope_convergence_relative']
    pass_gate=all(x<=tolerance for x in (leading_residual,convergence,
                                         slope_convergence,slope_residual))
    hashes={name:source_record(root,path) for name,path in paths.items()}
    extra_script=root/'scripts/DirectVsSequentialOneLoop.experimental.py'
    hashes['direct_parser']=source_record(root,extra_script)
    hashes['batch10_comparator']=source_record(root,root/'scripts/CompareT3BCBatch9ToTriangle.py')
    return {
        'model':model,
        'status':('ConditionalLeadingTriangleMatrixAgreement' if pass_gate else
                  'IndependentTriangleFiniteOrLogMatrixDisagreement_NotFrozen'),
        'comparison_is_direct_matchete_all_heavy_reference':True,
        'historical_direct_vs_triangle_v22_audit_reverified_this_run':False,
        'independent_Feynman_rule_vertex_sign_certified':False,
        'correct_EFT_sequence':'UV -> all F -> EFT1(SM+X+Y) -> all X,Y -> SMEFT',
        'direct_full_mass_and_batch9_have_different_MF_power_counting':True,
        'relative_frobenius':{
            'batch9_vs_direct_full_mass':full_residual,
            'batch9_vs_direct_leading_1_over_MF':leading_residual,
            'direct_full_vs_direct_leading_power_correction':truncation_residual,
            'independent_direct_large_MF_log_slope_vs_minus_2_F_pole':slope_residual,
            'leading_1_over_MF_extrapolation_matrix_convergence':convergence,
            'leading_1_over_MF_extrapolation_log_slope_convergence':slope_convergence},
        'matrices_c5_gev_inverse':{
            'batch9':encode_matrix(batch9),'direct_full_mass':encode_matrix(full),
            'direct_leading_1_over_MF':encode_matrix(leading),
            'batch9_minus_direct_full':encode_matrix(batch9-full),
            'batch9_minus_direct_leading':encode_matrix(batch9-leading),
            'direct_log_slope_per_log_MF':encode_matrix(asym['slope']),
            'minus_two_archived_F_pole':encode_matrix(slope_expected)},
        'asymptotic_rescaling_MF_factors':asym['samples'],
        'asymptotic_extraction_ansatz':asym['fit_ansatz'],
        'subleading_log_over_t2_removed':asym['subleading_log_over_t2_removed'],
        'independent_direct_source_terms':len(compiled),
        'independent_direct_yukawa_shape_counts':shapes,
        'independent_direct_quartic_counts':quartics,
        'inputs':{
            'all_couplings_evaluated_at':'muF',
            'muF_gev':float(inp['muF']),
            'muS_gev':float(b9['comparison_scale_gev']),
            'MF_gev':np.asarray(inp['mf'],float).tolist(),
            'MX_MY_gev':np.asarray(inp['ms'],float).tolist(),
            'U':encode_matrix(inp['U'],(2,2)),
            'y1':encode_matrix(inp['y1']),'y2':encode_matrix(inp['y2']),
            'native_quartics':{k:{'real':float(complex(v).real),'imag':float(complex(v).imag)}
                               for k,v in sorted(inp['couplings'].items())},
            'quartic_map_certified':False,
            'C5_normalization':'half-normalized symmetric C5=A5+A5^T; no fitted factor',
            'heavy_mass_lift_native_three_MF_verified':False},
        'source_hashes':hashes,
        'conditional_batch9_singular_values_ev':b9['conditional_common_scale_sum_singular_values_ev'],
        'conditional_direct_leading_singular_values_ev':
            np.sort(np.linalg.svd(leading,compute_uv=False)*246.22**2/2*1e9).tolist(),
        'passed_conditional_finite_and_pole_source_comparison':pass_gate,
        'complete_physical_C5':None,'physical_C5_certified':False,
        'UV_EFT_finite_field_EOM_counterterm_certified':False,
        'EFT1_running_performed':False,
    }

def self_test(root: Path) -> None:
    from math import log
    # Test asymptotic log extraction on an *independent* complex-symmetric
    # matrix function with a known 1/M, log(M)/M and 1/M^3 correction.
    a=np.array([[1+3j,2-1j,0],[2-1j,-4j,1],[0,1,-2]],complex)
    b=np.array([[2j,-1,0],[-1,2,0.5j],[0,0.5j,1j]],complex)
    c=np.array([[2,0,0],[0,1j,0],[0,0,-3]],complex)
    D=lambda t: a+b*log(t)+c*log(t)/t**2
    values={float(t):D(float(t)) for t in (128,256,512,1024,2048)}
    x,s=fit_log_plus_power((128.,256.,512.,1024.),values)
    x2,s2=fit_log_plus_power((256.,512.,1024.,2048.),values)
    require(relative_frobenius(x2,a)<5e-12,'Asymptotic log intercept test failed')
    require(relative_frobenius(s2,b)<5e-12,'Asymptotic log slope test failed')
    require(relative_frobenius(x,x2)<5e-12 and
            relative_frobenius(s,s2)<5e-12,
            'Four-point extrapolation windows inconsistent')
    symmetry_gate(a,'test complex Majorana matrix')
    require(relative_frobenius(a,2*a)>.4,'Incorrect factor-of-two accepted')
    # Run existing v22 safe Matchete parser and v17 triangle/quadrature tests;
    # these do NOT require any missing v14--v17 audit JSON files.
    v22=import_source(root/'scripts/DirectVsSequentialOneLoop.experimental.py',
                      't3_batch10_v3_direct_selftest')
    v22.self_test()
    v17=import_source(root/'scripts/CompareT3BCFiniteEnvelopeToTriangle.experimental.py',
                      't3_batch10_v3_master_selftest')
    v17.check_three_mass_integral()
    print('PASS: Batch 10 v3 four-point power/log extraction, direct parser, independent triangle integral')

def main(argv: list[str]|None=None) -> int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=ROOT)
    p.add_argument('--run-id',default='mixed_triangle_direct_03')
    p.add_argument('--relative-tolerance',type=float,default=2e-8)
    p.add_argument('--asymptotic-MF-factor',type=float,default=128.)
    p.add_argument('--self-test',action='store_true')
    a=p.parse_args(argv)
    root=a.root.resolve()
    if str(root) not in sys.path:sys.path.insert(0,str(root))
    if a.self_test:
        self_test(root)
        return 0
    require(re.fullmatch('[A-Za-z0-9_-]+',a.run_id) is not None,'Invalid run ID')
    require(math.isfinite(a.relative_tolerance) and 0<a.relative_tolerance<0.01,
            'Invalid relative matrix tolerance')
    out=root/OUTPUT_BASE/a.run_id
    require(not out.exists(),'Refusing to overwrite prior run: '+str(out))
    v22=import_source(root/'scripts/DirectVsSequentialOneLoop.experimental.py',
                      't3_batch10_v3_native_direct')
    models={}; errors={}
    for model in STEMS:
        try:
            row=compare_model(root,model,v22,a.relative_tolerance,a.asymptotic_MF_factor)
            models[model]=row
            r=row['relative_frobenius']
            print(f'{model}: {row["status"]}; leading={r["batch9_vs_direct_leading_1_over_MF"]:.4e}; '
                  f'full={r["batch9_vs_direct_full_mass"]:.4e}; '
                  f'log_slope={r["independent_direct_large_MF_log_slope_vs_minus_2_F_pole"]:.4e}; '
                  f'asym_conv={r["leading_1_over_MF_extrapolation_matrix_convergence"]:.4e}',
                  flush=True)
        except Exception as exc:
            errors[model]=f'{type(exc).__name__}: {exc}'
            print(model+': BLOCKED: '+errors[model],file=sys.stderr,flush=True)
    passed=(len(models)==2 and all(x['passed_conditional_finite_and_pole_source_comparison']
                                   for x in models.values()))
    result={
        'format':'t3-bc-batch9-direct-all-heavy-triangle-comparison-v3',
        'status':'ConditionalBCTriangleKernelBenchmarkFrozen' if passed else
                 'UncertifiedOrSourceMatrixDisagreement',
        'branch_reference':'GitHub main source read-only; local SHA256s recorded per input',
        'comparison_scope':'native all-heavy direct 480-term one-loop reference (v22 previously related to triangle; v22 audits unavailable), exact same Batch9 inputs; leading-1/MF asymptotic log projected independently',
        'relative_tolerance':a.relative_tolerance,
        'normalization_factors_fitted':False,
        'historical_v14_v15_v16_v17_audits_required':False,
        'historical_frozen_v18_yukawas_used':False,
        'ordinary_models_rerun':False,
        'fermion_matching_rerun':False,
        'scalar_matching_rerun':False,
        'EFT_order':'UV -> all F -> EFT1(SM+X+Y) -> all X,Y -> SMEFT',
        'physical_C5':None,'physical_C5_certified':False,
        'conditional_rough_benchmark_frozen':passed,
        'models':models,'errors':errors}
    out.mkdir(parents=True,exist_ok=False)
    (out/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('Evidence:',str(out/'summary.json'),flush=True)
    return 0 if passed else 2

if __name__=='__main__':
    try:
        raise SystemExit(main())
    except (ComparisonError,ValueError,KeyError,TypeError,OSError,ImportError) as exc:
        print('BLOCKED: '+str(exc),file=sys.stderr)
        raise SystemExit(2)
