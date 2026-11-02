"""Batch 3: final conditional T3 release and two provisional direct-UV audits.

This reads the validated Batch 2 release, reruns the EXISTING direct_one_loop
core on frozen source-mapped UV inputs for A_-4 and B_-1, and records both.
It never adjusts an EFT coefficient or fits an unknown normalization.

UV -> integrate ALL BSM fermions -> EFT1(S1,S2) -> integrate BOTH scalars
at ONE threshold -> SMEFT. The existing sequential response is unchanged.

The direct UV inputs are frozen at 10^7 GeV and DO NOT define the same
renormalized parameter point as the illustrative EFT1 threshold at 10^3 GeV.
Do not interpret their matrix difference as an EFT truncation error.
"""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
import math
from pathlib import Path
import shutil
from typing import Any, Mapping

import numpy as np

from Numerical.direct_one_loop import core as direct_core

ROOT = Path(__file__).resolve().parents[1]
SAFE = ROOT / 'output/general_observable/lo_joint_scalar_matching'
INPUT2 = SAFE / 'conditional_all16_batch2_v1'
UV_INPUT = ROOT / 'input/lo_direct_one_loop'
UV_PREVIOUS = ROOT / 'output/general_observable/lo_direct_one_loop/provisional'
OUTPUT = SAFE / 'final_conditional_all16_batch3_v1'

MODELS = ('A_-4', 'A_-2', 'A_+0', 'B_-3', 'B_-1', 'B_+1',
          'C_-3', 'C_-1', 'C_+1', 'D_-2', 'D_+0', 'D_+2',
          'E_-4', 'E_-2', 'E_+0', 'E_+2')
TARGETS = ('A_-4', 'B_-1')
BATCH2_FORMAT = 't3-conditional-all16-consolidated-batch2-v1'
SUCCESS = 'conditional_C5_per_unknown_NGamma_exploratory'
LABEL = 'EXPLORATORY / CONDITIONAL — NOT A COMPLETE PHYSICAL ONE-LOOP NEUTRINO-MASS PREDICTION'
VEV_TO_EV = 0.5 * 246.22**2 * 1.0e9
COMPARISON_COLUMNS = (
    'model', 'direct_status', 'eft_status', 'cross_method_comparison_status',
    'direct_scale_gev', 'eft_mu_f_gev', 'eft_mu_s_gev',
    'direct_m1_ev', 'direct_m2_ev', 'direct_m3_ev',
    'eft_m1_ev_per_abs_N', 'eft_m2_ev_per_abs_N', 'eft_m3_ev_per_abs_N',
    'uv_max_ms2_over_mf2', 'uv_leading_heavy_B0_tail_norm_ev',
    'uv_leading_heavy_B0_tail_rel_to_direct',
    'uv_pole_relative', 'uv_source_conventions_certified',
    'cross_method_relative_matrix_difference', 'blocker'
)


def require(ok: Any, message: str) -> None:
    if not bool(ok):
        raise ValueError(message)


def read_json(path: Path) -> dict:
    value = json.loads(Path(path).read_text(encoding='utf-8-sig'))
    require(isinstance(value,dict),f'expected JSON object: {path}')
    return value


def digest(path: Path) -> str:
    return sha256(Path(path).read_bytes()).hexdigest()


def matrix(obj: Any, name: str) -> np.ndarray:
    require(isinstance(obj,dict) and 'real' in obj and 'imag' in obj,
            f'{name}: missing complex real/imag encoding')
    a = np.asarray(obj['real'],float) + 1j*np.asarray(obj['imag'],float)
    require(a.shape == (3,3) and np.isfinite(a).all(),f'{name}: nonfinite/wrong shape')
    return a


def relative(a: np.ndarray, b: np.ndarray, floor:float=1e-100) -> float:
    return float(np.linalg.norm(a-b)/max(np.linalg.norm(a),np.linalg.norm(b),floor))


def assert_close(a: np.ndarray, b: np.ndarray, name:str,
                 *, rtol:float=5e-10, atol:float=1e-26) -> None:
    scale=max(float(np.linalg.norm(a)),float(np.linalg.norm(b)),1e-30)
    require(float(np.linalg.norm(a-b))<=atol+rtol*scale,
            f'{name}: mismatch (absolute {np.linalg.norm(a-b):.4e})')


def validate_batch2(b: dict) -> dict[str,dict]:
    require(b.get('format')==BATCH2_FORMAT,'wrong Batch 2 format')
    require(b.get('label')==LABEL,'Batch 2 conditional label changed')
    require(b.get('models_checked')==16 and b.get('models_with_conditional_C5')==16
            and b.get('models_blocked')==0,'Batch 2 is not validated 16/16 release')
    require(b.get('physical_C5_certified') is False and
            b.get('absolute_GammaCC_to_Weyl_factor') is None and
            b.get('complete_one_loop_neutrino_mass_prediction') is False,
            'Batch 2 physical-normalization disclaimer missing')
    require(float(b.get('mu_f_gev',0))==1e5 and float(b.get('mu_s_gev',0))==1e3,
            'Batch 2 threshold scales changed')
    rows=b.get('models')
    require(isinstance(rows,list) and len(rows)==16,'need sixteen model entries')
    by={r.get('model'):r for r in rows}
    require(len(by)==16 and set(by)==set(MODELS),'duplicate/missing model in Batch 2')
    for m in MODELS:
        r=by[m]
        require(r.get('status')==SUCCESS and r.get('physical_C5_certified') is False,
                f'{m}: Batch 2 model is not conditional success')
        c=matrix(r.get('conditional_C5_per_N_GeV_inv'),f'{m}.C5')
        mm=matrix(r.get('conditional_mnu_per_N_eV'),f'{m}.mnu')
        assert_close(mm,c*VEV_TO_EV,f'{m}.vev conversion',rtol=3e-10,atol=1e-14)
        assert_close(c,c.T,f'{m}.C5 symmetry')
        sv=np.sort(np.asarray(r['mass_singular_values_eV_per_abs_N'],float))
        require(sv.shape==(3,) and np.isfinite(sv).all() and (sv>=0).all(),
                f'{m}: invalid conditional singular values')
        require(np.allclose(np.sort(np.linalg.svd(mm,compute_uv=False)),sv,
                               rtol=3e-9,atol=1e-12),f'{m}: inconsistent singular values')
        diag=r.get('diagnostics',{})
        require(diag.get('all_entries_finite') is True and
                diag.get('SU2_source_beta_parity') is True and
                float(diag.get('weinberg_projection_residual',1))<=1e-8 and
                float(diag.get('matching_scale_cancellation_residual',1))<=2e-9,
                f'{m}: conditional numerical diagnostics not passed')
    return by


def validate_direct_input(data:dict,model:str,root:Path) -> dict:
    require(data.get('model')==model,f'{model}: direct source model mismatch')
    prov=data.get('provenance')
    require(isinstance(prov,dict),f'{model}: direct input lacks provenance')
    # The current repository's component adapter is explicitly NOT certified.
    require(prov.get('validated_uv_component_conventions') is False,
            f'{model}: unexpected source-convention status; re-audit before comparison')
    require(prov.get('cg_generator_covariance_checked') is True,
            f'{model}: native CG covariance check absent')
    bench=prov.get('benchmark_path')
    require(isinstance(bench,str) and bench and '..' not in Path(bench.replace('\\','/')).parts,
            f'{model}: invalid frozen benchmark path')
    frozen=(root/bench.replace('\\','/')).resolve()
    require(frozen.is_relative_to(root.resolve()),f'{model}: benchmark path escapes root')
    require(frozen.is_file(),f'{model}: frozen source benchmark not present: {frozen}')
    require(digest(frozen)==prov.get('benchmark_sha256'),
            f'{model}: frozen source benchmark hash differs from UV input')
    frozen_state=read_json(frozen).get('base_state',{})
    require(float(frozen_state.get('mu_gev',0))==1e7 and
            float(data['renormalization_scale_gev'])==1e7,
            f'{model}: UV input no longer at frozen 1e7 GeV base_state; cannot use cached provenance')
    require(isinstance(data.get('vertices'),list) and len(data['vertices'])>=2,
            f'{model}: missing source-derived UV vertices')
    return {'benchmark_path':str(frozen.relative_to(root.resolve())),
            'benchmark_sha256':prov['benchmark_sha256'],
            'source_component_mapping_certified':False,
            'unverified_items':prov.get('unverified_items',[])}


def tail_of_uv_B0_expansion(direct:dict) -> dict:
    """Compute the exact vs leading-mS^2/MF^2 B0 tail at the SAME UV point.

    B0_exact=1-ln(M^2/mu^2)+r*ln(r)/(1-r), r=mS^2/M^2.
    B0_leading=1-ln(M^2/mu^2)+r*ln(r).
    Delta=B0_exact-B0_leading = r^2 ln(r)/(1-r).

    Uses the direct engine's individual per-eigenstate terms, NOT EFT data.
    This is only the B0 heavy-mass tail, not a full UV-vs-EFT correction.
    """
    terms=direct.get('individual_terms')
    require(isinstance(terms,list) and terms,'direct engine did not retain individual terms')
    mu=float(direct['renormalization_scale_gev'])
    excess=np.zeros((3,3),complex)
    count=0
    max_r=0.
    for term in terms:
        mf=float(term['fermion_mass_gev'])
        ms=float(term['scalar_mass_gev'])
        require(mf>0 and ms>0 and math.isfinite(mf) and math.isfinite(ms),
                'UV eigenmass not finite/positive')
        r=(ms/mf)**2
        require(0<r<0.1,
                'leading-heavy B0 expansion not controlled for this UV point')
        exact=direct_core.b0_finite(mf**2,ms**2,mu**2)
        require(abs(exact)>1e-7,'B0 zero invalidates term-based tail reconstruction')
        omitted=(r*r*math.log(r)/(1-r))
        excess += matrix(term['matrix_gev'],'individual UV term')*(omitted/exact)
        count+=1
        max_r=max(max_r,r)
    exact_mass=matrix(direct['total']['complex_symmetric_matrix_gev'],'direct UV matrix')
    leading=exact_mass-excess
    return {'expansion_parameter_max_ms2_over_mf2':max_r,
            'B0_terms':count,
            'B0_difference_exact_minus_leading_mass_matrix_gev':direct_core.encoded(excess),
            'direct_uv_leading_heavy_B0_mass_matrix_gev':direct_core.encoded(leading),
            'tail_frobenius_ev':float(np.linalg.norm(excess)*1e9),
            'tail_over_exact_direct_mass_frobenius':float(np.linalg.norm(excess)/max(np.linalg.norm(exact_mass),1e-100)),
            'meaning':'Remainder of expansion of the SAME UV loop function; not missing EFT hard matching or physical UV/EFT discrepancy'}


def evaluate_uv(model:str, input_path:Path, root:Path,
                previous_path:Path|None=None) -> dict:
    raw=read_json(input_path)
    provenance=validate_direct_input(raw,model,root)
    direct=direct_core.evaluate(raw,retain_individual_terms=True)
    require(direct.get('status')=='uv_mass_evaluated_unverified_source_mapping' and
            direct.get('model')==model,f'{model}: unexpected UV engine status')
    d=direct.get('diagnostics',{})
    require(d.get('uv_finite') is True and d.get('neutrino_mass_matrix_symmetric') is True,
            f'{model}: UV pole/symmetry test failed')
    require(float(d.get('relative_uv_pole_coefficient',1))<1e-8,
            f'{model}: nonzero uncancelled UV pole')
    U=matrix(direct['total']['complex_symmetric_matrix_gev'],'direct UV matrix')
    assert_close(U,U.T,f'{model}: UV symmetric matrix',rtol=1e-10)
    masses=np.sort(np.asarray(direct['total']['neutrino_masses_ev'],float))
    require(masses.shape==(3,) and np.all(np.isfinite(masses)) and
            np.all(masses>=0),f'{model}: invalid direct spectrum')
    require(np.allclose(masses,np.sort(np.linalg.svd(U,compute_uv=False))*1e9,
                               rtol=3e-8,atol=1e-8),
            f'{model}: direct mass singular values inconsistent')
    previous=None
    if previous_path is not None and previous_path.is_file():
        old=read_json(previous_path)
        require(old.get('model')==model and
                old.get('source_provenance',{}).get('benchmark_sha256')==provenance['benchmark_sha256'],
                f'{model}: archived direct result does not match frozen input')
        before=matrix(old['total']['complex_symmetric_matrix_gev'],'previous direct matrix')
        assert_close(U,before,f'{model}: direct UV regression',rtol=2e-6,atol=1e-16)
        previous={'archived_direct_result_sha256':digest(previous_path),
                  'recomputed_vs_existing_relative_matrix_residual':relative(U,before)}
    tail=tail_of_uv_B0_expansion(direct)
    return {'model':model,'status':direct['status'],
            'source_provenance':provenance,
            'input_file_sha256':digest(input_path),
            'direct_scale_gev':float(direct['renormalization_scale_gev']),
            'loop_order':1,'loop_momentum':'p=0',
            'scheme':direct['subtraction'],
            'mass_matrix_gev':direct_core.encoded(U),
            'mass_singular_values_ev':masses.tolist(),
            'fermion_eigenmasses_gev':direct['fermion_eigenmasses_gev'],
            'scalar_eigenmasses_gev':direct['scalar_eigenmasses_gev'],
            'max_scalar_squared_over_fermion_squared':direct['mass_hierarchy']['max_scalar_squared_to_fermion_squared_ratio'],
            'relative_uv_pole_coefficient':d['relative_uv_pole_coefficient'],
            'uv_finite_in_declared_input':True,
            'uv_scalar_diagonalization_residual':direct['scalar_mass_diagonalization_residual'],
            'uv_fermion_takagi_residual':direct['fermion_takagi_residual'],
            'existing_uv_result_crosscheck':previous,
            'leading_heavy_mass_B0_tail':tail,
            'physical_source_component_mapping_certified':False}


def compare(model:str, eft:dict, uv:dict) ->dict:
    require(model==eft['model']==uv['model'],'mismatched comparison labels')
    # Distinct reference datasets: source UV frozen fitted_reference base_state
    # at 1e7 GeV, conditional EFT from illustrative running @1e3 GeV.
    e=matrix(eft['conditional_mnu_per_N_eV'],'EFT mnu per N')
    u=matrix(uv['mass_matrix_gev'],'direct UV mnu')
    masses_e=np.sort(np.linalg.svd(e,compute_uv=False))
    masses_u=np.sort(np.linalg.svd(u,compute_uv=False)*1e9)
    require(np.allclose(masses_e,np.sort(eft['mass_singular_values_eV_per_abs_N']),
                               rtol=3e-9,atol=1e-12),'EFT singular value regression')
    require(np.allclose(masses_u,uv['mass_singular_values_ev'],
                               rtol=3e-8,atol=1e-8),'UV singular value regression')
    blockers=[
        'Different renormalized parameter points: frozen UV base_state vs illustrative EFT1 @muS',
        'UV direct input scale is not the EFT scalar matching scale',
        'Direct UV two-component vertex phases/normalization not source-certified',
        'EFT absolute signed N_GammaCC-to-Weyl normalization unknown',
        'Matchete-to-RGBeta quartic aliases and special scalar mass basis not independently certified',
        'EFT lacks hard heavy-fermion one-loop matching, additional EFT1 mixing and finite counterterms',
        'Scalar inputs differ (direct EWSB mixed masses vs EFT unbroken illustrative endpoint)',
        'Common charged-lepton flavour basis and renormalization convention not independently certified'
    ]
    return {'model':model,
            'status':'two_method_records_available_but_not_same_parameter_point',
            'cross_method_comparison_status':'not_comparable_unmatched_parameter_point_and_conventions',
            'direct_uv':uv,
            'eft':{
                'source':'Batch 2 validated conditional all16 release',
                'scalar_case':eft['scalar_case'],
                'mu_f_gev':eft['mu_f_gev'],'mu_s_gev':eft['mu_s_gev'],
                'mass_matrix_eV_per_N':eft['conditional_mnu_per_N_eV'],
                'C5_matrix_GeV_inverse_per_N':eft['conditional_C5_per_N_GeV_inv'],
                'mass_singular_values_eV_per_abs_N_ascending':masses_e.tolist(),
                'scalar_input_masses_squared_gev2':eft['scalar_mass_inputs_gev2'],
                'scalar_eigenstate_masses_gev':eft['scalar_eigenstate_masses_gev'],
                'physical_C5_certified':False,
                'absolute_GammaCC_normalization':None},
            'direct_uv_mass_singular_values_eV_ascending':masses_u.tolist(),
            'eFT_conditional_mass_singular_values_eV_per_abs_N_ascending':masses_e.tolist(),
            'matrix_difference_gev':None,
            'relative_uv_vs_eft_matrix_difference':None,
            'ratio_uv_vs_eft_mass_singular_values':None,
            'physical_equivalence_claim':False,
            'finite_heavy_mass_effects':uv['leading_heavy_mass_B0_tail'],
            'blockers_to_direct_physical_comparison':blockers,
            'terms_in_direct_uv':['one-loop p=0 all supplied Weyl/scalar eigenstates',
                                  'exact finite B0 function and UV mass mixing at frozen point'],
            'terms_in_eft':['source-resolved LLSS provisional boundary',
                            'fixed-order LLSS-induced EFT1 quartic leading log',
                            'simultaneous finite two-multiplet scalar bubble at muS'],
            'omitted_from_eft':['hard one-loop fermion threshold',
                                'full EFT1 Wilson operator mixing and finite subtractions',
                                'independently certified source normalization/physical couplings'],
            'interpretation':LABEL}


def row_for(model:str,record:dict) ->dict:
    out=dict.fromkeys(COMPARISON_COLUMNS,'')
    out['model']=model
    out['eft_status']='conditional_C5_per_unknown_NGamma_exploratory'
    out['cross_method_comparison_status']=record['cross_method_comparison_status']
    if record.get('direct_uv') is None:
        out['direct_status']='ERROR_undelivered_UV'
        out['blocker']=record.get('blocker','UV input blocked')
        return out
    direct=record['direct_uv'];eft=record['eft'];tail=direct['leading_heavy_mass_B0_tail']
    out.update({'direct_status':direct['status'],
                'direct_scale_gev':direct['direct_scale_gev'],
                'eft_mu_f_gev':eft['mu_f_gev'],'eft_mu_s_gev':eft['mu_s_gev'],
                'uv_max_ms2_over_mf2':direct['max_scalar_squared_over_fermion_squared'],
                'uv_leading_heavy_B0_tail_norm_ev':tail['tail_frobenius_ev'],
                'uv_leading_heavy_B0_tail_rel_to_direct':tail['tail_over_exact_direct_mass_frobenius'],
                'uv_pole_relative':direct['relative_uv_pole_coefficient'],
                'uv_source_conventions_certified':'False',
                'cross_method_relative_matrix_difference':'',
                'blocker':' | '.join(record['blockers_to_direct_physical_comparison'])})
    for prefix,vals in [('direct',record['direct_uv_mass_singular_values_eV_ascending']),
                        ('eft',record['eFT_conditional_mass_singular_values_eV_per_abs_N_ascending'])]:
        for i,val in enumerate(vals,1):
            out[f'{prefix}_m{i}_ev'+('_per_abs_N' if prefix=='eft' else '')]=float(val)
    return out


def build_release(batch2_dir:Path,uv_input_dir:Path,previous_dir:Path|None,
                  output_dir:Path,project_root:Path=ROOT) -> dict:
    batch2_path=batch2_dir/'summary.json'
    batch2=read_json(batch2_path)
    by=validate_batch2(batch2)
    comparisons=[]
    for model in TARGETS:
        try:
            previous=(previous_dir/(model+'_direct_one_loop.json')) if previous_dir else None
            uv=evaluate_uv(model,uv_input_dir/(model+'.json'),project_root,previous)
            record=compare(model,by[model],uv)
        except (OSError,KeyError,ValueError,TypeError,RuntimeError,IndexError) as exc:
            record={'model':model,'status':'ERROR_direct_uv_input_or_evaluation',
                    'cross_method_comparison_status':'not_compared_UV_blocked',
                    'blocker':f'{type(exc).__name__}: {exc}',
                    'direct_uv':None,'physical_equivalence_claim':False,
                    'relative_uv_vs_eft_matrix_difference':None,
                    'interpretation':LABEL}
        comparisons.append(record)
    n=sum(r.get('direct_uv') is not None for r in comparisons)
    return {
        'format':'t3-final-conditional-all16-plus-two-direct-uv-batch3-v1',
        'label':LABEL,
        'status':'final_exploratory_numerical_release' if n==2 else 'partial_with_explicit_direct_uv_blockers',
        'all16_batch2_input_summary_sha256':digest(batch2_path),
        'all16_batch2_input_comparison_csv_sha256':digest(batch2_dir/'all16_conditional_mass_summary.csv'),
        'all16_models_checked':16,'all16_conditional_models':16,
        'all16_blocked':0,'direct_uv_models_requested':list(TARGETS),
        'direct_uv_evaluated':n,'direct_uv_blocked':len(TARGETS)-n,
        'physically_same_point_UV_EFT_comparisons':0,
        'cross_method_relative_differences_reported':0,
        'complete_one_loop_physical_neutrino_mass_predictions':0,
        'matching_sequence':'UV -> ALL BSM fermions at 100000 GeV -> EFT1(S1,S2) -> BOTH BSM scalar multiplets at 1000 GeV -> SMEFT',
        'all16_conditional_results':[by[k] for k in MODELS],
        'direct_vs_eft_comparisons':comparisons,
        'validation':{
            'batch2_all16_matrix_symmetry_vev_singular_values_and_ward_checked':True,
            'direct_core_is_existing_repository_engine':True,
            'direct_finite_B0_exact_and_large_mass_tail_evaluated':n,
            'frozen_UV_benchmark_sha256_validated_for_successful_direct_models':n,
            'direct_archived_regression_checked_when_available':True,
            'direct_UV_component_source_conventions_certified':False,
            'Matchete_to_RGBeta_quartic_map_certified':False,
            'signed_GammaCC_normalization_known':False
        },
        'next_step':'STOP: three-batch exploratory numerical milestone reached; deferred certification is not part of this release'
    }


def main(argv:list[str]|None=None)->int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--batch2-dir',type=Path,default=INPUT2)
    p.add_argument('--uv-input-dir',type=Path,default=UV_INPUT)
    p.add_argument('--previous-direct-dir',type=Path,default=UV_PREVIOUS)
    p.add_argument('--output-dir',type=Path,default=OUTPUT)
    args=p.parse_args(argv)
    output=args.output_dir.resolve();safe=SAFE.resolve()
    if output==safe or safe not in output.parents:
        p.error(f'new output must be under {safe}')
    if output.exists() and any(output.iterdir()):
        p.error('destination is nonempty; use a fresh release directory')
    if not (args.batch2_dir/'summary.json').is_file():
        p.error('validated Batch 2 summary missing')
    if not (args.batch2_dir/'all16_conditional_mass_summary.csv').is_file():
        p.error('Batch 2 consolidated CSV missing')
    final=build_release(args.batch2_dir,args.uv_input_dir,args.previous_direct_dir,
                        output,ROOT)
    output.mkdir(parents=True,exist_ok=True)
    # Retain the final all-16 table verbatim. Cross-UV comparisons never modify it.
    shutil.copyfile(args.batch2_dir/'all16_conditional_mass_summary.csv',
                    output/'all16_conditional_mass_summary.csv')
    for record in final['direct_vs_eft_comparisons']:
        name=record['model'].replace('-','m').replace('+','p')
        (output/f'direct_uv_vs_eft_{name}.json').write_text(
            json.dumps(record,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    with (output/'direct_uv_comparison.csv').open('w',newline='',encoding='utf-8') as fh:
        writer=csv.DictWriter(fh,fieldnames=COMPARISON_COLUMNS)
        writer.writeheader()
        writer.writerows(row_for(record['model'],record) for record in final['direct_vs_eft_comparisons'])
    (output/'summary.json').write_text(json.dumps(final,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    validation=dict(final['validation'])
    validation.update({'label':LABEL,
        'batch1_tests_reported_by_user':'58/58 passed',
        'batch2_tests_reported_by_user':'4/4 passed',
        'batch2_models_reported_by_user':'16/16 conditional C5, 0 blocked',
        'batch3_unit_tests':'run separately using supplied PowerShell command',
        'direct_uv_models':[{
            'model':x['model'],'status':x['status'],
            'input_sha256':x['direct_uv']['input_file_sha256'] if x.get('direct_uv') else None,
            'uv_finite':x['direct_uv']['uv_finite_in_declared_input'] if x.get('direct_uv') else False,
            'max_ms2_over_mf2':x['direct_uv']['max_scalar_squared_over_fermion_squared'] if x.get('direct_uv') else None,
            'source_mapping_certified':False
        } for x in final['direct_vs_eft_comparisons']],
        'scientifically_comparable_UV_EFT_points':0})
    (output/'validation.json').write_text(json.dumps(validation,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    for record in final['direct_vs_eft_comparisons']:
        print(record['model']+': '+record['status'])
        if 'blocker' in record: print('  '+record['blocker'])
    print(LABEL)
    print(f"Batch 3: all16=16/16; provisional direct UV={final['direct_uv_evaluated']}/2; same-point physical comparisons=0")
    print('JSON:',output/'summary.json')
    print('CSV:',output/'all16_conditional_mass_summary.csv')
    print('UV comparisons:',output/'direct_uv_comparison.csv')
    return 0 if final['direct_uv_evaluated']==2 else 2

if __name__=='__main__':
    raise SystemExit(main())
