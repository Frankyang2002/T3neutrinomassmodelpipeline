"""All-16 conditional native Matchete finite LLHH branch, strict EFT ordering.

Uses saved native PL constant+pole+log groups at MF, the all-16 running state,
V3 conditional quartic coordinate transformation, and existing Batch-2 strict
Wilson/Workstream-C solver. The already-included analytic fermion-hard triangle
is REPLACED, never added a second time. This native F-only MatchedEFT pole-drop
is a reproducible comparison SCENARIO, not certified MSbar hard subtraction.

UV -> ALL F integrated out at common MF -> EFT1 with BOTH scalars ->
BOTH scalar multiplets integrated out at common MS -> SMEFT -> masses.
"""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
import shutil
import sys

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
import numpy as np

from Numerical.parallel_t3.wilson.basis import ALL_MODELS,build_basis
from Numerical.parallel_t3.background.contracts import tag
from Numerical.parallel_t3.wilson.native_finite_batch4 import (
    native_all16_boundary,relative,digest,NATIVE_SCHEMA)
from Numerical.parallel_t3.wilson.coordinates import embed_c5,project_c5,symmetry_residual
from Numerical.parallel_t3.threshold.preflight import preflight_pair
from Numerical.parallel_t3.threshold.pipeline import run_case
from Numerical.parallel_t3.threshold.rough_numerics import neutrino_observables

SCHEMA='t3-native-matchete-pole-drop-conditional-all16-batch4-v1'
STATUS='NUMERICAL_NATIVE_MATCHeTE_FINITE_SCENARIO_NOT_PHYSICAL_MSbar'
BASE=Path('output/parallel_t3/strict_one_loop_all16_batch2_v1')
OUTPUT=Path('output/parallel_t3/native_finite_scenario_all16_batch4_v1')


def read_json(path:Path)->dict:
    d=json.loads(path.read_text(encoding='utf-8-sig'))
    if not isinstance(d,dict):raise ValueError('Expected JSON object '+str(path))
    return d


def read_npz(path:Path)->dict:
    with np.load(path,allow_pickle=False) as f:
        return {k:np.asarray(f[k]).copy() for k in f.files}


def save_npz(path:Path,d:dict)->None:
    path.parent.mkdir(parents=True,exist_ok=True)
    p=path.with_suffix('.partial.npz')
    np.savez_compressed(p,**d)
    p.replace(path)


def save_json(path:Path,d:dict)->None:
    path.parent.mkdir(parents=True,exist_ok=True)
    p=path.with_suffix('.partial.json')
    p.write_text(json.dumps(d,sort_keys=True,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    p.replace(path)


def jmat(x:np.ndarray)->dict:
    a=np.asarray(x,complex)
    return {'real':a.real.tolist(),'imag':a.imag.tolist(),
            'norm_gev_inverse':float(np.linalg.norm(a))}


def replace_strict_hard(strict_wilson:np.ndarray,old_hard:np.ndarray,
                        native_hard:np.ndarray,basis, *,running_delta:np.ndarray,
                        zero_order_full:np.ndarray) -> tuple[np.ndarray,np.ndarray,dict]:
    """Replace, never sum, the alternative MF finite one-loop C5.

    This operation changes the complete full tensor's Weinberg HH block but
    preserves all order-zero LLSS, other sectors, and the precomputed O(hbar)
    beta on C0. In strict counting, the new hard boundary is not re-evolved.
    """
    old=np.asarray(old_hard,complex)
    new=np.asarray(native_hard,complex)
    if old.shape!=(3,3) or new.shape!=(3,3) or not np.isfinite(new).all():
        raise ValueError('Invalid finite hard C5 input')
    if not np.allclose(new,new.T,rtol=1e-10,atol=1e-25):
        raise ValueError('Native C5 has wrong flavour symmetry')
    s=np.asarray(strict_wilson,complex)
    c0=np.asarray(zero_order_full,complex)
    if s.shape!=c0.shape or c0.shape!=(basis.nf,basis.nf,basis.ns,basis.ns):
        raise ValueError('Wrong full dimension-five Wilson tensor')
    original,_=project_c5(s,basis)
    initial,_=project_c5(c0,basis)
    if np.linalg.norm(initial)>1e-21:
        raise ValueError('Reference zero-order EFT1 has a separate C5 seed')
    d=embed_c5(new-old,basis)
    new_full=s+d
    c5,_=project_c5(new_full,basis)
    expected=np.asarray(running_delta,complex)+new
    err=relative(c5,expected)
    if err>1e-9 or symmetry_residual(new_full)>max(1e-24,1e-9*np.max(abs(new_full))):
        raise ValueError('Alternative hard C5 plus EFT1 one-loop mixing closure failed')
    # Entire non-HH Wilson sector must remain exactly as before.
    outside=np.ones((basis.ns,basis.ns),bool);outside[:4,:4]=False
    if np.linalg.norm((new_full-s)[:,:,outside])>1e-10*max(np.linalg.norm(d),1e-100):
        raise ValueError('Alternative native hard matching altered non-Weinberg operators')
    return new_full,c5,{
        'old_analytic_triangle_MF_C5':jmat(old),
        'replacement_native_pole_drop_MF_C5':jmat(new),
        'relative_replacement_shift_to_old':relative(new,old),
        'strict_total_C5_at_MS':jmat(c5),
        'strict_C5_running_delta_unchanged':jmat(running_delta),
        'strict_C5_closure_residual':err,
        'other_EFT1_operators_unchanged':True,
        'hard_triangles_added_to_native_source':False}


def run_model(root:Path,base:Path,out:Path,model:str,*,
              native_spinor_factor:float,mu_low:float,vev:float,
              audit_only:bool)->dict:
    slug=tag(model)
    row={'model':model,'status':'BLOCKED','stage':'input'}
    try:
        prefix=base/'inputs'/slug
        bw=read_json(prefix/'wilson.json')
        bb=read_json(prefix/'background.json')
        wz=read_npz(prefix/'wilson_at_MS.npz')
        sb=read_npz(prefix/'strict_wilson_channels.npz')
        if (bw.get('model')!=model or bb.get('model')!=model or
                bw.get('strict_1loop_endpoint_at_ms') is not True or
                bw.get('tensor_npz_sha256')!=digest(prefix/'wilson_at_MS.npz') or
                bb.get('tensor_npz_sha256')!=digest(prefix/'background.npz')):
            raise ValueError('Batch-2 strict Wilson/scalar input invalid')
        oldhard=sb['C5_MF_hard_fermion_triangle']
        if oldhard.shape!=(3,3):raise ValueError('Missing Batch-2 triangle')
        if relative(sb['C_psi2phi2_ms_strict_1loop'],wz['C_psi2phi2_ms_strict_1loop'])>1e-11:
            raise ValueError('Batch-2 full strict tensor does not match its saved channel ledger')
        row['stage']='source_native_matchete_one_loop'
        native,per,source=native_all16_boundary(root,model,native_spinor_factor=native_spinor_factor)
        if source.get('schema')!=NATIVE_SCHEMA or source['native_PL_kernel_count']<=0:
            raise ValueError('Native source did not produce a complete candidate')
        row['stage']='strict_one_loop_hard_replacement'
        basis=build_basis(model)
        full,c5,comparison=replace_strict_hard(
            sb['C_psi2phi2_ms_strict_1loop'],oldhard,native,basis,
            running_delta=sb['C5_EFT1_strict_running_delta'],
            zero_order_full=sb['C_psi2phi2_tree0_mf'])
        folder=out/'inputs'/slug
        folder.mkdir(parents=True,exist_ok=True)
        save_npz(folder/'native_MF_LLHH_finite_source_candidate.npz',{
            'C5_native_finite_MF':native,
            'C5_native_by_kernel_and_heavy':per,
            'C5_analytic_triangle_replaced_MF':oldhard,
            'C5_native_minus_triangle_MF':native-oldhard,
            'C5_native_strict_EFT1_MS':c5,
            'C5_EFT1_running_delta_unchanged':sb['C5_EFT1_strict_running_delta']})
        save_json(folder/'native_MF_LLHH_finite_source_candidate.json',source)
        save_json(folder/'strict_hard_replacement_ledger.json',comparison)
        z=dict(wz)
        z['C_psi2phi2_ms_strict_1loop']=full
        z['C5_EFT1_MS_strict_1loop']=c5
        wp=folder/'wilson_at_MS.npz'
        save_npz(wp,z)
        wm=dict(bw)
        wm.update({'tensor_npz_sha256':digest(wp),
           'strict_1loop_origin':('tree0_MF + native_Matchete_FORMAL_pole_drop_C5_MF '
                  '+ ln(MS/MF)/(16pi2)*candidate_full_EFT1_beta(tree0_MF)'),
           'strict_native_F_finite_candidate_replaces_F_triangle':True,
           'original_Batch2_resummed_C5_is_DIAGNOSTIC_NOT_this_strict_branch':True,
           'native_source_MSbar_and_soft_subtraction_certified':False,
           'native_to_master_signed_spinor_factor':native_spinor_factor,
           'strict_one_loop_with_included_topologies_only':True,
           'full_one_loop_physical_matching_certified':False})
        save_json(folder/'wilson.json',wm)
        for name in ('background.npz','background.json'):
            shutil.copy2(prefix/name,folder/name)
        row.update({'native_F_PL_kernels':source['native_PL_kernel_count'],
             'native_MF_C5':jmat(native),'strict_EFT1_MS_C5':jmat(c5),
             'hard_replacement':comparison,
             'original_batch2_strict_summary_sha256':digest(base/'summary.json'),
             'source_shas':source['source_sha256'],
             'native_F_one_loop_is_independently_renormalized':False,
             'strict_EFT1_full_operator_running_retained':True,
             'both_scalar_multiplets_active_between_MF_MS':True})
        row['stage']='C_preflight'
        pf=preflight_pair(folder/'wilson.json',folder/'background.json',
            mode='strict_one_loop',mu_low=mu_low,expected_model=model)
        row['C_preflight_status']=pf.get('status')
        if pf.get('status')!='ready_conditional':
            raise ValueError('Workstream-C strict preflight failed: '+str(pf.get('blockers')))
        if audit_only:
            row['status']='PREPARED_NATIVE_FINITE_SCENARIO'
            row['stage']='prepared'
            return row
        row['stage']='S_joint_threshold_and_SM_Weinberg_running'
        results=out/'results'/slug
        payload=run_case(wilson_manifest=folder/'wilson.json',
              background_manifest=folder/'background.json',
              output_dir=results,mode='strict_one_loop',
              mu_low=mu_low,vev=vev,
              allow_projection_residual=False,match_only=False,overwrite=False)
        if payload.get('SMEFT_running',{}).get('success') is not True:
            raise ValueError('SM+Weinberg RGE did not report numerical success')
        at_ms=read_npz(results/'C5_smeft_at_MS.npz')
        at_low=read_npz(results/'C5_smeft_low.npz')
        if relative(at_ms['C5_SMEFT_MS'],c5+at_ms['C5_scalar_bubble_MS']+
                 at_ms['C5_independent_finite_MS'])>1e-9:
            raise ValueError('Native + strict Weinberg mixing + scalar bubble closure failed')
        if np.linalg.norm(at_ms['C5_independent_finite_MS'])>1e-24:
            raise ValueError('Extra independently matched scalar C5 term unexpectedly included')
        b=read_npz(prefix/'background.npz')
        obs,matrices=neutrino_observables(at_low['m_nu_GeV'],b['ye'])
        if not np.allclose(obs['masses_eV'],at_low['mass_singular_values_eV'],rtol=1e-8,atol=1e-13):
            raise ValueError('PMNS/Takagi masses mismatch with existing SMEFT solver')
        save_npz(results/'PMNS_matrices.npz',matrices)
        obs.update({'status':STATUS,'conditional_not_physical':True,
               'source_matching':'native_finite_formal_pole_drop_NOT_MSbar_hard',
               'native_F_finite_replaces_analytic_triangle':True,
               'one_common_MF_all_fermions':True,'one_common_MS_all_scalars':True})
        save_json(results/'observables.json',obs)
        row.update({'status':STATUS,'stage':'complete',
               'conditional_masses_eV':obs['masses_eV'],
               'conditional_dm21_sq_eV2':obs['dm21_sq_eV2'],
               'conditional_dm31_sq_eV2':obs['dm31_sq_eV2'],
               'C5_smeft_at_MS':jmat(at_ms['C5_SMEFT_MS']),
               'C5_scalar_bubble_at_MS':jmat(at_ms['C5_scalar_bubble_MS']),
               'finite_new_F_matching_MSbar_certified':False,
               'full_one_loop_matching_completed':False})
    except Exception as exc:
        row['status']='BLOCKED'
        row['blocker']=type(exc).__name__+': '+str(exc)
        print(model,'BLOCKED',row['stage'],row['blocker'],file=sys.stderr)
    return row


def comparisons(out:Path,base:Path,report:dict)->None:
    src=read_json(base/'summary.json')
    rows=[]
    for model in ALL_MODELS:
        x=report['results'].get(model,{})
        if x.get('status')!=STATUS:continue
        old=src.get('results',{}).get(model,{})
        a=x['conditional_masses_eV']
        b=old.get('masses_eV_CONDITIONAL')
        row={'model':model,'native_PL_kernels':x['native_F_PL_kernels'],
             'native_C5_MF_norm':x['native_MF_C5']['norm_gev_inverse'],
             'Batch2_triangle_C5_MF_norm':x['hard_replacement']['old_analytic_triangle_MF_C5']['norm_gev_inverse']}
        for i in range(3):
            row[f'm{i+1}_native_formal_eV']=a[i]
            row[f'm{i+1}_batch2_strict_eV']=b[i] if b else ''
            row[f'm{i+1}_delta_eV']=a[i]-b[i] if b else ''
        rows.append(row)
    if rows:
        with (out/'batch4_native_vs_batch2_strict_all16.csv').open('w',newline='',encoding='utf-8') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def run_batch(root:Path,base:Path,output:Path,models:list[str],*,
              native_spinor_factor:float=1.,mu_low:float=173.,vev:float=246.22,
              overwrite:bool=False,resume:bool=False,audit_only:bool=False)->dict:
    root=Path(root).resolve()
    base=base if base.is_absolute() else root/base
    out=output if output.is_absolute() else root/output
    if out.resolve() in (root,base.resolve()):raise ValueError('Output must not replace baseline')
    if overwrite and resume:raise ValueError('Use one of --overwrite or --resume')
    baseline=read_json(base/'summary.json')
    if baseline.get('format')!='t3-strict-one-loop-all16-hard-triangle-and-scalar-bubble-v1':
        raise ValueError('Need completed all-16 strict Batch 2 source')
    for model in models:
        if (baseline.get('results',{}).get(model,{}).get('status') !=
                'NUMERICAL_STRICT_ONE_LOOP_INCLUDED_TOPOLOGIES_CANDIDATE_NOT_CERTIFIED'):
            raise ValueError('Model not complete in saved Batch 2: '+model)
    if out.exists() and any(out.iterdir()):
        if not resume and not overwrite:raise FileExistsError('Output exists: --resume or --overwrite')
        old=read_json(out/'summary.json')
        if old.get('format')!=SCHEMA:raise ValueError('Unrecognized output directory; will not overwrite')
        if overwrite:
            shutil.rmtree(out);summary=None
        else:summary=old
    else:summary=None
    out.mkdir(parents=True,exist_ok=True)
    if summary is None:
        summary={'format':SCHEMA,'results':{},'baseline_batch2_summary_sha256':digest(base/'summary.json'),
            'candidate_native_spinor_factor':native_spinor_factor,
            'mu_low_gev':mu_low,'vev_gev':vev,
            'muF_gev':100_000.,'muS_gev':1_000.,
            'EFT_order':'UV -> all F out -> EFT1 with S1,S2 -> all S out -> SMEFT',
            'native_matchete_F_one_loop_formal_pole_drop':True,
            'old_Batch1_triangle_replaced_not_added':True,
            'native_F_EFT1_soft_and_counterterms_subtracted':False,
            'GammaCC_Weyl_normalization_physical':False,
            'strict_EFT1_one_loop_Wilson_mixing':True,
            'scalar_pair_bubble_one_loop':True,
            'complete_one_loop_matching':False,
            'physical_neutrino_masses':False}
    else:
        if summary['baseline_batch2_summary_sha256']!=digest(base/'summary.json') or \
          summary['candidate_native_spinor_factor']!=native_spinor_factor or \
          summary['mu_low_gev']!=mu_low or summary['vev_gev']!=vev:
            raise ValueError('Cannot resume with changed source or benchmark assumptions')
    done={STATUS,'PREPARED_NATIVE_FINITE_SCENARIO'} if audit_only else {STATUS}
    selected=[m for m in models if not resume or summary['results'].get(m,{}).get('status') not in done]
    if resume:print('Resume selected:',', '.join(selected) if selected else '(none)')
    save_json(out/'summary.json',summary)
    for model in selected:
        for dirname in ('inputs','results'):
            p=out/dirname/tag(model)
            if p.is_dir():shutil.rmtree(p)
        result=run_model(root,base,out,model,native_spinor_factor=native_spinor_factor,
                         mu_low=mu_low,vev=vev,audit_only=audit_only)
        summary['results'][model]=result
        save_json(out/'summary.json',summary)
        if result['status']==STATUS:
            print(model,result['status'],'masses_eV:',result['conditional_masses_eV'])
        elif result['status']=='PREPARED_NATIVE_FINITE_SCENARIO':
            print(model,'NATIVE ONE-LOOP SOURCE EVALUATED; STRICT C PREFLIGHT READY')
    summary['complete_count']=sum(x.get('status')==STATUS for x in summary['results'].values())
    summary['prepared_count']=sum(x.get('status')=='PREPARED_NATIVE_FINITE_SCENARIO' for x in summary['results'].values())
    summary['blocked_count']=sum(x.get('status')=='BLOCKED' for x in summary['results'].values())
    summary['executed_this_call']=selected
    if not audit_only:comparisons(out,base,summary)
    save_json(out/'summary.json',summary)
    return summary


def show(output:Path)->None:
    s=read_json(output/'summary.json')
    if s.get('format')!=SCHEMA:raise ValueError('Not a Batch 4 native scenario')
    print('Model    Native PL   Status                           Conditional masses (eV)')
    print('-'*114)
    for model in ALL_MODELS:
        row=s.get('results',{}).get(model,{})
        ms=row.get('conditional_masses_eV')
        value=', '.join(f'{a:.8g}' for a in ms) if ms else row.get('blocker','not run')
        print(f"{model:<8} {row.get('native_F_PL_kernels','?')!s:<11} {row.get('status','NOT_RUN'):<33} {value}")
    print('Native pole-drop FORMAL scenario: NOT UV-EFT1-soft-subtracted MSbar matching.')


def main(argv=None)->int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root',type=Path,default=ROOT)
    p.add_argument('--batch2',type=Path,default=BASE)
    p.add_argument('--output',type=Path,default=OUTPUT)
    p.add_argument('--model',choices=ALL_MODELS)
    p.add_argument('--native-spinor-factor',type=float,default=1.)
    p.add_argument('--mu-low',type=float,default=173.)
    p.add_argument('--vev',type=float,default=246.22)
    p.add_argument('--overwrite',action='store_true')
    p.add_argument('--resume',action='store_true')
    p.add_argument('--audit-only',action='store_true')
    p.add_argument('--show',action='store_true')
    a=p.parse_args(argv)
    try:
        out=a.output if a.output.is_absolute() else a.project_root/a.output
        if a.show:show(out);return 0
        s=run_batch(a.project_root,a.batch2,a.output,[a.model] if a.model else list(ALL_MODELS),
            native_spinor_factor=a.native_spinor_factor,mu_low=a.mu_low,vev=a.vev,
            overwrite=a.overwrite,resume=a.resume,audit_only=a.audit_only)
        print('Batch summary:',out/'summary.json','complete=',s['complete_count'],
              'blocked=',s['blocked_count'])
        return 0 if s['blocked_count']==0 else 2
    except Exception as e:
        print('BATCH4 ERROR:',type(e).__name__,str(e),file=sys.stderr)
        return 1


if __name__=='__main__':raise SystemExit(main())
