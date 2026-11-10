"""T3 Batch 2: all-16 strict O(hbar) Wilson + scalar-matching execution.

Loads immutable Batch-1 hard-triangle outputs and original saved Workstream B
UV/EFT1 trajectories. Does NOT rerun or change any previous result.

FULL first-order endpoint for the INCLUDED topologies:
 C(MS) = C_tree(MF) + C_F-hard_triangle(MF)
         + ln(MS/MF)/(16pi²) beta_EFT1[C_tree(MF); g(MF)].
The existing Workstream-C strict_one_loop mode uses C_tree(MF), not C(MS),
for the scalar-bubble term. It then uses the existing SM+Weinberg ODE.

Both MF hard triangles and MS scalar bubbles are decomposed by the same
conditional source-named RGBeta quartic channels. This does not compute the
as-yet-unknown additional fermion-hard diagrams or scalar finite diagrams.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import shutil
import sys
import traceback

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0,str(ROOT))

import numpy as np

from Numerical.parallel_t3.background.contracts import tag
from Numerical.parallel_t3.wilson.basis import ALL_MODELS
from Numerical.parallel_t3.wilson.source_audited_all16 import hash_file
from Numerical.parallel_t3.wilson.rough_stream2 import load_canonical
from Numerical.parallel_t3.wilson.strict_one_loop_batch2 import (
    construct_strict_eft1,SCOPE)
from Numerical.parallel_t3.wilson.quartic_channels_batch2 import (
    named_quartics,split_existing_f_hard)
from Numerical.parallel_t3.threshold.preflight import preflight_pair
from Numerical.parallel_t3.threshold.pipeline import run_case
from Numerical.parallel_t3.threshold.rough_numerics import neutrino_observables

SCHEMA='t3-strict-one-loop-all16-hard-triangle-and-scalar-bubble-v1'
STATUS='NUMERICAL_STRICT_ONE_LOOP_INCLUDED_TOPOLOGIES_CANDIDATE_NOT_CERTIFIED'
BASE=Path('output/parallel_t3/fhard_triangle_all16_batch1_v1')
BASE_INPUT=Path('output/parallel_t3/rough_inputs/fhard_triangle_all16_batch1_v1/a')
OUTPUT=Path('output/parallel_t3/strict_one_loop_all16_batch2_v1')


def read_json(path:Path)->dict:
    value=json.loads(Path(path).read_text(encoding='utf-8-sig'))
    if not isinstance(value,dict):
        raise ValueError('Expected JSON object: '+str(path))
    return value


def put_json(path:Path,value:dict)->None:
    path.parent.mkdir(parents=True,exist_ok=True)
    p=path.with_suffix(path.suffix+'.partial')
    p.write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n',encoding='utf-8')
    p.replace(path)


def read_npz(path:Path)->dict:
    with np.load(path,allow_pickle=False) as f:
        return {name:np.asarray(f[name]).copy() for name in f.files}


def put_npz(path:Path,arrays:dict)->None:
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix('.partial.npz')
    np.savez_compressed(tmp,**arrays)
    tmp.replace(path)


def cmat(z:np.ndarray)->dict:
    z=np.asarray(z,complex)
    return {'real':z.real.tolist(),'imag':z.imag.tolist(),
            'frobenius_gev_inverse':float(np.linalg.norm(z))}


def _max_rel(a:np.ndarray,b:np.ndarray)->float:
    return float(np.linalg.norm(a-b)/max(np.linalg.norm(b),1e-100))


def _verify_first_batch(b:Path,models)->dict:
    s=read_json(b/'summary.json')
    if s.get('format')!='t3-fhard-triangle-conditional-all16-batch1-v1':
        raise ValueError('Input is not the original completed hard-triangle Batch 1')
    records=s.get('results')
    if not isinstance(records,dict):
        raise ValueError('No per-model Batch-1 records')
    for model in models:
        record=records.get(model,{})
        if record.get('status')!='NUMERICAL_FHARD_TRIANGLE_CANDIDATE_NOT_CERTIFIED':
            raise ValueError(model+': Batch-1 calculation not completed')
    return s


def run_model(root:Path,batch1:Path,out:Path,model:str,*,
              mu_low:float,vev:float,preflight_only:bool=False)->dict:
    """Source-backed strict Wilson endpoint and actual C strict matching/SMEFT."""
    record={'model':model,'slug':tag(model),'status':'BLOCKED','stage':'loading'}
    try:
        slug=tag(model)
        prefix=batch1/'adapted'/slug
        wpath=prefix/'wilson.json'
        bpath=prefix/'background.json'
        wm=read_json(wpath)
        bm=read_json(bpath)
        wz=read_npz(prefix/'wilson_at_MS.npz')
        bz=read_npz(prefix/'background.npz')
        if (wm.get('model')!=model or bm.get('model')!=model or
            wm.get('tensor_npz_sha256')!=hash_file(prefix/'wilson_at_MS.npz') or
            bm.get('tensor_npz_sha256')!=hash_file(prefix/'background.npz')):
            raise ValueError('Batch-1 adapted pair source or SHA-256 mismatch')
        if not (wm.get('includes_eft1_weinberg') is True and
                wm.get('adm_variant')=='fonseca_eq22_six_term'):
            raise ValueError('Batch-1 does not have the required full Weinberg-mixed candidate RGE')
        case=load_canonical(root,root/BASE_INPUT/slug,model,c5_boundary='require')
        hardpath=batch1/'hard_fermion'/slug/'fermion_hard_triangle.npz'
        hardmeta=read_json(hardpath.with_name('fermion_hard_triangle.json'))
        if hardmeta.get('model')!=model or hardmeta.get('hard_npz_sha256')!=hash_file(hardpath):
            raise ValueError('Batch-1 hard matching manifest/NPZ SHA-256 mismatch')
        hp=read_npz(hardpath)
        hardc5=np.asarray(hp['C5_hard_F_triangle_MF'],complex)
        if _max_rel(wz['C_psi2phi2_tree_mf'],case.initial)>1e-10:
            raise ValueError('Batch-1 full MF Wilson differs from reconstructed saved source')
        record['stage']='constructing_strict_wilson'
        strict,book=construct_strict_eft1(case,hardc5)
        if _max_rel(wz['C_psi2phi2_ms'],strict['C_psi2phi2_ms_strict_1loop'])>0.20:
            # Do not impose proximity as a hard gate: non-small shifts might be
            # genuine and are recorded in the summary below instead.
            record['large_resummed_vs_strict_difference']=True
        # Decompose named quartics directly from the source’s saved native
        # RGBeta coordinates in the same selected tensor convention as Batch 1.
        q_at_mf=np.asarray(case.background.quartic[-1],float)
        q_at_ms=np.asarray(bz['quartic_real_full_ms'],float)
        namesf,qf,audit_f=named_quartics(root,model,'MF',q_at_mf)
        names_s,qs,audit_s=named_quartics(root,model,'MS',q_at_ms)
        hard_parts,hardaudit=split_existing_f_hard(root,model,batch1,namesf,qf)
        record['stage']='saving_perturbative_inputs'
        folder=out/'inputs'/slug
        folder.mkdir(parents=True,exist_ok=True)
        outw=dict(wz)
        # Strict C requires the *order-zero* full Wilson source for the bubble,
        # not the historical Batch-1 array that included the O(hbar) hard C5.
        outw['C_psi2phi2_tree_mf']=strict['C_psi2phi2_tree0_mf']
        outw['C_psi2phi2_ms_strict_1loop']=strict['C_psi2phi2_ms_strict_1loop']
        outw['C5_EFT1_MS_strict_1loop']=strict['C5_EFT1_MS_strict']
        outwp=folder/'wilson_at_MS.npz'
        put_npz(outwp,outw)
        outwm=dict(wm)
        outwm.update({'tensor_npz_sha256':hash_file(outwp),
                      'strict_1loop_endpoint_at_ms':True,
                      'strict_1loop_origin':'tree0_MF + sourced_F_triangle_hard1_MF + ln(MS/MF)/(16pi2)*full_candidate_beta(tree0_MF)',
                      'strict_1loop_source':'Numerical.parallel_t3.wilson.strict_one_loop_batch2.construct_strict_eft1',
                      'strict_f_hard_loop_inserted_once':True,
                      'strict_MF_tree_in_wilson_npz_excludes_hard_C5':True,
                      'strict_full_operator_tensor_not_only_C5':True,
                      'strict_ADM_candidate':'fonseca_eq22_six_term',
                      'one_loop_physical_prediction_certified':False,
                      'one_loop_missing_hard_and_scalar_topologies':True})
        wout=folder/'wilson.json'
        put_json(wout,outwm)
        outb=dict(bz)
        outb['quartic_real_channels_ms']=qs
        outbp=folder/'background.npz'
        put_npz(outbp,outb)
        outbm=dict(bm)
        outbm.update({'tensor_npz_sha256':hash_file(outbp),
                      'quartic_channel_names':list(names_s),
                      'named_source_quartic_decomposition_at_MS':audit_s,
                      'full_quartic_same_as_Batch1':True,
                      'physical_quartic_normalisation_certified':False})
        bout=folder/'background.json'
        put_json(bout,outbm)
        put_npz(folder/'strict_wilson_channels.npz',strict)
        put_json(folder/'strict_wilson_channels.json',book)
        put_npz(folder/'MF_hard_triangle_by_quartic.npz',hard_parts)
        put_json(folder/'quartic_channel_provenance.json',{
            'MF_quartic_source':audit_f,'MS_quartic_source':audit_s,
            'MF_hard_source_closure':hardaudit,
            'MF_channel_order':list(namesf),'MS_channel_order':list(names_s),
            'hard_triangle_omitted_topologies_not_supplied':True,
            'missing_extra_quartic_source_directions_not_silently_injected':True})
        record['source_inputs']={
            'batch1_summary_sha256':hash_file(batch1/'summary.json'),
            'batch1_wilson_npz_sha256':hash_file(prefix/'wilson_at_MS.npz'),
            'batch1_scalar_npz_sha256':hash_file(prefix/'background.npz'),
            'batch1_hard_npz_sha256':hash_file(hardpath)}
        record['strict_EFT1']=book
        record['MF_hard_quartic_channels']=list(namesf)
        record['MS_scalar_quartic_channels']=list(names_s)
        record['MF_hard_C5']=cmat(hardc5)
        record['MF_hard_split']=hardaudit
        record['strict_C5_at_MS']=cmat(strict['C5_EFT1_MS_strict'])
        record['EFT1_resummed_vs_strict_C5_relative']=_max_rel(
            wz['C5_EFT1_ms'],strict['C5_EFT1_MS_strict'])
        record['stage']='C_preflight'
        pf=preflight_pair(wout,bout,mode='strict_one_loop',mu_low=mu_low,expected_model=model)
        record['preflight_status']=pf.get('status')
        if pf.get('status')!='ready_conditional':
            raise ValueError('Strict Workstream C preflight blocked: '+str(pf.get('blockers')))
        if preflight_only:
            record['status']='PREPARED_STRICT_ONE_LOOP'
            record['stage']='prepared_only'
            return record
        record['stage']='C_joint_scalar_matching_and_SMEFT'
        results=out/'results'/slug
        results.parent.mkdir(parents=True,exist_ok=True)
        payload=run_case(wilson_manifest=wout,background_manifest=bout,
                         output_dir=results,mode='strict_one_loop',mu_low=mu_low,
                         vev=vev,allow_projection_residual=False,
                         match_only=False,overwrite=False)
        if payload.get('SMEFT_running',{}).get('success') is not True:
            raise ValueError('Existing SM+Weinberg one-loop numerical ODE not successful')
        low=read_npz(results/'C5_smeft_low.npz')
        mat=read_npz(results/'C5_smeft_at_MS.npz')
        obs,matrices=neutrino_observables(low['m_nu_GeV'],outb['ye'])
        if not np.allclose(obs['masses_eV'],low['mass_singular_values_eV'],rtol=1e-8,atol=1e-13):
            raise ValueError('Low-energy Takagi masses disagree with existing SMEFT output')
        # The strict-mode C result must be the exact sum of hard-MF boundary,
        # EFT1 first-order mixing, and the simultaneous scalar bubble at MS.
        expected=(strict['C5_EFT1_MS_strict']+mat['C5_scalar_bubble_MS']+
                  mat['C5_independent_finite_MS'])
        err=_max_rel(mat['C5_SMEFT_MS'],expected)
        if err>1e-9:
            raise ValueError('MF hard + EFT1 mixing + scalar matching decomposition did not close')
        if np.linalg.norm(mat['C5_independent_finite_MS'])>1e-24:
            raise ValueError('Unexpected independent finite hard term at MS')
        # Ensure existing scalar-matching code genuinely created a channel for
        # each source-named quartic (including zero-producing channels).
        channels=read_npz(results/'C5_per_quartic_channel.npz')
        if len(channels)!=len(names_s):
            raise ValueError('Existing Workstream C did not export all named scalar quartic channels')
        split_sum=sum(channels.values(),np.zeros((3,3),complex))
        if _max_rel(split_sum,mat['C5_scalar_bubble_MS'])>1e-8:
            raise ValueError('MS scalar bubble named quartics do not sum to C5 bubble')
        # Every nonzero/zero source-named quartic keeps its own hard-MF and
        # scalar-MS C5 response, and the seven EFT1 ADM terms are separate.
        soft_parts={name:channels[f'channel_{i:03d}'] for i,name in enumerate(names_s)}
        beta_parts={name:strict['C5_EFT1_strict_running_'+name]
                    for name in ('yukawa_wavefunction','scalar_pair',
                                 'mixed_yukawa_gauge','crossed_yukawa',
                                 'conjugate_coefficient_yukawa',
                                 'scalar_anomalous_dimension',
                                 'fermion_anomalous_dimension')}
        matrix_ledger={
            'MF_fermion_hard_triangle_by_quartic':{n:cmat(v) for n,v in hard_parts.items()},
            'EFT1_running_by_ADM_channel':{n:cmat(v) for n,v in beta_parts.items()},
            'MS_scalar_bubble_by_quartic':{n:cmat(v) for n,v in soft_parts.items()},
            'full_MF_hard':cmat(hardc5),
            'full_EFT1_running':cmat(strict['C5_EFT1_strict_running_delta']),
            'full_MS_scalar_bubble':cmat(mat['C5_scalar_bubble_MS']),
            'full_SMEFT_C5_at_MS':cmat(mat['C5_SMEFT_MS']),
            'remaining_uncomputed_hard_MS_finite':None,
            'hard_and_scalar_matching_full_one_loop_complete':False,
            'strict_quartic_running_freeze_point':'MF for EFT1 beta, MS for scalar bubble'}
        put_json(results/'C5_one_loop_channel_ledger.json',matrix_ledger)
        # The metadata ledger is a complete decomposition of the conditional
        # calculation, NOT a claim of a complete UV matching graph inventory.
        # Do not rewrite Workstream C outputs: preserve their own manifest SHA.
        put_npz(results/'PMNS_matrices.npz',matrices)
        obs.update({'model':model,'status':STATUS,'conditional_not_physical':True,
                    'matching_scope':'F-hard triangle + full candidate EFT1 ADM + S-pair bubble',
                    'mu_f_gev':float(case.mf),'mu_s_gev':float(case.ms),
                    'mu_low_gev':mu_low,'vev_gev':vev})
        put_json(results/'observables.json',obs)
        record.update({
            'status':STATUS,'stage':'complete',
            'masses_eV_CONDITIONAL':obs['masses_eV'],
            'dm21_sq_eV2_CONDITIONAL':obs['dm21_sq_eV2'],
            'dm31_sq_eV2_CONDITIONAL':obs['dm31_sq_eV2'],
            'mixing_defined':obs['mixing_defined'],
            'strict_C5_scalar_bubble_MS':cmat(mat['C5_scalar_bubble_MS']),
            'strict_C5_SMEFT_at_MS':cmat(mat['C5_SMEFT_MS']),
            'strict_C5_additive_closure_relative':err,
            'strict_MF_hard_fraction':float(np.linalg.norm(hardc5)/max(np.linalg.norm(mat['C5_SMEFT_MS']),1e-100)),
            'strict_MS_bubble_source':'order_zero_MF_LLSS_both_scalars_common_MS',
            'SMEFT_RGE_implemented':True,
            'finite_one_loop_other_topologies_still_not_computed':True,
        })
    except Exception as exc:
        record['status']='BLOCKED'
        record['blocker']=type(exc).__name__+': '+str(exc)
        print(model,'BLOCKED',record['stage'],record['blocker'],file=sys.stderr)
    return record


def _mass_csv(out:Path,batch1_summary:dict,summary:dict)->None:
    rows=[]
    for model in ALL_MODELS:
        item=summary['results'].get(model,{})
        if item.get('status')!=STATUS:
            continue
        baseline=batch1_summary['results'].get(model,{})
        b=baseline.get('masses_eV_CONDITIONAL')
        now=item['masses_eV_CONDITIONAL']
        row={'model':model}
        for i in range(3):
            row[f'm{i+1}_strict_eV']=now[i]
            row[f'm{i+1}_resummed_batch1_eV']=b[i] if b else ''
            row[f'm{i+1}_difference_eV']=now[i]-b[i] if b else ''
        row['dm21_sq_strict_eV2']=item['dm21_sq_eV2_CONDITIONAL']
        row['dm31_sq_strict_eV2']=item['dm31_sq_eV2_CONDITIONAL']
        row['hard_C5_norm_gev_inv']=item['MF_hard_C5']['frobenius_gev_inverse']
        row['scalar_bubble_C5_norm_gev_inv']=item['strict_C5_scalar_bubble_MS']['frobenius_gev_inverse']
        row['EFT1_C5_strict_vs_resummed_rel']=item['EFT1_resummed_vs_strict_C5_relative']
        rows.append(row)
    if rows:
        path=out/'strict_vs_resummed_all16.csv'
        with path.open('w',newline='',encoding='utf-8') as file:
            w=csv.DictWriter(file,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def run_batch(root:Path,*,batch1:Path=BASE,output:Path=OUTPUT,
              models=ALL_MODELS,resume:bool=False,overwrite:bool=False,
              preflight_only:bool=False,mu_low:float=173.,vev:float=246.22)->dict:
    root=Path(root).resolve()
    batch1=Path(batch1);out=Path(output)
    if not batch1.is_absolute():batch1=root/batch1
    if not out.is_absolute():out=root/out
    if out.resolve()==batch1.resolve() or out.resolve()==root:
        raise ValueError('Batch 2 output must be a separate directory')
    if overwrite and resume:
        raise ValueError('Select only --resume or --overwrite')
    source=_verify_first_batch(batch1,models)
    if out.exists() and any(out.iterdir()):
        if not resume and not overwrite:
            raise FileExistsError('Existing Batch 2 output: pass --resume or --overwrite')
        old=read_json(out/'summary.json')
        if old.get('format')!=SCHEMA:
            raise ValueError('Cannot overwrite or resume unrecognized output')
        if overwrite:
            shutil.rmtree(out)
            summary=None
        else:
            summary=old
    else:summary=None
    out.mkdir(parents=True,exist_ok=True)
    if summary is None:
        summary={'format':SCHEMA,'results':{},'mf_gev':100_000.,'ms_gev':1_000.,
                 'mu_low_gev':mu_low,'vev_gev':vev,
                 'EFT_order':'UV -> ALL F out -> EFT1 with BOTH S -> ALL S out -> SMEFT',
                 'strict_one_loop_MF_to_MS':True,
                 'RGE_improved_SMEFT_below_MS':True,
                 'full_eft1_weinberg_operator_mixing':True,
                 'MF_hard_triangle_included':True,
                 'MF_other_hard_topologies_not_computed':True,
                 'MS_scalar_pair_bubble_included':True,
                 'MS_other_scalar_finite_topologies_not_computed':True,
                 'full_ADM_not_certified':True,
                 'physical_masses_certified':False,
                 'baseline_batch1_summary_sha256':hash_file(batch1/'summary.json'),
                 'scope':SCOPE}
    else:
        if (summary.get('baseline_batch1_summary_sha256')!=hash_file(batch1/'summary.json') or
                summary.get('mu_low_gev')!=mu_low or summary.get('vev_gev')!=vev):
            raise ValueError('Resuming with changed benchmark sources or SM parameters')
    skipped=(STATUS,'PREPARED_STRICT_ONE_LOOP') if preflight_only else (STATUS,)
    selected=[model for model in models if not resume or
              summary['results'].get(model,{}).get('status') not in skipped]
    if resume:print('Resume selected:',', '.join(selected) if selected else '(none)')
    put_json(out/'summary.json',summary)
    for model in selected:
        for root_dir in ('inputs','results'):
            path=out/root_dir/tag(model)
            if path.is_dir():shutil.rmtree(path)
        record=run_model(root,batch1,out,model,mu_low=mu_low,vev=vev,
                         preflight_only=preflight_only)
        summary['results'][model]=record
        put_json(out/'summary.json',summary)
        if record['status']==STATUS:
            print(model,record['status'],'masses_eV:',record['masses_eV_CONDITIONAL'])
        elif record['status']=='PREPARED_STRICT_ONE_LOOP':
            print(model,'STRICT WILSON AND C PREFLIGHT READY')
    summary['complete_count']=sum(r.get('status')==STATUS for r in summary['results'].values())
    summary['blocked_count']=sum(r.get('status')=='BLOCKED' for r in summary['results'].values())
    summary['prepared_count']=sum(r.get('status')=='PREPARED_STRICT_ONE_LOOP' for r in summary['results'].values())
    summary['executed_this_call']=selected
    _mass_csv(out,source,summary)
    put_json(out/'summary.json',summary)
    return summary


def show(path:Path)->None:
    summary=read_json(path/'summary.json')
    if summary.get('format')!=SCHEMA:
        raise ValueError('Invalid Batch 2 summary')
    print('Model    Status                                     Strict conditional masses (eV)')
    print('-'*118)
    for model in ALL_MODELS:
        record=summary.get('results',{}).get(model,{})
        masses=record.get('masses_eV_CONDITIONAL')
        label=', '.join(f'{v:.7g}' for v in masses) if masses else record.get('blocker','not run')
        print(f'{model:<8} {record.get("status","NOT_RUN"):<42} {label}')
    print('Strict O(hbar) for included topologies; not complete physical one-loop matching.')


def main(argv=None)->int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root',type=Path,default=ROOT)
    parser.add_argument('--batch1',type=Path,default=BASE)
    parser.add_argument('--output',type=Path,default=OUTPUT)
    parser.add_argument('--resume',action='store_true')
    parser.add_argument('--overwrite',action='store_true')
    parser.add_argument('--preflight-only',action='store_true')
    parser.add_argument('--show',action='store_true')
    parser.add_argument('--model',choices=ALL_MODELS)
    parser.add_argument('--mu-low',type=float,default=173.)
    parser.add_argument('--vev',type=float,default=246.22)
    args=parser.parse_args(argv)
    out=args.output if args.output.is_absolute() else args.project_root/args.output
    try:
        if args.show:
            show(out)
            return 0
        models=[args.model] if args.model else ALL_MODELS
        summary=run_batch(args.project_root,batch1=args.batch1,
                          output=args.output,models=models,resume=args.resume,
                          overwrite=args.overwrite,preflight_only=args.preflight_only,
                          mu_low=args.mu_low,vev=args.vev)
        print('Batch summary:',out/'summary.json','complete=',summary['complete_count'],
              'blocked=',summary['blocked_count'])
        return 0 if summary['blocked_count']==0 else 2
    except Exception as exc:
        print('ALL16 STRICT BATCH2 ERROR:',type(exc).__name__,str(exc),file=sys.stderr)
        return 1

if __name__=='__main__':
    raise SystemExit(main())
