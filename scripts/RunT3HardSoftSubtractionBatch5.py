"""All-16 UV triangle/EFT1-soft subtraction + conditional native-UV hypothesis.

Reads immutable completed Batches 1, 2, 4. Keeps existing strict
EFT1 running, one common scalar threshold and SMEFT running unchanged.

Two distinct objects MUST NOT be confused:
(1) source-triangle full UV - EFT1 soft: an actual subtraction identity for
    an already implemented topology, with dimension-7+ remainder separated;
(2) native-as-UV minus EFT1 soft: a hypothetical numerical sensitivity branch.
    The saved native Matchete MatchedEFT is not a bare UV amplitude.

No data are written into existing Batch 1-4 output directories.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
import json
import shutil
import sys

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
import numpy as np

from Numerical.parallel_t3.background.contracts import tag
from Numerical.parallel_t3.wilson.basis import ALL_MODELS,build_basis
from Numerical.parallel_t3.wilson.source_audited_all16 import hash_file
from Numerical.parallel_t3.wilson.coordinates import project_c5
from Numerical.parallel_t3.wilson.matching_soft_subtraction_batch5 import (
    SCHEMA as LEDGER_SCHEMA,triangle_soft_ledger,native_laurent_from_evaluated,
    comparison_with_native,_rel)
from scripts.RunT3NativeFiniteBatch4 import replace_strict_hard
from Numerical.parallel_t3.wilson.native_finite_batch4 import NATIVE_SCHEMA
from Numerical.parallel_t3.threshold.preflight import preflight_pair
from Numerical.parallel_t3.threshold.pipeline import run_case
from Numerical.parallel_t3.threshold.rough_numerics import neutrino_observables

SCHEMA='t3-native-UV-hypothesis-and-analytic-hard-soft-all16-batch5-v1'
STATUS='NUMERICAL_NATIVE_AS_UV_MINUS_EFT1_SOFT_HYPOTHESIS_NOT_PHYSICAL'
PREPARED='PREPARED_ANALYTIC_HARD_SOFT_AND_NATIVE_POLE_LEDGER'
BASE1=Path('output/parallel_t3/fhard_triangle_all16_batch1_v1')
BASE2=Path('output/parallel_t3/strict_one_loop_all16_batch2_v1')
BASE4=Path('output/parallel_t3/native_finite_scenario_all16_batch4_v1')
OUTPUT=Path('output/parallel_t3/hard_soft_matching_all16_batch5_v1')


def read_json(path:Path)->dict:
    obj=json.loads(Path(path).read_text(encoding='utf-8-sig'))
    if not isinstance(obj,dict):raise ValueError('Expected JSON object '+str(path))
    return obj


def read_npz(path:Path)->dict:
    with np.load(path,allow_pickle=False) as f:
        return {name:np.asarray(f[name]).copy() for name in f.files}


def save_json(path:Path,obj:dict)->None:
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.partial')
    tmp.write_text(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+'\n',encoding='utf-8')
    tmp.replace(path)


def save_npz(path:Path,**arrays)->None:
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix('.partial.npz')
    np.savez_compressed(tmp,**arrays)
    tmp.replace(path)


def serial(a:np.ndarray):
    x=np.asarray(a,complex)
    return {'real':x.real.tolist(),'imag':x.imag.tolist(),
            'frobenius_gev_inverse':float(np.linalg.norm(x))}


def _validate_input(root:Path,b1:Path,b2:Path,b4:Path,model:str)->dict:
    slug=tag(model)
    hp=b1/'hard_fermion'/slug/'fermion_hard_triangle.npz'
    hm=read_json(hp.with_suffix('.json'))
    b1z=read_npz(hp)
    if hm.get('model')!=model or hm.get('hard_npz_sha256')!=hash_file(hp):
        raise ValueError('Batch-1 hard triangle source hash/model mismatch')
    src=b2/'inputs'/slug
    bs=read_npz(src/'strict_wilson_channels.npz')
    b2w=read_npz(src/'wilson_at_MS.npz')
    b2m=read_json(src/'wilson.json')
    if b2m.get('tensor_npz_sha256')!=hash_file(src/'wilson_at_MS.npz'):
        raise ValueError('Batch-2 strict Wilson manifest/NPZ checksum mismatch')
    native=b4/'inputs'/slug
    nw=read_npz(native/'wilson_at_MS.npz')
    nm=read_json(native/'wilson.json')
    nsource=native/'native_MF_LLHH_finite_source_candidate.npz'
    nz=read_npz(nsource)
    nmeta=read_json(native/'native_MF_LLHH_finite_source_candidate.json')
    if (nm.get('tensor_npz_sha256')!=hash_file(native/'wilson_at_MS.npz') or
            nmeta.get('schema')!=NATIVE_SCHEMA or nmeta.get('model')!=model):
        raise ValueError('Batch-4 native source format/model or Wilson NPZ hash mismatch')
    if (_rel(nz['C5_native_finite_MF'],nz['C5_native_by_kernel_and_heavy'].sum(axis=(0,1)))>1e-10 or
            _rel(b1z['C5_hard_F_triangle_MF'],bs['C5_MF_hard_fermion_triangle'])>1e-10 or
            _rel(bs['C_psi2phi2_ms_strict_1loop'],b2w['C_psi2phi2_ms_strict_1loop'])>1e-10):
        raise ValueError('Saved native/hard/strict tensor inputs fail algebraic closure')
    if nm.get('strict_1loop_endpoint_at_ms') is not True:
        raise ValueError('Batch-4 strict endpoint not declared')
    if (_rel(nz['C5_analytic_triangle_replaced_MF'],b1z['C5_hard_F_triangle_MF'])>1e-10 or
            _rel(bs['C5_EFT1_strict_running_delta']+nz['C5_native_finite_MF'],
                 nw['C5_EFT1_MS_strict_1loop'])>1e-9):
        raise ValueError('Batch-4 native replacement and strict running do not close')
    if (b1z['C_LLSS_per_fermion_canonical'].shape[0]!=3 or
            len(b1z['heavy_masses_gev'])!=3):
        raise ValueError('All three fermion eigenstates must be present')
    return {'base1':b1z,'base2':bs,'base4':nz,'native_wilson':nw,
            'native_manifest':nm,'native_source':nmeta,
            'native_inputs':native,'source_hashes':{
            'batch1_triangle_npz':hash_file(hp),
            'batch2_strict_npz':hash_file(src/'strict_wilson_channels.npz'),
            'batch4_native_source_npz':hash_file(nsource),
            'batch4_native_source_json':hash_file(native/'native_MF_LLHH_finite_source_candidate.json')}}


def run_one(root:Path,b1:Path,b2:Path,b4:Path,out:Path,model:str,*,
            mu_low:float=173.,vev:float=246.22,audit_only:bool=False)->dict:
    slug=tag(model)
    row={'model':model,'status':'BLOCKED','stage':'input_validation'}
    try:
        data=_validate_input(root,b1,b2,b4,model)
        b=data['base1'];s=data['base2'];n=data['base4'];bw=data['native_wilson']
        mu_f=100_000.
        masses=np.asarray(b['heavy_masses_gev'],float)
        row['stage']='full_UV_minus_EFT1_soft_exact_integral'
        triangle=triangle_soft_ledger(
            b['C_LLSS_per_fermion_canonical'],masses,
            b['HHSS_quartic_MF'],b['mass_squared_heavy_MS'],mu_f)
        p=triangle['scales'][1.]
        if _rel(p['hard_local_dim5'],b['C5_hard_F_triangle_MF'])>1e-9:
            raise ValueError('Original Batch-1 hard MF kernel is not reproduced')
        if _rel(p['full_UV_triangle'],b['C5_full_UV_triangle_diagnostic'])>1e-9:
            raise ValueError('Full massive triangle differs from Batch-1 source')
        row['stage']='native_PL_Laurent_pole_counterterm_ledger'
        native=native_laurent_from_evaluated(n['C5_native_by_kernel_and_heavy'],masses,mu_f)
        if _rel(native['formal_pole_dropped'],n['C5_native_finite_MF'])>1e-11:
            raise ValueError('Native PL source constant/log reassembly differs from Batch 4')
        comparison=comparison_with_native(triangle,native)
        row['stage']='alternative_strict_Wilson_endpoint'
        hypo=comparison['hypothetical_hard_MF']
        basis=build_basis(model)
        full,c5,replace=replace_strict_hard(
            bw['C_psi2phi2_ms_strict_1loop'],n['C5_native_finite_MF'],hypo,basis,
            running_delta=s['C5_EFT1_strict_running_delta'],
            zero_order_full=s['C_psi2phi2_tree0_mf'])
        dest=out/'inputs'/slug
        dest.mkdir(parents=True,exist_ok=True)
        compact_pairs=[{key:value for key,value in entry.items()
                        if key not in ('projected_C5_base_not_independently_gauge_invariant',)}
                       for entry in _pair_reference(b['C_LLSS_per_fermion_canonical'],
                          masses,b['HHSS_quartic_MF'],b['mass_squared_heavy_MS'])]
        save_json(dest/'heavy_scalar_pair_kernels.json',{
            'model':model,'pair_count':len(compact_pairs),'pairs':compact_pairs,
            'one_common_fermion_threshold_gev':mu_f,
            'couplings_and_masses_from_saved_batch1':True})
        arrays={}
        for factor,z in triangle['scales'].items():
            lab={0.5:'m05',1.:'m10',2.:'m20'}[factor]
            for name,val in z.items():
                if isinstance(val,np.ndarray):arrays[f'{name}_mu_{lab}']=val
        arrays.update({
            'C5_native_raw_pole_residue':native['pole_residue'],
            'C5_native_raw_finite_constant':native['finite_constant'],
            'C5_native_raw_finite_scale_log':native['finite_scale_log'],
            'C5_native_formal_pole_drop_MF':native['formal_pole_dropped'],
            'C5_EFT1_scalar_loop_pole_residue':triangle['eft_soft_pole_at_MF'],
            'C5_residual_pole_before_counterterms':comparison['unaccounted_pole_or_counterterm_residual'],
            'C5_hypothetical_native_UV_minus_EFT1_soft_at_MF':hypo,
            'C5_hypothetical_EFT1_at_MS':c5,
            'C5_analytic_local_hard_MF':p['hard_local_dim5'],
            'C5_full_UV_triangle':p['full_UV_triangle'],
            'C5_power_suppressed_dim7plus':p['power_suppressed_dim7plus_remainder'],
            'fermion_masses_GeV':masses,
            'heavy_scalar_mass2_GeV2':b['mass_squared_heavy_MS']})
        ap=dest/'hard_soft_and_native_Laurent.npz'
        save_npz(ap,**arrays)
        wp=dest/'wilson_at_MS.npz'
        w=dict(bw)
        w['C_psi2phi2_ms_strict_1loop']=full
        w['C5_EFT1_MS_strict_1loop']=c5
        save_npz(wp,**w)
        meta=dict(data['native_manifest'])
        meta.update({'tensor_npz_sha256':hash_file(wp),
            'strict_1loop_origin':('tree0_MF + HYPOTHETICAL native Matchete raw UV finite '
                '- analytic EFT1 scalar soft(MF) + strict EFT1 one-loop beta(tree0)'),
            'native_is_UNSUBTRACTED_UV_amplitude_ASSUMED_NOT_VERIFIED':True,
            'native_matchete_actual_source_is_MatchedEFT_not_bare_UV':True,
            'analytic_UV_triangle_exact_mass_dependence_only_diagnostic':True,
            'missing_counterterms_computed':False,
            'full_one_loop_hard_matching_completed':False,
            'physical_MSbar_hard_matching_claimed':False,
            'C5_hypothesis_source':'native_as_UV_minus_EFT1_soft_conditional_only'})
        wmeta=dest/'wilson.json'
        save_json(wmeta,meta)
        for name in ('background.npz','background.json'):
            shutil.copy2(data['native_inputs']/name,dest/name)
        ledger={
            'schema':LEDGER_SCHEMA,'model':model,'status':'AnalyticSubtractionAndHypotheticalNativeUVOnly',
            'source_files_SHA256':data['source_hashes'],
            'one_common_MF_gev':mu_f,'one_common_MS_gev':1_000.,
            'all_three_F_removed_at_MF':True,'both_S_active_between_MF_MS':True,
            'massive_triangle_full_minus_soft_analytic_identity':True,
            'scalar_mass_hessian_from_saved_MS_used_as_IR_regulator':True,
            'matching_mu_over_MF_scans':list(triangle['scales']),
            'matching_scale_cancellation_relative':triangle['matching_scale_cancellation_relative'],
            'full_UV_scale_invariance_relative':triangle['full_UV_scale_invariance_relative'],
            'dimension7plus_relative_to_full_triangle':triangle['dimension7plus_relative_to_full_triangle'],
            'native_source_kernel_count':int(n['C5_native_by_kernel_and_heavy'].shape[0]),
            'source_native_constant_and_log_pole_equality':True,
            'native_pole_cancellation_relative_no_CT':comparison['pole_cancellation_relative'],
            'native_vs_analytic_local_hard_relative':comparison['native_vs_analytic_hard_relative'],
            'physical_MSbar_source_matching_certified':False,
            'native_as_UV_hypothesis_certified':False,
            'C5_raw_native_formal_MF':serial(native['formal_pole_dropped']),
            'C5_hypothetical_source_MF':serial(hypo),
            'C5_EFT1_strict_MS_alternative':serial(c5),
            'C5_EFT1_scalar_soft_finite_MF':serial(p['EFT1_soft_finite']),
            'C5_local_analytic_hard_MF':serial(p['hard_local_dim5']),
            'C5_native_raw_pole_MF':serial(native['pole_residue']),
            'C5_unaccounted_pole_residue':serial(comparison['unaccounted_pole_or_counterterm_residual']),
            'hypothesis_limitation':comparison['hypothesis_not_certified_reason'],
            'strict_Wilson_replacement_audit':replace,
            'input_triangle_projector_residuals':p['projection_residuals'],
            'full_one_loop_matching_completed':False}
        save_json(dest/'hard_soft_and_native_Laurent.json',ledger)
        row.update({'stage':'C_strict_preflight','ledger':str(dest/'hard_soft_and_native_Laurent.json'),
            'triangle_subtraction_scale_residual':ledger['matching_scale_cancellation_relative'],
            'dimension7plus_to_UV_triangle':ledger['dimension7plus_relative_to_full_triangle'],
            'native_pole_plus_soft_pole_relative':ledger['native_pole_cancellation_relative_no_CT'],
            'native_finite_to_analytic_triangle_relative':ledger['native_vs_analytic_local_hard_relative'],
            'C5_hypothetical_MF_norm':float(np.linalg.norm(hypo)),
            'C5_strict_MS_hypothetical_norm':float(np.linalg.norm(c5)),
            'native_PL_kernels':ledger['native_source_kernel_count']})
        pre=preflight_pair(wmeta,dest/'background.json',
            mode='strict_one_loop',mu_low=mu_low,expected_model=model)
        if pre.get('status')!='ready_conditional':
            raise ValueError('Strict C preflight: '+str(pre.get('blockers')))
        if audit_only:
            row.update({'status':PREPARED,'stage':'prepared_only'})
            return row
        row['stage']='one_common_scalar_threshold_and_SMEFT'
        results=out/'results'/slug
        payload=run_case(wilson_manifest=wmeta,
            background_manifest=dest/'background.json',output_dir=results,
            mode='strict_one_loop',mu_low=mu_low,vev=vev,
            match_only=False,overwrite=False,allow_projection_residual=False)
        if payload.get('SMEFT_running',{}).get('success') is not True:
            raise ValueError('Existing SM/Weinberg running did not succeed')
        at_ms=read_npz(results/'C5_smeft_at_MS.npz')
        at_low=read_npz(results/'C5_smeft_low.npz')
        if _rel(at_ms['C5_SMEFT_MS'],c5+at_ms['C5_scalar_bubble_MS']+
                at_ms['C5_independent_finite_MS'])>1e-9:
            raise ValueError('F-hard hypothesis + EFT1 running + S matching does not close')
        if np.linalg.norm(at_ms['C5_independent_finite_MS'])>1e-24:
            raise ValueError('Unrelated C5 finite term unexpectedly included at MS')
        b=read_npz(dest/'background.npz')
        obs,matrices=neutrino_observables(at_low['m_nu_GeV'],b['ye'])
        if not np.allclose(obs['masses_eV'],at_low['mass_singular_values_eV'],rtol=1e-8,atol=1e-13):
            raise ValueError('Existing SMEFT neutrino masses disagree with independent observables')
        save_npz(results/'PMNS_matrices.npz',**matrices)
        obs.update({'status':STATUS,'not_physical_MSbar_matching':True,
            'native_as_UV_assumption_UNVERIFIED':True,
            'no_missing_counterterm_filled_in':True})
        save_json(results/'observables.json',obs)
        row.update({'status':STATUS,'stage':'complete',
            'conditional_masses_eV':obs['masses_eV'],
            'conditional_dm21_sq_eV2':obs['dm21_sq_eV2'],
            'conditional_dm31_sq_eV2':obs['dm31_sq_eV2'],
            'C5_SMEFT_at_MS_norm':float(np.linalg.norm(at_ms['C5_SMEFT_MS']))})
    except Exception as exc:
        row.update({'status':'BLOCKED','blocker':type(exc).__name__+': '+str(exc)})
        print(model,'BLOCKED',row['stage'],row['blocker'],file=sys.stderr)
    return row


def _pair_reference(per,masses,q,mass2):
    """Read one unmodified library result to store the full scalar pair kernel grid."""
    from Numerical.lo_certification.completeness.hard_matching import evaluate_hard_scalar_channel
    r=evaluate_hard_scalar_channel(per,masses,q,mass2,1e5,1e5,
                                   C_total_muF=np.asarray(per).sum(axis=0))
    return r['per_pair_terms']


def export_comparison(out:Path,s:dict,b2:dict,b4:dict)->None:
    rows=[]
    for model in ALL_MODELS:
        x=s['results'].get(model,{})
        if x.get('status')!=STATUS:continue
        two=b2['results'].get(model,{})
        four=b4['results'].get(model,{})
        a=two.get('masses_eV_CONDITIONAL')
        b=four.get('conditional_masses_eV')
        c=x['conditional_masses_eV']
        row={'model':model,'native_PL_kernels':x['native_PL_kernels'],
            'UV_soft_exact_closure':x['triangle_subtraction_scale_residual'],
            'raw_native_pole_residual_before_CT':x['native_pole_plus_soft_pole_relative'],
            'native_to_analytic_hard_relative':x['native_finite_to_analytic_triangle_relative']}
        for i in range(3):
            row[f'm{i+1}_analytic_strict_eV']=a[i] if a else ''
            row[f'm{i+1}_native_poledrop_eV']=b[i] if b else ''
            row[f'm{i+1}_native_as_UV_soft_subtracted_eV']=c[i]
        rows.append(row)
    if rows:
        p=out/'batch5_analytic_vs_native_vs_native_as_UV_all16.csv'
        with p.open('w',newline='',encoding='utf-8') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def run_batch(root:Path,b1:Path,b2:Path,b4:Path,output:Path,models:list[str],*,
              mu_low:float=173.,vev:float=246.22,
              overwrite:bool=False,resume:bool=False,audit_only:bool=False)->dict:
    root=Path(root).resolve()
    a=[x if x.is_absolute() else root/x for x in (b1,b2,b4,output)]
    b1,b2,b4,out=a
    if out.resolve() in (root,b1.resolve(),b2.resolve(),b4.resolve()) or (
         out.resolve() in (x.resolve().parent for x in (b1,b2,b4))):
        raise ValueError('Batch5 output may not overwrite project or a previous batch')
    if overwrite and resume:raise ValueError('Choose --overwrite or --resume')
    sums=[read_json(x/'summary.json') for x in (b1,b2,b4)]
    expected=('t3-fhard-triangle-conditional-all16-batch1-v1',
              't3-strict-one-loop-all16-hard-triangle-and-scalar-bubble-v1',
              't3-native-matchete-pole-drop-conditional-all16-batch4-v1')
    if tuple(x.get('format') for x in sums)!=expected:
        raise ValueError('Require completed original Batch1/Batch2/Batch4 format versions')
    for model in models:
        if (sums[0]['results'].get(model,{}).get('status')!='NUMERICAL_FHARD_TRIANGLE_CANDIDATE_NOT_CERTIFIED' or
            sums[1]['results'].get(model,{}).get('status')!='NUMERICAL_STRICT_ONE_LOOP_INCLUDED_TOPOLOGIES_CANDIDATE_NOT_CERTIFIED' or
            sums[2]['results'].get(model,{}).get('status')!='NUMERICAL_NATIVE_MATCHeTE_FINITE_SCENARIO_NOT_PHYSICAL_MSbar'):
            raise ValueError('Earlier completed source batch missing model '+model)
    summary=None
    if out.exists() and any(out.iterdir()):
        if not overwrite and not resume:raise FileExistsError('Batch5 output exists; use --resume or --overwrite')
        old=read_json(out/'summary.json')
        if old.get('format')!=SCHEMA:raise ValueError('Refuse overwriting unknown output directory')
        old_shas=old.get('source_summary_SHA256',{})
        expected_shas={f'batch{k}':hash_file(p/'summary.json')
                       for k,p in ((1,b1),(2,b2),(4,b4))}
        if old_shas!=expected_shas:
            raise ValueError('Existing Batch5 source hashes differ; choose a NEW --output')
        if old.get('mu_low_gev')!=mu_low or old.get('vev_gev')!=vev:
            raise ValueError('Cannot change low scales/vev under same Batch5 output path')
        if overwrite:shutil.rmtree(out)
        else:summary=old
    out.mkdir(parents=True,exist_ok=True)
    if summary is None:
        summary={'format':SCHEMA,'results':{},
          'source_summary_SHA256':{f'batch{k}':hash_file(p/'summary.json')
                                   for k,p in ((1,b1),(2,b2),(4,b4))},
          'mu_low_gev':mu_low,'vev_gev':vev,
          'one_common_fermion_threshold_GeV':1e5,
          'one_common_scalar_threshold_GeV':1e3,
          'EFT_order':'UV -> ALL F out -> EFT1 BOTH S active -> BOTH S out -> SMEFT',
          'analytic_full_triangle_minus_EFT1_soft_calculated':True,
          'native_as_unsubtracted_UV_assumption_verified':False,
          'UV_counterterms_calculated':False,'missing_F_hard_topologies_calculated':False,
          'strict_full_Wilson_candidate_ADM':True,
          'complete_physical_one_loop_matching':False,
          'neutrino_masses_physical':False}
    requested=list(models)
    if resume:
        skip=(STATUS,PREPARED) if audit_only else (STATUS,)
        requested=[m for m in requested if summary['results'].get(m,{}).get('status') not in skip]
        print('Resume selected:',', '.join(requested) if requested else 'none')
    for model in requested:
        slug=tag(model)
        for rel in ('inputs','results'):
            p=out/rel/slug
            if p.exists():shutil.rmtree(p)
        row=run_one(root,b1,b2,b4,out,model,
                    mu_low=mu_low,vev=vev,audit_only=audit_only)
        summary['results'][model]=row
        if row['status']==STATUS:
            print(model,row['status'],'masses_eV:',row['conditional_masses_eV'],
                  'scale_residual:',f"{row['triangle_subtraction_scale_residual']:.2e}")
        elif row['status']==PREPARED:
            print(model,'ANALYTIC UV-SOFT LEDGER AND NATIVE POLE PREPARED')
        # Save each completed model atomically: --resume survives interrupted batches.
        summary['complete_count']=sum(x.get('status')==STATUS for x in summary['results'].values())
        summary['prepared_count']=sum(x.get('status')==PREPARED for x in summary['results'].values())
        summary['blocked_count']=sum(x.get('status')=='BLOCKED' for x in summary['results'].values())
        summary['executed_this_call']=requested
        save_json(out/'summary.json',summary)
    summary['complete_count']=sum(x.get('status')==STATUS for x in summary['results'].values())
    summary['prepared_count']=sum(x.get('status')==PREPARED for x in summary['results'].values())
    summary['blocked_count']=sum(x.get('status')=='BLOCKED' for x in summary['results'].values())
    summary['executed_this_call']=requested
    if not audit_only:export_comparison(out,summary,sums[1],sums[2])
    save_json(out/'summary.json',summary)
    return summary


def show(output:Path)->None:
    s=read_json(output/'summary.json')
    if s.get('format')!=SCHEMA:raise ValueError('Wrong Batch5 summary')
    print('Model    Status                        Scale cancellation  Raw pole residual   Hypothesis masses (eV)')
    print('-'*120)
    for model in ALL_MODELS:
        z=s.get('results',{}).get(model,{})
        mass=z.get('conditional_masses_eV')
        val=', '.join(f'{v:.8g}' for v in mass) if mass else z.get('blocker','not run')
        print(f"{model:<8} {z.get('status','NOT_RUN'):<29} "
              f"{z.get('triangle_subtraction_scale_residual',float('nan')):>17.3e} "
              f"{z.get('native_pole_plus_soft_pole_relative',float('nan')):>18.3e}   {val}")
    print('Native-as-UV is HYPOTHETICAL, not a certified MSbar subtraction or physical C5.')


def main(argv=None)->int:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--project-root',type=Path,default=ROOT)
    ap.add_argument('--batch1',type=Path,default=BASE1)
    ap.add_argument('--batch2',type=Path,default=BASE2)
    ap.add_argument('--batch4',type=Path,default=BASE4)
    ap.add_argument('--output',type=Path,default=OUTPUT)
    ap.add_argument('--model',choices=ALL_MODELS)
    ap.add_argument('--mu-low',type=float,default=173.)
    ap.add_argument('--vev',type=float,default=246.22)
    ap.add_argument('--overwrite',action='store_true')
    ap.add_argument('--resume',action='store_true')
    ap.add_argument('--audit-only',action='store_true')
    ap.add_argument('--show',action='store_true')
    a=ap.parse_args(argv)
    try:
        out=a.output if a.output.is_absolute() else a.project_root/a.output
        if a.show:show(out);return 0
        s=run_batch(a.project_root,a.batch1,a.batch2,a.batch4,a.output,
                    [a.model] if a.model else list(ALL_MODELS),
                    mu_low=a.mu_low,vev=a.vev,overwrite=a.overwrite,
                    resume=a.resume,audit_only=a.audit_only)
        print('Batch summary:',out/'summary.json','complete=',s['complete_count'],
              'blocked=',s['blocked_count'])
        return 0 if s['blocked_count']==0 else 2
    except Exception as exc:
        print('BATCH5 ERROR:',type(exc).__name__,str(exc),file=sys.stderr)
        return 1

if __name__=='__main__':raise SystemExit(main())
