#!/usr/bin/env python3
"""All-16 explicit EFT1 scalar MSbar CT, native pole and matching-scale sweep.

Adds, never overwrites, a new workstream output. Uses immutable completed
Batch 1/2/4/5 saved tensors, native Laurent kernels and Workstream B quartic
channels. All original matching/running engines remain unchanged.

Strict O(hbar) scan makes no new threshold: mu=0.5/1/2 MF is the SUBTRACTION
scale of the same simultaneous F matching and the MF couplings stay frozen.
The resulting scan is diagnostic; the non-triangle native UV 1PI amplitude,
UV counterterms and signed GammaCC factor are not derived here.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import shutil
import sys
from typing import Mapping

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path: sys.path.insert(0, str(ROOT))
import numpy as np

from Numerical.parallel_t3.background.contracts import tag
from Numerical.parallel_t3.wilson.basis import ALL_MODELS
from Numerical.parallel_t3.wilson.source_audited_all16 import hash_file
from Numerical.parallel_t3.wilson.quartic_channels_batch2 import named_quartics
from Numerical.parallel_t3.wilson.msbar_counterterm_batch6 import (
    SCHEMA, SCENARIOS, ADM_CHANNELS, SCALE_FACTORS, relative, symmetric3,
    source_scalar_counterterm, candidate_adm_counterterms,
    native_pole_by_quartic, analytic_ct_from_soft, matching_scale_scenarios)
from Numerical.parallel_t3.threshold.pipeline import evolve_existing_smeft
from Numerical.parallel_t3.threshold.rough_numerics import neutrino_observables

STATUS = 'NUMERICAL_EXPLICIT_EFT1_SCALAR_MSbar_CT_AND_NATIVE_SCALE_SCAN_CONDITIONAL'
BASE1 = Path('output/parallel_t3/fhard_triangle_all16_batch1_v1')
BASE2 = Path('output/parallel_t3/strict_one_loop_all16_batch2_v1')
BASE4 = Path('output/parallel_t3/native_finite_scenario_all16_batch4_v1')
BASE5 = Path('output/parallel_t3/hard_soft_matching_all16_batch5_v1')
DEST = Path('output/parallel_t3/eft1_msbar_counterterms_all16_batch6_v1')
MF, MS = 1e5, 1e3


def read_json(p: Path) -> dict:
    val = json.loads(Path(p).read_text('utf-8-sig'))
    if not isinstance(val,dict): raise ValueError('Expected JSON object: '+str(p))
    return val


def read_npz(p: Path) -> dict[str, np.ndarray]:
    with np.load(p, allow_pickle=False) as z:
        return {k: np.asarray(z[k]).copy() for k in z.files}


def save_json(p: Path, d: dict) -> None:
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(p.suffix+'.partial')
    tmp.write_text(json.dumps(d, indent=2, sort_keys=True, allow_nan=False)+'\n', encoding='utf-8')
    tmp.replace(p)


def save_npz(p: Path, z: Mapping[str, np.ndarray]) -> None:
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp=p.with_suffix('.partial.npz')
    np.savez_compressed(tmp, **z)
    tmp.replace(p)


def cm(z: np.ndarray) -> dict:
    x=symmetric3(z, 'manifest C5')
    return {'real': x.real.tolist(), 'imag': x.imag.tolist(),
            'frobenius_GeV_inverse': float(np.linalg.norm(x))}


def _source_filenames(root: Path, slug: str) -> dict[str, Path]:
    return {
      'B1_hard': root/BASE1/'hard_fermion'/slug/'fermion_hard_triangle.npz',
      'B1_hard_json': root/BASE1/'hard_fermion'/slug/'fermion_hard_triangle.json',
      'B1_MF_quartic': root/'output/parallel_t3/rough_inputs/fhard_triangle_all16_batch1_v1/a'/slug/'background.npz',
      'B2_strict': root/BASE2/'inputs'/slug/'strict_wilson_channels.npz',
      'B2_strict_json': root/BASE2/'inputs'/slug/'strict_wilson_channels.json',
      'B2_MS_C5': root/BASE2/'results'/slug/'C5_smeft_at_MS.npz',
      'B4_native': root/BASE4/'inputs'/slug/'native_MF_LLHH_finite_source_candidate.npz',
      'B4_native_json': root/BASE4/'inputs'/slug/'native_MF_LLHH_finite_source_candidate.json',
      'B4_MS_C5': root/BASE4/'results'/slug/'C5_smeft_at_MS.npz',
      'B5_pole': root/BASE5/'inputs'/slug/'hard_soft_and_native_Laurent.npz',
      'B5_pole_json': root/BASE5/'inputs'/slug/'hard_soft_and_native_Laurent.json',
      'B5_MS_C5': root/BASE5/'results'/slug/'C5_smeft_at_MS.npz',
      'B5_background': root/BASE5/'inputs'/slug/'background.npz',
      'B5_background_json': root/BASE5/'inputs'/slug/'background.json',
    }


def _inputs(root: Path, model: str) -> dict:
    slug = tag(model)
    paths = _source_filenames(root, slug)
    missing = [name for name,p in paths.items() if not p.is_file()]
    if missing:
        raise FileNotFoundError('Missing completed Batch 1/2/4/5 source files: '+repr(missing))
    h={name: hash_file(p) for name,p in paths.items()}
    hard=read_npz(paths['B1_hard'])
    hm=read_json(paths['B1_hard_json'])
    if hm.get('model')!=model or hm.get('hard_npz_sha256')!=h['B1_hard']:
        raise ValueError('Batch-1 hard triangle hash/model mismatch')
    strict=read_npz(paths['B2_strict'])
    sm=read_json(paths['B2_strict_json'])
    native=read_npz(paths['B4_native'])
    nm=read_json(paths['B4_native_json'])
    pole=read_npz(paths['B5_pole'])
    pm=read_json(paths['B5_pole_json'])
    if (nm.get('model')!=model or pm.get('model')!=model or
            sm.get('mode')!='strict_one_loop' or pm.get('full_one_loop_matching_completed') is not False):
        raise ValueError('Source topology/strict/native provenance changed')
    if (abs(float(sm['mf_gev'])-MF)>1e-5 or
            abs(float(sm['ms_gev'])-MS)>1e-5):
        raise ValueError('Require one F threshold 1e5 and one scalar threshold 1e3 GeV')
    if relative(hard['C5_hard_F_triangle_MF'],strict['C5_MF_hard_fermion_triangle'])>1e-9:
        raise ValueError('Batch-1 and Batch-2 hard MF contributions disagree')
    if relative(native['C5_native_finite_MF'],pole['C5_native_formal_pole_drop_MF'])>1e-9:
        raise ValueError('Batch-4 and Batch-5 native pole-dropped finite disagree')
    if relative(pole['C5_analytic_local_hard_MF'],hard['C5_hard_F_triangle_MF'])>1e-9:
        raise ValueError('Batch-5 analytic hard kernel does not reproduce Batch-1')
    pms=read_npz(paths['B5_MS_C5'])
    bms=read_npz(paths['B5_background'])
    bmeta=read_json(paths['B5_background_json'])
    if bmeta.get('model')!=model or bmeta.get('tensor_npz_sha256')!=h['B5_background']:
        raise ValueError('Batch-5 SMEFT background source hash/model mismatch')
    # First-order matching uses exactly the same source C0 and scalar bubble.
    rb=read_npz(paths['B2_MS_C5'])
    rn=read_npz(paths['B4_MS_C5'])
    for name,z in [('Batch2',rb),('Batch4',rn),('Batch5',pms)]:
        if not all(k in z for k in ('C5_SMEFT_MS','C5_scalar_bubble_MS')):
            raise ValueError(name+' lacks complete strict scalar-matching endpoint')
        if relative(z['C5_scalar_bubble_MS'],rb['C5_scalar_bubble_MS'])>1e-9:
            raise ValueError(name+' scalar bubble changed between strict branches')
    return {'h':h, 'hard':hard,'strict':strict,'strict_meta':sm,'native':native,
            'native_meta':nm,'pole':pole,'pole_meta':pm,
            'B2_MS':rb,'B4_MS':rn,'B5_MS':pms,
            'background_arrays':bms,'background_meta':bmeta,
            'MF_quartic':read_npz(paths['B1_MF_quartic']),
            'all_paths':paths}


def _scan_soft_samples(pole:dict) -> dict[float,np.ndarray]:
    return {0.5: pole['EFT1_soft_finite_mu_m05'],
            1.0: pole['EFT1_soft_finite_mu_m10'],
            2.0: pole['EFT1_soft_finite_mu_m20']}


def run_one(root:Path, out:Path, model:str, *, mu_low:float, vev:float,
            skip_running:bool=False)->dict:
    row={'model':model,'status':'BLOCKED','stage':'loading_saved_batches'}
    try:
        src=_inputs(root,model)
        q_original=src['MF_quartic']
        # Same canonical MF quartic as the original all-16 calculation.
        scales=np.asarray(q_original['mu_gev'],float)
        indices=np.flatnonzero(np.isclose(scales,MF,atol=1e-6,rtol=1e-11))
        if len(indices)!=1: raise ValueError('No unique source MF quartic background')
        qfull=np.asarray(q_original['lambda4'][indices[0]],float)
        names,channels,source_info=named_quartics(root,model,'MF',qfull)
        q=channels[:, :4, :4, 4:, 4:]
        hard=src['hard'];strict=src['strict'];pole=src['pole'];native=src['native']
        if relative(q.sum(axis=0),hard['HHSS_quartic_MF'])>1e-9:
            raise ValueError('All MF quartic source channels do not reconstruct original hard loop')
        row['stage']='construct_explicit_EFT1_MSbar_pole_counterterm'
        ct=source_scalar_counterterm(hard['C_LLSS_per_fermion_canonical'],q,names)
        native_k=native['C5_native_by_kernel_and_heavy']
        masses=np.asarray(hard['heavy_masses_gev'],float)
        factor=1+2*np.log(MF/masses)
        source_constants=native_k/factor[np.newaxis,:,np.newaxis,np.newaxis]
        if relative(source_constants.sum(axis=(0,1)),pole['C5_native_raw_pole_residue'])>1e-9:
            raise ValueError('Native raw source pole/residue disagrees with Batch 5')
        records=src['native_meta']['kernel_records']
        if len(records)!=len(source_constants):
            raise ValueError('Native saved kernel records mismatch source tensor')
        labels=[entry['native_quartic'] for entry in records]
        nsource=native_pole_by_quartic(source_constants, labels)
        ct_other=candidate_adm_counterterms(strict,float(src['strict_meta']['log_ms_over_mf']))
        ct_extracted=-symmetric3(pole['C5_EFT1_scalar_loop_pole_residue'], 'Batch5 soft pole')
        ct_integral=analytic_ct_from_soft(_scan_soft_samples(pole),SCALE_FACTORS)
        if (relative(ct['CT_scalar_LLHH'],ct_extracted)>1e-8 or
                relative(ct['CT_scalar_LLHH'],ct_integral)>1e-8):
            raise ValueError('Explicit MSbar B0 counterterm differs from independent saved loop integral')
        # The selected candidate ADM may use a different normalization. Record
        # that disagreement; do NOT fit an arbitrary counterterm to cancel it.
        rg_discrepancy=relative(ct['CT_scalar_LLHH'],ct_other['CT_candidate_scalar_pair'])
        native_discrepancy=relative(nsource_sum(nsource),ct['CT_scalar_LLHH'])
        row['stage']='propagate_subtraction_scale_through_strict_EFT1_and_SMEFT'
        scan=matching_scale_scenarios(
            native_constant_by_kernel_fermion=source_constants,masses_gev=masses,
            c5_analytic_hard_mf=hard['C5_hard_F_triangle_MF'],
            c5_native_formal_mf=native['C5_native_finite_MF'],
            c5_running_delta_mf_to_ms=strict['C5_EFT1_strict_running_delta'],
            c5_scalar_bubble_ms=src['B5_MS']['C5_scalar_bubble_MS'],
            soft_finite=_scan_soft_samples(pole),mu_f=MF,mu_s=MS,
            factors=SCALE_FACTORS)
        refs={'analytic_triangle':src['B2_MS']['C5_SMEFT_MS'],
              'native_formal':src['B4_MS']['C5_SMEFT_MS'],
              'native_as_UV_minus_soft_hypothesis':src['B5_MS']['C5_SMEFT_MS']}
        nominal_reproduction={name:relative(scan['scales'][1.][name]['C5_SMEFT_MS'],ref)
                              for name,ref in refs.items()}
        if max(nominal_reproduction.values())>2e-8:
            raise ValueError('Nominal matching scale does not reproduce completed Batches 2,4,5: '
                             +repr(nominal_reproduction))
        base=out/tag(model)
        matrix_scan={}
        diagnostics={}
        all_rows=[]
        for factor,states in scan['scales'].items():
            diagnostics[str(factor)]={}
            for scenario,z in states.items():
                key=f'{scenario}_scale_{factor:g}'
                matrix_scan[key+'_C5_at_MS']=z['C5_SMEFT_MS']
                record={'model':model,'scale_factor':factor,'subtraction_mu_gev':MF*factor,
                        'scenario':scenario,'C5_at_MS_norm':float(np.linalg.norm(z['C5_SMEFT_MS'])),
                        'not_physical_MSbar':True,'couplings_frozen_at_muF':True}
                if not skip_running:
                    sol=evolve_existing_smeft(z['C5_SMEFT_MS'],
                         src['background_arrays'],src['background_meta'],
                         mu_s=MS,mu_low=mu_low)
                    if sol.get('solver_success') is not True:
                        raise RuntimeError('Existing one-loop SM/Weinberg running failed')
                    c5=symmetric3(sol['c5_low'],key+' low C5')
                    mass=-0.5*vev**2*c5
                    obs,_=neutrino_observables(mass,src['background_arrays']['ye'])
                    record.update({'masses_eV':obs['masses_eV'],
                      'dm21_sq_eV2':obs['dm21_sq_eV2'],
                      'dm31_sq_eV2':obs['dm31_sq_eV2'],
                      'mixing_defined':obs['mixing_defined'],
                      'sin2_theta12':obs['sin2_theta12'],
                      'sin2_theta13':obs['sin2_theta13'],
                      'sin2_theta23':obs['sin2_theta23'],
                      'SMEFT_nfev':int(sol['nfev'])})
                    matrix_scan[key+'_C5_at_low']=c5
                    all_rows.append(record)
                diagnostics[str(factor)][scenario]=record
        row['stage']='write_counterterm_results'
        counter_npz=base/'counterterms_and_native_sources.npz'
        save_npz(counter_npz,{
           'C5_CT_per_quartic_per_fermion':ct['CT_per_channel_per_heavy'],
           'C5_CT_per_quartic':ct['CT_per_channel'],
           'C5_CT_per_heavy_fermion':ct['CT_per_heavy'],
           'C5_CT_EFT1_scalar_bubble_1over_epsbar':ct['CT_scalar_LLHH'],
           'C5_EFT1_bare_scalar_soft_loop_pole':ct['EFT1_soft_bare_pole'],
           'C5_EFT1_renormalized_soft_loop_plus_CT_pole':ct['renormalized_soft_pole'],
           'C5_native_raw_pole':nsource_sum(nsource),
           'C5_native_minus_explicit_EFT1_CT_residual':nsource_sum(nsource)-ct['CT_scalar_LLHH'],
           'C5_ADM_candidate_seven_channel_CT':ct_other['CT_candidate_total_7channel'],
           'C5_ADM_candidate_scalar_CT':ct_other['CT_candidate_scalar_pair'],
           **{'C5_candidate_MSbar_CT_'+n:mat for n,mat in ct_other['CT_per_ADM_channel'].items()},
           **{'C5_native_raw_pole_'+n:mat for n,mat in nsource.items()}})
        scan_npz=base/'matching_scale_scenario_matrices.npz'
        save_npz(scan_npz,matrix_scan)
        counter_summary={
           'schema':SCHEMA,'model':model,'status':STATUS,
           'EFT_order':['UV','ALL heavy fermions out at one MF','EFT1 with BOTH S1 and S2',
                        'BOTH scalars out at one MS','SMEFT Weinberg running'],
           'mu_f_gev':MF,'mu_s_gev':MS,'mu_low_gev':mu_low,'vev_gev':vev,
           'source_SHA256':src['h'],
           'EFT1_scalar_loop_bare_pole_source':'C_LLSS_tree * lambda_HHSS_full /(2*16pi2)',
           'EFT1_MSbar_counterterm_calculated':True,
           'EFT1_scalar_bubble_pole_cancelled_by_its_own_CT':True,
           'EFT1_scalar_bubble_CT_matrix':cm(ct['CT_scalar_LLHH']),
           'EFT1_scalar_pole_projection_residual':ct['projected_tensor_purity_residual'],
           'original_triangle_B0_CT_relative':relative(ct['CT_scalar_LLHH'],ct_integral),
           'original_Batch5_soft_pole_CT_relative':relative(ct['CT_scalar_LLHH'],ct_extracted),
           'candidate_ADM_scalar_pair_CT_relative':rg_discrepancy,
           'native_raw_pole_vs_explicit_EFT1_CT_relative':native_discrepancy,
           'native_raw_pole_matrix':cm(nsource_sum(nsource)),
           'native_unassigned_pole_matrix':cm(nsource_sum(nsource)-ct['CT_scalar_LLHH']),
           'CT_per_RGBeta_quartic_source_names':{name:cm(value) for name,value in zip(names,ct['CT_per_channel'])},
           'native_pole_per_Matchete_quartic_source_name':{name:cm(value) for name,value in nsource.items()},
           'candidate_MSbar_CT_7_ADM_channels':{name:cm(value) for name,value in ct_other['CT_per_ADM_channel'].items()},
           'one_common_subtraction_mu_variations':[MF*f for f in SCALE_FACTORS],
           'nominal_mu_f_reproductions_Batch2_4_5':nominal_reproduction,
           'renormalizable_scalar_quartic_input_scope':source_info,
           'status_of_native':('Saved Matchete one-loop MatchedEFT expressions are NOT established '
                              'as bare UV amplitudes; their pole is NOT cancelled by fitting a CT'),
           'native_full_UV_counterterms_calculated':False,
           'complete_one_loop_UV_EFT1_matching':False,
           'signed_GammaCC_normalization_certified':False,
           'mass_results_physical':False,
           'scalar_threshold_existing_strict_source_unchanged':True,
           'one_loop_matching_scale_variation_interpretation':(
               'Fixed UV couplings at physical MF; varied subtraction scale only. '
               'Does not match a new UV trajectory at a different fermion threshold.'),
           'counterterm_npz_sha256':hash_file(counter_npz),
           'scale_npz_sha256':hash_file(scan_npz),
           'scale_scan':diagnostics}
        save_json(base/'counterterm_and_scale_ledger.json',counter_summary)
        if all_rows:
            columns=list(all_rows[0])
            csvpath=base/'mass_and_mixing_scan.csv'
            tmp=csvpath.with_suffix('.partial.csv')
            with tmp.open('w',newline='',encoding='utf-8') as stream:
                writer=csv.DictWriter(stream,fieldnames=columns)
                writer.writeheader();writer.writerows(all_rows)
            tmp.replace(csvpath)
        return {'model':model,'status':STATUS,'stage':'complete',
                'EFT1_scalar_CT_residual_B0':counter_summary['original_triangle_B0_CT_relative'],
                'scalar_ADM_CT_relative':rg_discrepancy,
                'native_pole_unassigned_relative':native_discrepancy,
                'nominal_masses_eV_by_scenario':{
                  s:diagnostics['1.0'][s].get('masses_eV') for s in SCENARIOS},
                'source_quartic_channels':len(names),
                'native_PL_kernels':len(records),
                'source_hash':hash_file(base/'counterterm_and_scale_ledger.json'),
                'output':str((base/'counterterm_and_scale_ledger.json').relative_to(root))}
    except Exception as exc:
        row['blocker']=type(exc).__name__+': '+str(exc)
        print(model, 'BLOCKED', row['stage'], row['blocker'],file=sys.stderr)
        return row


def nsource_sum(d:Mapping[str,np.ndarray])->np.ndarray:
    return symmetric3(sum(d.values(),np.zeros((3,3),complex)), 'native original pole')


def write_summary(out:Path,results:dict,executed:list,*,skip_running:bool,mu_low:float,vev:float)->dict:
    completed=sum(z.get('status')==STATUS for z in results.values())
    manifest={'schema':SCHEMA,'status':'CONDITIONAL_EFT1_SCALAR_MSbar_CT_CALCULATED_NATIVE_UV_STILL_UNRESOLVED',
              'completed_count':completed,'blocked_count':16-completed,'results':results,
              'executed_in_this_call':executed,
              'EFT_order':'UV -> ALL F out -> EFT1 BOTH S active -> BOTH S out -> SMEFT',
              'MSbar_scalar_loop_CT_calculated':True,
              'native_UV_counterterms_calculated':False,
              'mass_scan_run':not skip_running,
              'mu_low_gev':mu_low,'vev_gev':vev,
              'full_one_loop_matching_complete':False,
              'neutrino_masses_physical':False}
    save_json(out/'summary.json',manifest)
    tables=[]
    for model in ALL_MODELS:
        r=results.get(model,{})
        if r.get('status') != STATUS:continue
        for scenario,mass in r['nominal_masses_eV_by_scenario'].items():
            tables.append({'model':model,'scenario':scenario,
                'm1_eV':mass[0] if mass is not None else None,
                'm2_eV':mass[1] if mass is not None else None,
                'm3_eV':mass[2] if mass is not None else None,
                'EFT1_CT_scalar_relative':r['EFT1_scalar_CT_residual_B0'],
                'candidate_ADM_CT_relative':r['scalar_ADM_CT_relative'],
                'native_unassigned_pole_relative':r['native_pole_unassigned_relative']})
    if tables:
        tmp=out/'model_scenario_comparison.partial.csv'
        with tmp.open('w',newline='',encoding='utf-8') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(tables[0]))
            writer.writeheader();writer.writerows(tables)
        tmp.replace(out/'model_scenario_comparison.csv')
    return manifest


def main(argv=None)->int:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root',type=Path,default=ROOT)
    ap.add_argument('--out',type=Path,default=DEST)
    ap.add_argument('--model',choices=ALL_MODELS)
    ap.add_argument('--resume',action='store_true')
    ap.add_argument('--overwrite',action='store_true')
    ap.add_argument('--show',action='store_true')
    ap.add_argument('--skip-low-energy-running',action='store_true')
    ap.add_argument('--mu-low',type=float,default=173.)
    ap.add_argument('--vev',type=float,default=246.22)
    args=ap.parse_args(argv)
    root=args.root.resolve();out=(root/args.out).resolve() if not args.out.is_absolute() else args.out.resolve()
    if out==root or not out.is_relative_to(root/'output'):
        raise ValueError('Counterterm outputs must remain inside the local project output tree')
    if not 0<args.mu_low<=MS or not 0<args.vev<1000:
        raise ValueError('Invalid low scale or Higgs VEV')
    summary_path=out/'summary.json'
    old=read_json(summary_path) if summary_path.exists() else None
    if old and old.get('schema')!=SCHEMA:
        raise ValueError('Batch6 summary schema changed; refuse to overwrite unrelated output')
    if args.show:
        if old is None:
            print('Batch6 output summary not found:',summary_path);return 1
        for model in ALL_MODELS:
            d=old['results'].get(model,{})
            print(f'{model:6} {d.get("status","MISSING"):65} '
                  f'CT/B0={d.get("EFT1_scalar_CT_residual_B0","-")} '
                  f'native-pole-residual={d.get("native_pole_unassigned_relative","-")} '
                  f'nominal_analytic_masses={d.get("nominal_masses_eV_by_scenario",{}).get("analytic_triangle")}')
        print('complete',old['completed_count'],'blocked',old['blocked_count']);return 0
    if args.overwrite and args.resume:raise ValueError('Select overwrite OR resume, not both')
    if args.overwrite and out.exists():shutil.rmtree(out)
    if out.exists() and any(out.iterdir()) and not args.resume and not args.overwrite:
        raise FileExistsError('Batch6 outputs exist. Use --resume or --overwrite')
    out.mkdir(parents=True,exist_ok=True)
    results=dict(old.get('results',{})) if old and not args.overwrite else {}
    requested=[args.model] if args.model else list(ALL_MODELS)
    selected=[]
    for model in requested:
        row=results.get(model,{})
        output=out/tag(model)/'counterterm_and_scale_ledger.json'
        good=(row.get('status')==STATUS and output.is_file()
              and row.get('source_hash')==hash_file(output))
        if args.resume and good: continue
        selected.append(model)
    print('Selected:', ', '.join(selected) if selected else 'none (already completed)')
    for model in selected:
        row=run_one(root,out,model,mu_low=args.mu_low,vev=args.vev,
                    skip_running=args.skip_low_energy_running)
        results[model]=row
        manifest=write_summary(out,results,selected,
                               skip_running=args.skip_low_energy_running,
                               mu_low=args.mu_low,vev=args.vev)
        if row.get('status')==STATUS:
            print(model,STATUS,'EFT1_CT_vs_B0',f'{row["EFT1_scalar_CT_residual_B0"]:.3g}',
                  'native_pole_residual',f'{row["native_pole_unassigned_relative"]:.3g}',
                  'masses_eV',row['nominal_masses_eV_by_scenario']['analytic_triangle'])
    manifest=write_summary(out,results,selected,
                           skip_running=args.skip_low_energy_running,
                           mu_low=args.mu_low,vev=args.vev)
    print('Batch 6 summary:',summary_path,'completed=',manifest['completed_count'],
          'blocked=',manifest['blocked_count'])
    return 0 if all(results.get(model,{}).get('status')==STATUS for model in requested) else 1


if __name__=='__main__':sys.exit(main())
