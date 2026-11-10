"""T3 Batch 3: all-16 one-loop heavy-scalar Coleman-Weinberg Higgs threshold.

Works from the completed Batch 2 *strict* numerical sources. Numerically
matches the SM Higgs quadratic and quartic couplings when BOTH scalar
multiplets are removed at the same MS. This is genuine one-loop renormalisable
potential matching within the declared real-field convention, NOT an additional
Weinberg hard diagram. Optional matched-lambdaH SMEFT re-evolution illustrates
higher-order feedback on the already-one-loop C5 and is explicitly flagged.

Also inventories ALL saved native F-only LLHH Laurent kernel groups so that
non-renormalized Matchete terms are not silently treated as missing zeros.
They are NOT added to the Batch 1 triangle (which could double count).
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

from Numerical.parallel_t3.wilson.basis import ALL_MODELS
from Numerical.parallel_t3.background.contracts import tag
from Numerical.parallel_t3.threshold.scalar_cw_batch3 import one_loop_higgs_potential
from Numerical.parallel_t3.threshold.native_source_inventory_batch3 import source_inventory
from Numerical.parallel_t3.threshold.pipeline import evolve_existing_smeft
from Numerical.parallel_t3.threshold.rough_numerics import neutrino_observables

SCHEMA='t3-scalar-cw-one-loop-all16-batch3-v1'
STATUS='NUMERICAL_SCALAR_CW_ONE_LOOP_CANDIDATE_NOT_CERTIFIED'
BATCH2=Path('output/parallel_t3/strict_one_loop_all16_batch2_v1')
OUTPUT=Path('output/parallel_t3/scalar_cw_all16_batch3_v1')


def digest(path:Path)->str:
    return sha256(Path(path).read_bytes()).hexdigest()


def load_json(path:Path)->dict:
    d=json.loads(Path(path).read_text(encoding='utf-8-sig'))
    if not isinstance(d,dict):raise ValueError('Expected JSON object at '+str(path))
    return d


def load_npz(path:Path)->dict:
    with np.load(path,allow_pickle=False) as z:
        return {k:np.asarray(z[k]).copy() for k in z.files}


def save_json(path:Path,data:dict)->None:
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix('.partial.json')
    tmp.write_text(json.dumps(data,indent=2,sort_keys=True,allow_nan=False)+'\n',encoding='utf-8')
    tmp.replace(path)


def save_npz(path:Path,**arrays)->None:
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix('.partial.npz')
    np.savez_compressed(tmp,**arrays)
    tmp.replace(path)


def complex_summary(x:np.ndarray)->dict:
    x=np.asarray(x,complex)
    return {'real':x.real.tolist(),'imag':x.imag.tolist(),'frobenius_norm':float(np.linalg.norm(x))}


def one_model(root:Path,batch2:Path,out:Path,model:str,*,mu_low:float,vev:float,
              transfer_factor:float,skip_low_feedback:bool)->dict:
    slug=tag(model)
    item={'model':model,'slug':slug,'status':'BLOCKED','stage':'loading'}
    try:
        base=batch2/'inputs'/slug
        bmeta=load_json(base/'background.json')
        source=load_npz(base/'background.npz')
        if (bmeta.get('model')!=model or bmeta.get('tensor_npz_sha256')!=digest(base/'background.npz') or
                bmeta.get('single_common_scalar_threshold') is not True):
            raise ValueError('Saved Batch-2 scalar source hash/model/common-threshold invalid')
        if bmeta.get('quartic_potential_convention')!='V4_lambda_real_symmetric_over_4factorial_v1':
            raise ValueError('Changed real-scalar potential normalization')
        atms=batch2/'results'/slug/'C5_smeft_at_MS.npz'
        atlow=batch2/'results'/slug/'C5_smeft_low.npz'
        ms=load_npz(atms)
        low=load_npz(atlow)
        if ms['C5_SMEFT_MS'].shape!=(3,3):
            raise ValueError('Missing original strict SMEFT Weinberg coefficient')
        names=tuple(bmeta.get('quartic_channel_names',()))
        if not names or 'quartic_real_channels_ms' not in source:
            raise ValueError('Batch 2 did not retain named scalar quartic channels')
        item['stage']='Coleman_Weinberg_one_loop_matching'
        cw=one_loop_higgs_potential(
            source['mass_squared_gev2_ms'],source['quartic_real_full_ms'],1000.,
            channels=source['quartic_real_channels_ms'],names=names)
        # These invariants should reduce to gauge singlets, since the source
        # matrices already passed the all-16 mass/quartic Ward gate.
        if cw['quadratic_singlet_residual']>1e-7 or cw['quartic_singlet_residual']>1e-7:
            raise ValueError('One-loop Higgs threshold is not SU2-invariant: '
                             f"{cw['quadratic_singlet_residual']}, {cw['quartic_singlet_residual']}")
        bsm=bmeta.get('sm_couplings_at_ms',{})
        if not isinstance(bsm,dict) or 'lambdaH' not in bsm:
            raise ValueError('Source saved SM Higgs quartic at MS missing')
        before=float(bsm['lambdaH'])
        delta=float(transfer_factor)*cw['delta_lambdaH_source_convention']
        after=before+delta
        if not np.isfinite(after):raise ValueError('Nonfinite threshold-matched SM quartic')
        item['stage']='native_F_source_Laurent_inventory'
        bfolder=root/'output/parallel_t3/background/B_all16_qualitative_v1'/slug
        bsource=load_json(bfolder/'background_manifest.json')
        inv=source_inventory(root,model,slug,set(bsource['beta_symbol_order_eft1']))
        folder=out/'results'/slug
        folder.mkdir(parents=True,exist_ok=True)
        save_json(folder/'native_F_Laurent_not_injected.json',inv)
        save_npz(folder/'scalar_CW_one_loop_MS.npz',
            heavy_masses_squared_gev2=cw['mass2_gev2'],
            heavy_mass_rotation=cw['mass_eigenvectors'],
            delta_mH2_real_4x4_gev2=cw['higgs_mass2_shift_gev2'],
            delta_mH2_by_channel_gev2=cw['higgs_mass2_channel_shifts_gev2'],
            delta_lambda4_H_real_4x4x4x4=cw['higgs_quartic_shift_tensor'],
            delta_lambda4_H_by_ordered_channel_pair=cw['higgs_quartic_channel_pairs'],
            delta_lambdaH_by_ordered_channel_pair=cw['quartic_projected_channel_pairs'],
            C5_SMEFT_MS_unchanged=ms['C5_SMEFT_MS'])
        if not np.allclose(cw['quartic_projected_channel_pairs'].sum(),
                           cw['delta_lambdaH_source_convention'],rtol=1e-10,atol=1e-15):
            raise ValueError('Quartic source-pair decomposition does not close')
        item['CW_threshold']={
            'muS_gev':1000.,
            'delta_mH2_source_gev2':cw['delta_mH_squared_source_convention_gev2'],
            'delta_lambdaH_source':cw['delta_lambdaH_source_convention'],
            'quadratic_singlet_residual':cw['quadratic_singlet_residual'],
            'quartic_singlet_residual':cw['quartic_singlet_residual'],
            'num_real_heavy_scalar_modes':cw['num_heavy_real_modes'],
            'named_portal_channels':list(names),
            'lambdaH_before':before,
            'lambdaH_after_candidate':after,
            'lambdaH_source_to_SM_beta_transfer_factor':transfer_factor,
            'mass2_gev2':cw['mass2_gev2'].tolist(),
            'all_scalar_modes_integrated_out_together':True,
            'scalar_CW_is_one_loop_renormalisable_matching':True,
            'Weinberg_C5_matching_diagram_added':False,
            'other_scalar_gauge_and_wf_matching_not_computed':True}
        item['native_F_Laurent']={
            'PL_kernel_count':inv['num_native_PL_kernels'],
            'PR_kernel_count':inv['num_native_PR_kernels'],
            'native_quartic_symbols':inv['native_quartics'],
            'source_file_sha256':inv['source_sha256'],
            'raw_one_loop_terms_added':False,
            'EFT1_soft_and_counterterm_subtraction_still_missing':True}
        item['input_source_hashes']={'Batch2_scalar_input_npz':digest(base/'background.npz'),
                                     'Batch2_strict_C5_at_MS':digest(atms),
                                     'Batch2_low':digest(atlow)}
        item['unchanged_strict_C5_SMEFT_MS']=complex_summary(ms['C5_SMEFT_MS'])
        item['stage']='optional_SM_matching_feedback'
        if not skip_low_feedback:
            # Multiplying the already O(hbar) C5 with a lambdaH O(hbar)
            # shift through its one-loop RGE is O(hbar^2).  This is a
            # *higher-order exploratory comparison*, never called strict O(hbar).
            matched_meta=dict(bmeta)
            matched_meta['sm_couplings_at_ms']={**bsm,'lambdaH':after}
            evolved=evolve_existing_smeft(ms['C5_SMEFT_MS'],source,matched_meta,
                                          mu_s=1000.,mu_low=mu_low)
            if not evolved['solver_success']:
                raise ValueError('Existing SM+Weinberg solver failed after scalar CW threshold')
            mass=-.5*vev*vev*evolved['c5_low']
            observables,pmns=neutrino_observables(mass,source['ye'])
            oldmass=np.asarray(low['mass_singular_values_eV'],float)
            newmass=np.asarray(observables['masses_eV'],float)
            save_npz(folder/'CW_feedback_SM_Weinberg_low.npz',
                C5_low_with_CW_lambdaH=evolved['c5_low'],
                m_nu_low_GeV_with_CW_lambdaH=mass,
                masses_eV_with_CW_lambdaH=newmass,
                mass_difference_from_strict_Batch2_eV=newmass-oldmass,**pmns)
            save_json(folder/'CW_feedback_observables.json',observables)
            item['masses_eV_Batch2_strict_baseline']=oldmass.tolist()
            item['masses_eV_CW_feedback_higher_order']=newmass.tolist()
            item['delta_masses_eV_CW_feedback_higher_order']=(newmass-oldmass).tolist()
            item['CW_feedback_is_strict_one_loop_in_C5']=False
            item['SMEFT_feedback_ODE_ran']=True
        else:
            item['masses_eV_Batch2_strict_baseline']=np.asarray(low['mass_singular_values_eV'],float).tolist()
            item['SMEFT_feedback_ODE_ran']=False
        item['stage']='complete'
        item['status']=STATUS
        save_json(folder/'model_manifest.json',item)
    except Exception as exc:
        item['blocker']=type(exc).__name__+': '+str(exc)
        print(model,'BLOCKED',item['stage'],item['blocker'],file=sys.stderr)
    return item


def batch(root:Path,from_batch2:Path,output:Path,models,*,resume:bool,overwrite:bool,
          mu_low:float,vev:float,transfer_factor:float,skip_low_feedback:bool):
    root=Path(root).resolve()
    base=from_batch2 if from_batch2.is_absolute() else root/from_batch2
    out=output if output.is_absolute() else root/output
    if base.resolve()==out.resolve() or out.resolve()==root:
        raise ValueError('Output must be separate from the completed strict benchmark')
    if resume and overwrite:raise ValueError('Cannot request both --resume and --overwrite')
    prev=load_json(base/'summary.json')
    if prev.get('format')!='t3-strict-one-loop-all16-hard-triangle-and-scalar-bubble-v1':
        raise ValueError('Not the expected completed Batch2 numerical input')
    for model in models:
        if prev.get('results',{}).get(model,{}).get('status')!='NUMERICAL_STRICT_ONE_LOOP_INCLUDED_TOPOLOGIES_CANDIDATE_NOT_CERTIFIED':
            raise ValueError('Batch2 incomplete for '+model)
    old=None
    if out.exists() and any(out.iterdir()):
        if not resume and not overwrite:
            raise FileExistsError('Output already exists; pass --resume or --overwrite')
        old=load_json(out/'summary.json')
        if old.get('format')!=SCHEMA:raise ValueError('Unrecognized existing Batch3 output')
        if overwrite:
            shutil.rmtree(out)
            old=None
    out.mkdir(parents=True,exist_ok=True)
    basehash=digest(base/'summary.json')
    if old is None:
        summary={'format':SCHEMA,'results':{},'MF_gev':100000.,'MS_gev':1000.,
                 'mu_low_gev':mu_low,'vev_gev':vev,
                 'lambdaH_transfer_factor':transfer_factor,
                 'SMEFT_feedback_enabled':not skip_low_feedback,
                 'batch2_summary_sha256':basehash,
                 'EFT_order':'UV -> ALL F out -> EFT1(S1,S2) -> BOTH S out at MS -> SMEFT',
                 'scalar_one_loop_Higgs_potential_matching_included':True,
                 'native_Laurent_C5_terms_not_injected':True,
                 'previous_Batch2_strict_C5_unchanged_at_MS':True,
                 'full_physical_matching_certified':False}
    else:
        summary=old
        if (summary['batch2_summary_sha256']!=basehash or
            summary['mu_low_gev']!=mu_low or summary['vev_gev']!=vev or
            summary['lambdaH_transfer_factor']!=transfer_factor or
            summary['SMEFT_feedback_enabled']==skip_low_feedback):
            raise ValueError('Resume source/parameter mismatch')
    selected=[m for m in models if not resume or summary['results'].get(m,{}).get('status')!=STATUS]
    if resume:print('Resume selected:',', '.join(selected) if selected else '(none)')
    save_json(out/'summary.json',summary)
    for model in selected:
        path=out/'results'/tag(model)
        if path.exists():shutil.rmtree(path)
        result=one_model(root,base,out,model,mu_low=mu_low,vev=vev,
                         transfer_factor=transfer_factor,skip_low_feedback=skip_low_feedback)
        summary['results'][model]=result
        save_json(out/'summary.json',summary)
        if result['status']==STATUS:
            m=result['CW_threshold']
            print(model,STATUS,'delta_lambdaH:',f"{m['delta_lambdaH_source']:.9g}",
                  'delta_mH2_GeV2:',f"{m['delta_mH2_source_gev2']:.9g}")
    summary['complete_count']=sum(x.get('status')==STATUS for x in summary['results'].values())
    summary['blocked_count']=sum(x.get('status')=='BLOCKED' for x in summary['results'].values())
    save_json(out/'summary.json',summary)
    rows=[]
    for model in ALL_MODELS:
        x=summary['results'].get(model,{})
        if x.get('status')!=STATUS:continue
        q=x['CW_threshold']
        row={'model':model,'delta_lambdaH_MS':q['delta_lambdaH_source'],
             'delta_mH2_GeV2':q['delta_mH2_source_gev2'],
             'quadratic_singlet_residual':q['quadratic_singlet_residual'],
             'quartic_singlet_residual':q['quartic_singlet_residual'],
             'native_F_PL_Laurent_groups':x['native_F_Laurent']['PL_kernel_count']}
        if x.get('masses_eV_CW_feedback_higher_order'):
            for j in range(3):
                row[f'm{j+1}_baseline_strict_eV']=x['masses_eV_Batch2_strict_baseline'][j]
                row[f'm{j+1}_CW_feedback_eV']=x['masses_eV_CW_feedback_higher_order'][j]
                row[f'dm{j+1}_CW_feedback_eV']=x['delta_masses_eV_CW_feedback_higher_order'][j]
        rows.append(row)
    if rows:
        headers=list(dict.fromkeys(k for row in rows for k in row))
        with (out/'scalar_CW_all16.csv').open('w',newline='',encoding='utf-8') as fp:
            writer=csv.DictWriter(fp,fieldnames=headers);writer.writeheader();writer.writerows(rows)
    return summary


def show(out:Path)->None:
    s=load_json(out/'summary.json')
    if s.get('format')!=SCHEMA:raise ValueError('Incorrect Batch3 result schema')
    print('Model   Result                                   delta lambdaH       delta mH2 (GeV2)   native F PL')
    print('-'*120)
    for m in ALL_MODELS:
        x=s['results'].get(m,{})
        if x.get('status')==STATUS:
            q=x['CW_threshold']
            print(f"{m:<7} {x['status']:<40} {q['delta_lambdaH_source']:+.9e}  {q['delta_mH2_source_gev2']:+.9e}  {x['native_F_Laurent']['PL_kernel_count']:>3}")
        else:print(m,x.get('status','NOT_RUN'),x.get('blocker',''))
    print('Renormalizable CW potential threshold only; C5 unchanged at strict O(hbar).')


def main(argv=None)->int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root',type=Path,default=ROOT)
    p.add_argument('--batch2',type=Path,default=BATCH2)
    p.add_argument('--output',type=Path,default=OUTPUT)
    p.add_argument('--model',choices=ALL_MODELS)
    p.add_argument('--resume',action='store_true')
    p.add_argument('--overwrite',action='store_true')
    p.add_argument('--show',action='store_true')
    p.add_argument('--skip-low-feedback',action='store_true')
    p.add_argument('--lambdaH-transfer-factor',type=float,default=1.)
    p.add_argument('--mu-low',type=float,default=173.)
    p.add_argument('--vev',type=float,default=246.22)
    args=p.parse_args(argv)
    out=args.output if args.output.is_absolute() else args.project_root/args.output
    try:
        if args.show:
            show(out)
            return 0
        models=[args.model] if args.model else ALL_MODELS
        data=batch(args.project_root,args.batch2,args.output,models,
                   resume=args.resume,overwrite=args.overwrite,
                   mu_low=args.mu_low,vev=args.vev,
                   transfer_factor=args.lambdaH_transfer_factor,
                   skip_low_feedback=args.skip_low_feedback)
        print('Batch summary:',out/'summary.json','complete=',data['complete_count'],'blocked=',data['blocked_count'])
        return 0 if data['blocked_count']==0 else 2
    except Exception as exc:
        print('BATCH3 ERROR',type(exc).__name__,str(exc),file=sys.stderr)
        return 1

if __name__=='__main__':raise SystemExit(main())
