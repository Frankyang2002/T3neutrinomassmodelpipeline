"""Read-only all16 Ward gate, using actual saved Stream1 and Stream2 endpoints.

Default: inexpensive tensor/gauge tests. --with-beta calls the *existing*
full-tensor RGE beta at MF (no ODE re-running) and traces its seven terms.
"""
from __future__ import annotations
import argparse,csv,json,sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
from Numerical.parallel_t3.threshold.models import ALL_MODELS,normalize_model
from Numerical.parallel_t3.threshold.rough_ward_gate_v6 import audit_model,FORMAT


def flatten(r):
    row={'model':r['model'],'status':r['status']}
    if r['status']!='SourceWardFirstCauseDiagnosticOnly':
        row['blocker']=r.get('blocker');return row
    for scale in ('MF','MS'):
        for sector in ('Wilson','Yukawa','quartic'):
            part=r[scale][sector]
            row[scale+'_'+sector+'_norm']=part['tensor_norm']
            for group in ('SU2','U1','SU3'):
                row[scale+'_'+sector+'_'+group]=part['ward_relative'][group]
    units=r['source_quartic_units']
    row['bad_source_quartic_unit_names']=';'.join(units.get('noninvariant_unit_names',()))
    row['first_cause_candidates']=';'.join(r.get('first_cause_candidates',()))
    if r.get('saved_candidate_beta_MF'):
        for group in ('SU2','U1','SU3'):
            row['MF_saved_beta_'+group]=r['saved_candidate_beta_MF']['ward_relative'][group]
    if r.get('candidate_beta_MF'):
        b=r['candidate_beta_MF']
        row['MF_beta_recomputed_diff']=b['stored_beta_relative_disagreement']
        for group in ('SU2','U1','SU3'):
            row['MF_beta_'+group]=b['total']['ward_relative'][group]
    return row


def main(argv=None):
    p=argparse.ArgumentParser(description='Trace T3 Stream2 gauge violation to MF Wilson, Yukawa, quartic or one-loop beta channel')
    p.add_argument('--root',type=Path,default=Path('.'))
    p.add_argument('--models',nargs='+',default=None)
    p.add_argument('--run-id',default='rough_fermion_first_v1')
    p.add_argument('--with-beta',action='store_true',help='Recompute saved MF beta ONCE using project engine (potentially expensive)')
    p.add_argument('--skip-source-units',action='store_true',help='Skip individual source-coded quartic CG Ward tests')
    p.add_argument('--output',type=Path,default=None)
    a=p.parse_args(argv)
    models=[normalize_model(x) for x in (a.models or ALL_MODELS)]
    if len(models)!=len(set(models)):p.error('Duplicate model')
    destination=a.output or (a.root/'output/parallel_t3/rough_results'/a.run_id/'_diagnostics'/'first_cause_v6')
    destination.mkdir(parents=True,exist_ok=True)
    rows=[]
    for model in models:
        slug=model.replace('_-','_m').replace('_+','_p')
        try:
            result=audit_model(a.root,model,input_run_id=a.run_id,
                               with_beta=a.with_beta,with_units=not a.skip_source_units)
            row=flatten(result)
            print(f'{model:<5} MF(C,Y,Q) SU2= '
                  f'{row["MF_Wilson_SU2"]:.4g}, '
                  f'{row["MF_Yukawa_SU2"]:.4g}, '
                  f'{row["MF_quartic_SU2"]:.4g} ; '
                  f'MS SU2(W)={row["MS_Wilson_SU2"]:.4g} ; '
                  f'bad units={row["bad_source_quartic_unit_names"] or "none"} ; '
                  f'flags={row["first_cause_candidates"] or "none"}',flush=True)
        except Exception as exc:
            result={'format':FORMAT,'model':model,'status':'SOURCE_AUDIT_BLOCKED',
                    'blocker':f'{type(exc).__name__}: {exc}'}
            row=flatten(result)
            print(model,'SOURCE_AUDIT_BLOCKED',result['blocker'],flush=True)
        (destination/f'{slug}_first_cause.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        rows.append(row)
    cols=list(dict.fromkeys(key for row in rows for key in row))
    with (destination/'first_cause_summary.csv').open('w',encoding='utf-8',newline='') as fh:
        writer=csv.DictWriter(fh,fieldnames=cols);writer.writeheader();writer.writerows(rows)
    (destination/'first_cause_summary.json').write_text(json.dumps({
        'format':'t3-rough-first-cause-summary-v6',
        'status':'SourceWardFirstCauseDiagnosticOnly','models':rows,
        'count':len(rows),'blocked':sum(r['status']=='SOURCE_AUDIT_BLOCKED' for r in rows),
        'matching_or_RGE_not_recomputed':True,'read_only':True
    },indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('FIRST CAUSE:',destination/'first_cause_summary.csv',flush=True)
    return 1 if any(r['status']=='SOURCE_AUDIT_BLOCKED' for r in rows) else 0

if __name__=='__main__':sys.exit(main())
