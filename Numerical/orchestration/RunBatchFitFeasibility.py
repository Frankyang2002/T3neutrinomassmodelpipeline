"""Run the existing single-case feasibility pilot over every discovered T3 model.

Separate results per model, resumable, preserving model-specific RGE/matching.
"""
from __future__ import annotations
import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

from Numerical.orchestration.ResumePipelineNumericalResults import (
    PROJECT_ROOT, OUTPUT_ROOT, discover_cases,
)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scenario', default='smallY_smallL')
    p.add_argument('--points', type=int, default=256)
    p.add_argument('--local-evals', type=int, default=30)
    p.add_argument('--limit', type=int, default=0, help='Maximum models; 0 means all')
    p.add_argument('--rerun', action='store_true', help='Repeat completed cases')
    p.add_argument('--run', action='store_true', help='Execute numerical fits')
    args = p.parse_args(argv)
    cases = sorted((key,case) for scenario,key,case in discover_cases(OUTPUT_ROOT)
                   if scenario == args.scenario)
    if args.limit:
        cases = cases[:args.limit]
    if not cases:
        p.error('No comparison cases found for scenario '+args.scenario)
    root = PROJECT_ROOT/'Reports'/'output'/'full'/'fit_feasibility'/args.scenario
    print(f'Scenario: {args.scenario}; models: {len(cases)}; output: {root}', flush=True)
    print('Class E: original per-representation artifacts and evaluator are reused; no universal multiplicity factor applied.',flush=True)
    if not args.run:
        for key,_ in cases:
            print('  ',key)
        print('Dry run. Add --run to fit.')
        return 0
    rows=[]
    for index,(key,case) in enumerate(cases,1):
        path=root/key/'feasibility_summary.json'
        error=''
        if path.is_file() and not args.rerun:
            state='AlreadyComplete'
        else:
            print(f'[{index}/{len(cases)}] {key}',flush=True)
            command=[sys.executable,'-m','Numerical.orchestration.RunFitFeasibility',
                     '--scenario',args.scenario,'--model',key,
                     '--points',str(args.points),'--local-evals',str(args.local_evals),'--run']
            proc=subprocess.run(command,cwd=PROJECT_ROOT,check=False)
            state='Success' if proc.returncode==0 and path.is_file() else 'Failed'
            if state=='Failed':
                error=f'Exit code {proc.returncode}; see console output'
        row={'scenario':args.scenario,'model':key,'class':key.split('_')[1] if '_' in key else '',
             'status':state,'chi2':'','best_scale':'','local_status':'','error':error}
        if path.is_file():
            try:
                summary=json.loads(path.read_text(encoding='utf-8'))
                best=summary.get('best',{})
                row['chi2']=best.get('chi2','')
                row['best_scale']=summary.get('best_scale',{}).get('factor','')
                row['local_status']=summary.get('local',{}).get('status','')
                row.update({k: v for k,v in best.get('prediction',{}).items()})
            except Exception as exc:
                row['status']='Failed'
                row['error']=f'Invalid summary: {exc}'
        rows.append(row)
        root.mkdir(parents=True,exist_ok=True)
        fields=list(dict.fromkeys(k for item in rows for k in item))
        with (root/'batch_summary.csv').open('w',newline='',encoding='utf-8') as handle:
            writer=csv.DictWriter(handle,fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
        (root/'batch_summary.json').write_text(json.dumps(rows,indent=2)+'\n',encoding='utf-8')
        print(f'  {row["status"]}: chi2={row["chi2"]}',flush=True)
    print('Saved:',root/'batch_summary.csv')
    return 1 if any(row['status']=='Failed' for row in rows) else 0


if __name__=='__main__':
    raise SystemExit(main())
