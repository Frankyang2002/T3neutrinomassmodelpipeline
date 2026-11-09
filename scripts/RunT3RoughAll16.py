"""All-16 fermion-first rough benchmark: existing Stream1 -> Stream2 -> C.

This is a read-only consumer of upstream outputs unless explicit upstream
command templates are provided. It never runs old symbolic source generators.
Missing genuine data or numerical ambiguities produce BLOCKED rows, not numbers.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0,str(ROOT))

from Numerical.parallel_t3.threshold.models import ALL_MODELS,normalize_model
from Numerical.parallel_t3.threshold.rough_adapter import AdapterBlocked,prepare_pair
from Numerical.parallel_t3.threshold.rough_numerics import OMITTED,run_rough_case
from Numerical.parallel_t3.threshold.preflight import preflight_pair

RUN_ID='rough_fermion_first_v1'
B_DEFAULT=Path('output/parallel_t3/rough_inputs')/RUN_ID/'c'
W_DEFAULT=Path('output/parallel_t3/rough_wilson')/RUN_ID
RESULT_DEFAULT=Path('output/parallel_t3/rough_results')/RUN_ID


def slug(model:str)->str:
    return model.replace('_-','_m').replace('_+','_p')


def inputs(root:Path,model:str,kind:str)->tuple[Path,Path]:
    """The requested handoff layout takes priority; no speculative searches."""
    folder=root/slug(model)
    if kind=='background':
        normal=folder/'background.json'
        fallback=folder/'background_manifest.json'
        npz=folder/'background.npz'
    else:
        normal=folder/'wilson.json'
        fallback=folder/'manifest.json'
        npz=folder/'wilson_at_MS.npz'
    return normal if normal.exists() else fallback,npz


def invoke_command(template:str, *, model:str, cwd:Path) -> dict:
    """Explicit optional source stage; no shell, never Wolfram by default."""
    rendered=template.format(model=model,slug=slug(model),project_root=str(cwd))
    command=rendered if os.name=='nt' else shlex.split(rendered)
    result=subprocess.run(command,cwd=cwd,capture_output=True,text=True,
                          errors='replace',check=False)
    report={'command':rendered,'exit_code':result.returncode,
            'stdout_tail':result.stdout[-5000:],'stderr_tail':result.stderr[-5000:]}
    if result.returncode:
        raise RuntimeError(f'Upstream command failed (exit={result.returncode}): '
                           f'{rendered}; {report["stderr_tail"][-700:]}')
    return report


def _save_json(path:Path,obj:dict)->None:
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def _write_summary(root:Path,entries:list[dict],*,allow_synthetic:bool)->None:
    count=sum(e['status']=='RoughNumericalBenchmark' for e in entries)
    smoke=sum(e['status']=='SyntheticFixtureSmokeTest' for e in entries)
    doc={'format':'t3-all16-rough-fermion-first-summary-v1',
         'run_id':RUN_ID,'classification':'RoughNumericalBenchmark_NotCertified',
         'threshold_sequence':['UV','ALL_fermions_integrated_at_one_MF',
                               'EFT1_BOTH_scalars','ALL_scalars_integrated_at_one_MS',
                               'SMEFT_running_to_173_GeV'],
         'common_scalar_threshold_gev':1000.0,
         'mu_low_gev':173.0,
         'counting_mode':'rge_improved',
         'synthetic_smoke_tests_enabled':bool(allow_synthetic),
         'additional_finite_terms':{'default':'explicit_zero_only_for_rough_benchmark',
             'independently_provided_hard_terms':'retained_if_supplied_and_validated',
             'omitted':list(OMITTED)},
         'counts':{'requested':len(entries),'successful_real':count,
                   'successful_synthetic_smoke':smoke,
                   'blocked':sum(e['status']=='BLOCKED' for e in entries)},
         'models':entries}
    _save_json(root/'summary.json',doc)
    fields=['model','model_slug','status','stage','blocker',
            'mu_f_gev','mu_s_gev','mu_low_gev',
            'm1_eV','m2_eV','m3_eV','dm21_sq_eV2','dm31_sq_eV2','dm32_sq_eV2',
            'sin2_theta12','sin2_theta13','sin2_theta23','mixing_defined',
            'scalar_mass_kind','scalar_eigenstates_integrated_out',
            'source_background','source_wilson','output_manifest',
            'matching_conventions','omitted_terms','classification']
    with (root/'summary.csv').open('w',newline='',encoding='utf-8') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields,extrasaction='ignore')
        writer.writeheader()
        for r in entries:
            row=r.copy()
            masses=row.pop('masses_eV',None)
            if masses is not None:
                for i,m in enumerate(masses,1):row[f'm{i}_eV']=m
            row['omitted_terms']='; '.join(OMITTED)
            row['matching_conventions']='C_ijab full sum no factorial; V4=lambda*phi^4/4!; H=(R+iI)/sqrt2; MSbar'
            row['classification']=row['status']
            writer.writerow(row)


def parse_args(argv=None):
    p=argparse.ArgumentParser(description='All16 rough T3 fermion-first scalar threshold and neutrino results')
    p.add_argument('--project-root',type=Path,default=ROOT)
    p.add_argument('--background-root',type=Path,default=B_DEFAULT)
    p.add_argument('--wilson-root',type=Path,default=W_DEFAULT)
    p.add_argument('--output',type=Path,default=RESULT_DEFAULT)
    p.add_argument('--models',nargs='+',default=None,
                   help='Optional subset, default one A_-4 first followed by all 15')
    p.add_argument('--stream1-command',type=str,
                   help='Optional argv command template for missing B files; {model}, {slug}, {project_root}')
    p.add_argument('--stream2-command',type=str,
                   help='Optional argv command template for missing A files; {model}, {slug}, {project_root}')
    p.add_argument('--overwrite',action='store_true',help='Overwrite only C-owned model results')
    p.add_argument('--mu-low',type=float,default=173.)
    p.add_argument('--vev',type=float,default=246.22)
    p.add_argument('--allow-projection-residual',action='store_true',
                   help='Diagnostic benchmark only: project SU2-violating mass-split residual onto C5')
    p.add_argument('--allow-synthetic-smoke-test',action='store_true',
                   help='Tests ONLY: outputs labelled SyntheticFixtureSmokeTest, never counted as genuine')
    return p.parse_args(argv)


def main(argv=None)->int:
    args=parse_args(argv)
    project=args.project_root.resolve()
    def resolve(x:Path)->Path:
        return x if x.is_absolute() else project/x
    b_root,w_root,out=map(resolve,(args.background_root,args.wilson_root,args.output))
    models=[normalize_model(s) for s in (args.models or ALL_MODELS)]
    if len(set(models))!=len(models):
        raise SystemExit('Duplicate model')
    # A(-4) is intentionally attempted before any other model.
    if 'A_-4' in models:models=['A_-4']+[m for m in models if m!='A_-4']
    out.mkdir(parents=True,exist_ok=True)
    rows=[]
    for model in models:
        entry={'model':model,'model_slug':slug(model),'status':'BLOCKED',
               'stage':'inputs','blocker':None,
               'mu_s_gev':1000.0,'mu_low_gev':args.mu_low}
        rows.append(entry)
        b,b_npz=inputs(b_root,model,'background')
        w,w_npz=inputs(w_root,model,'wilson')
        entry.update(source_background=str(b),source_wilson=str(w))
        try:
            entry['stage']='stream1_background'
            if (not b.exists() or not b_npz.exists()) and args.stream1_command:
                entry['stream1_execution']=invoke_command(args.stream1_command,model=model,cwd=project)
                b,b_npz=inputs(b_root,model,'background')
            if not b.exists() or not b_npz.exists():
                raise AdapterBlocked(f'Stream 1 missing {b} or {b_npz}; no synthetic fallback')
            entry['stage']='stream2_wilson'
            if (not w.exists() or not w_npz.exists()) and args.stream2_command:
                entry['stream2_execution']=invoke_command(args.stream2_command,model=model,cwd=project)
                w,w_npz=inputs(w_root,model,'wilson')
            if not w.exists() or not w_npz.exists():
                raise AdapterBlocked(f'Stream 2 missing {w} or {w_npz}; no LLSS->C5 log regenerated')
            entry['source_background'],entry['source_wilson']=str(b),str(w)
            entry['stage']='component_adapter'
            pair=prepare_pair(model,w,b,out/'_adapted_inputs'/slug(model),
                              allow_synthetic=args.allow_synthetic_smoke_test)
            entry['stage']='existing_C_preflight'
            pf=preflight_pair(pair.wilson,pair.background,mode='rge_improved',
                expected_model=model,mu_low=args.mu_low,
                permit_synthetic=args.allow_synthetic_smoke_test,
                allow_projection_residual=args.allow_projection_residual)
            if pf['status']!='ready_conditional':
                raise AdapterBlocked('C preflight: '+' | '.join(pf['blockers'][:4]))
            entry['stage']='single_MS_match_and_SMEFT'
            info=run_rough_case(model,pair,out/slug(model),
                mu_low=args.mu_low,vev=args.vev,overwrite=args.overwrite,
                allow_synthetic=args.allow_synthetic_smoke_test,
                allow_projection_residual=args.allow_projection_residual)
            entry.update(info)
            entry['stage']='complete'
            entry['blocker']=None
            print(f'{model}: {entry["status"]} masses_eV={info["masses_eV"]}',flush=True)
        except Exception as exc:
            entry['status']='BLOCKED'
            entry['blocker']=f'{type(exc).__name__}: {exc}'
            _save_json(out/'_blocked'/f'{slug(model)}.json',entry)
            print(f'{model}: BLOCKED at {entry["stage"]}: {entry["blocker"]}',flush=True)
        # Durable summary after each model: interruption still leaves statuses.
        _write_summary(out,rows,allow_synthetic=args.allow_synthetic_smoke_test)
    counts=json.loads((out/'summary.json').read_text(encoding='utf-8'))['counts']
    print('All-model summary:',out/'summary.json','and',out/'summary.csv',counts,flush=True)
    return 0 if counts['blocked']==0 else 1


if __name__=='__main__':
    sys.exit(main())
