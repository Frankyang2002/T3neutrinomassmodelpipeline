"""Batch 14: export and evolve all 16 candidate scalar-only EFT1 systems.

Project path: scripts/RunBatch14EFT1All16.py

Usage:
  python scripts/RunBatch14EFT1All16.py --stage exports
  python scripts/RunBatch14EFT1All16.py --stage run
  python scripts/RunBatch14EFT1All16.py --stage all

Requires Batch13 UV JSON beta exports and completed UV numerical endpoints.
This is a renormalisable F-threshold projection, NOT Wilson matching.
"""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import traceback

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
from Numerical.experimental.Batch13GeneralUVRunning import REPRESENTATIONS
from Numerical.experimental.Batch14EFT1Running import (
    load_saved_uv_endpoint,run_eft1,serialise_eft1_state,validate_eft1_export,
)
MODELS=('A_+0','A_-2','A_-4','B_+1','B_-1','B_-3',
        'C_+1','C_-1','C_-3','D_+0','D_+2','D_-2',
        'E_+0','E_+2','E_-2','E_-4')


def model_command(model,output,wolframscript='wolframscript'):
    family,alpha_string=model.split('_')
    d1,d2,df=REPRESENTATIONS[family]
    alpha=int(alpha_string)
    return [wolframscript,'-file',str(ROOT/'RGE/running/rgbeta/RunGeneralT3QuadraticRGBeta.wl'),
            'EFT1',str(d1),str(d2),str(df),f'm{abs(alpha)}' if alpha<0 else str(alpha),
            model,str(output)]


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage',choices=('exports','run','all'),default='all')
    parser.add_argument('--model',choices=MODELS)
    parser.add_argument('--uv-exports',type=Path,default=ROOT/'output/batch13/all16_uv')
    parser.add_argument('--uv-endpoints',type=Path,default=ROOT/'output/batch13/general_uv_running')
    parser.add_argument('--out',type=Path,default=ROOT/'output/batch14/eft1')
    parser.add_argument('--mu-s',type=float,default=1e3)
    parser.add_argument('--wolframscript',default='wolframscript')
    parser.add_argument('--rtol',type=float,default=1e-7)
    parser.add_argument('--atol',type=float,default=1e-10)
    args=parser.parse_args(argv)
    args.out.mkdir(parents=True,exist_ok=True)
    models=(args.model,) if args.model else MODELS
    summary={'status':'CandidateEFT1RenormalisableOnly',
             'physical_neutrino_mass_prediction':False,
             'fermion_threshold_wilson_matching_performed':False,
             'both_scalars_active_in_eft1':True,
             'models':{}}
    for model in models:
        try:
            uv_payload=json.loads((args.uv_exports/f'{model}_UV.json').read_text(encoding='utf-8-sig'))
            eft_path=args.out/f'{model}_EFT1.json'
            if args.stage in ('exports','all'):
                # Never allow an old JSON to mask a fresh Wolfram failure.
                eft_path.unlink(missing_ok=True)
                env=os.environ.copy()
                env['T3_EXPERIMENTAL_GENERAL_FERMIONS']='0'
                cmd=model_command(model,eft_path,args.wolframscript)
                proc=subprocess.run(cmd,cwd=ROOT,text=True,stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT,env=env,check=False)
                (args.out/f'{model}_EFT1.log').write_text(proc.stdout,encoding='utf-8')
                if proc.returncode!=0 or not eft_path.is_file():
                    raise RuntimeError(f'EFT1 Wolfram export exit={proc.returncode}; see {model}_EFT1.log')
            eft_payload=json.loads(eft_path.read_text(encoding='utf-8-sig'))
            validate_eft1_export(eft_payload,uv_payload,model)
            if args.stage in ('run','all'):
                uv_state=load_saved_uv_endpoint(args.uv_endpoints/f'{model}_UV_running.json',uv_payload,model)
                result=run_eft1(uv_state,uv_payload,eft_payload,model,
                                mu_scalar_gev=args.mu_s,rtol=args.rtol,atol=args.atol)
                output={k:v for k,v in result.items() if k not in ('boundary','state')}
                output['initial']=serialise_eft1_state(result['boundary'])
                output['final']=serialise_eft1_state(result['state'])
                (args.out/f'{model}_EFT1_running.json').write_text(
                    json.dumps(output,indent=2,allow_nan=False),encoding='utf-8')
                print(f'[{model}] PASSED | EFT1 {result["n_couplings"]} couplings, '
                      f'{result["n_parameters_real"]} real ODE variables, nfev={result["nfev"]}',flush=True)
            else:
                print(f'[{model}] PASSED | EFT1 export',flush=True)
            summary['models'][model]={'passed':True}
        except Exception as exc:
            summary['models'][model]={'passed':False,'error':str(exc)}
            (args.out/f'{model}_error.log').write_text(traceback.format_exc(),encoding='utf-8')
            print(f'[{model}] FAILED: {exc}',flush=True)
    summary['passed']=sum(v['passed'] for v in summary['models'].values())
    summary['failed']=len(models)-summary['passed']
    (args.out/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(f'Summary: {args.out/"summary.json"} | passed={summary["passed"]} failed={summary["failed"]}')
    return 0 if summary['failed']==0 else 1


if __name__=='__main__':raise SystemExit(main())
