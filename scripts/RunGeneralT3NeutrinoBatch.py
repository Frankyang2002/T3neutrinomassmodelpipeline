"""Use per-model matched C5 and full-general running snapshots for numerical spectra.

Defaults to the 14 successfully matchable models, excluding B_-1 and C_-1.
If a model has no running JSON, the existing GeneralRunningBenchmark can be run
with --run-running (requires canonical symbolic RGE exports for each model).
"""
from __future__ import annotations
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import subprocess
import sys
from Numerical.matching.GeneralT3MultiModelNeutrinos import parser as numerical_parser, run

MODELS=['A_-4','A_-2','A_+0','B_-3','B_+1','C_-3','C_+1',
        'D_-2','D_+0','D_+2','E_-4','E_-2','E_+0','E_+2']
BLOCKED=['B_-1','C_-1']

def tag(m):return m.replace('+','p').replace('-','m')

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--model', choices=MODELS)
    p.add_argument('--matching-root',type=Path,default=Path('output/general_physical_matching_batch'))
    p.add_argument('--running-root',type=Path,default=Path('output/general_physical_matching_batch'))
    p.add_argument('--output-root',type=Path,default=Path('output/general_numerical_all14'))
    p.add_argument('--run-running',action='store_true')
    p.add_argument('--canonical-run-id',default='full16_general_v1')
    p.add_argument('--generate-symbolic',action='store_true')
    p.add_argument('--mu-uv',type=float,default=1e7)
    p.add_argument('--mu-f',type=float,default=1e5)
    p.add_argument('--mu-s',type=float,default=1e3)
    p.add_argument('--mu-low',type=float,default=91.1876)
    p.add_argument('--initial-json',type=Path)
    p.add_argument('--kappa-real',type=float)
    p.add_argument('--kappa-imag',type=float,default=0.)
    p.add_argument('--quartic-real',type=float,default=0.02)
    p.add_argument('--quartic-imag',type=float,default=0.)
    p.add_argument('--strict-mapping',action='store_true')
    p.add_argument('--wolframscript',default='wolframscript')
    args=p.parse_args()
    rows=[]
    models=[args.model] if args.model else MODELS
    args.output_root.mkdir(parents=True,exist_ok=True)
    for model in models:
        t=tag(model); folder=args.output_root/t
        folder.mkdir(parents=True,exist_ok=True)
        coeff=args.matching_root/t/'c5_coefficient.txt'
        running=args.running_root/t/'general_running.json'
        if model=='A_-2' and not running.is_file():
            other=Path('output/general_A_m2_running_better.json')
            if other.is_file():running=other
        row={'model':model,'coefficient':str(coeff),'running':str(running)}
        try:
            if not coeff.is_file():raise FileNotFoundError(f'No matching coefficient: {coeff}')
            if args.run_running:
                running=args.running_root/t/'general_running.json'
                running.parent.mkdir(parents=True,exist_ok=True)
                cmd=[sys.executable,'-X','utf8','-m','Numerical.orchestration.GeneralRunningBenchmark',
                     '--canonical-run-id',args.canonical_run_id,'--model',model,
                     '--output-json',str(running),'--mu-uv',str(args.mu_uv),
                     '--mu-f',str(args.mu_f),'--mu-s',str(args.mu_s)]
                if args.generate_symbolic:cmd.append('--generate-symbolic')
                if args.initial_json:cmd += ['--initial-json',str(args.initial_json)]
                cp=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
                (folder/'running.log').write_text(cp.stdout,encoding='utf-8')
                if cp.returncode:raise RuntimeError('Running failed: '+cp.stdout[-2000:])
            if not running.is_file():
                raise FileNotFoundError(f'No running output: {running}. Use --run-running to generate.')
            num=numerical_parser().parse_args(['--running-json',str(running),
                    '--coefficient',str(coeff),'--output-json',str(folder/'neutrinos.json'),
                    '--mu-low',str(args.mu_low),'--wolframscript',args.wolframscript,
                    '--quartic-real',str(args.quartic_real),'--quartic-imag',str(args.quartic_imag)])
            num.use_running_bilinear=True
            num.strict_mapping=args.strict_mapping
            num.kappa_real=args.kappa_real
            num.kappa_imag=args.kappa_imag
            result=run(num)
            row.update(status='Success',masses_ev=result['neutrino_masses_ev'],
                       low_energy_rge_status=result['low_energy_rge_status'],
                       output_json=str(folder/'neutrinos.json'))
        except Exception as exc:
            row.update(status='Failed',error=str(exc))
        rows.append(row)
        manifest={'generated_utc':datetime.now(timezone.utc).isoformat(),
                  'numerical_success_count':sum(x['status']=='Success' for x in rows),
                  'attempted':len(rows),'known_blocked':BLOCKED,
                  'physical_normalization_certified':False,'models':rows}
        (args.output_root/'summary.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
        print(model,row['status'],row.get('error',''),flush=True)
    return 0 if all(row['status']=='Success' for row in rows) else 1

if __name__=='__main__':raise SystemExit(main())
