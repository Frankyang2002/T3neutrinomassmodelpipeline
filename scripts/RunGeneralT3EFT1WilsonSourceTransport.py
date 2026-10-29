"""Stage18: numerical full-flavour additional EFT1 LLHH Wilson *source* insertion."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from Numerical.experimental.GeneralT3ExtraWilsonTransport import run_stage18_model


def _read(path:Path):
    if not path.is_file():
        raise FileNotFoundError(f'Missing input: {path}')
    value=json.loads(path.read_text(encoding='utf-8-sig'))
    if not isinstance(value,dict):
        raise ValueError(f'Expected JSON object in {path}')
    return value


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',type=Path,default=Path('docs/general_t3_threshold_benchmark.json'))
    p.add_argument('--rgbeta-root',type=Path,default=Path('output'))
    p.add_argument('--stage17',type=Path,default=Path('output/general_t3_stage17'))
    p.add_argument('--output',type=Path,default=Path('output/general_t3_stage18'))
    p.add_argument('--samples',type=int,default=17)
    args=p.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    failures=[]; models=[]
    for model, folder in (('B_m1','T3_B_alpha_m1'),('C_m1','T3_C_alpha_m1')):
        root=args.rgbeta_root/folder/'data'
        try:
            result=run_stage18_model(model,_read(args.config),
                _read(root/'uv_rgbeta_rge.json'),
                _read(root/'eft1_rgbeta_rge.json'),
                _read(args.stage17/f'{model}_extra_eft1_scalar_pair_beta.json'),
                n_save=args.samples)
            report=args.output/f'{model}_extra_eft1_wilson_source_integral.json'
            report.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
            norms=[x['integrated_delta_C_ijab_gev_inverse']['frobenius_norm_gev_inverse']
                   for x in result['higgs_real_components']]
            print(model,result['status'], 'UV nfev=',result['uv_nfev'],
                  'EFT1 nfev=',result['eft1_nfev'])
            print('  LLHH component Frobenius norms (GeV^-1):', norms)
            print('  Physical C5:', 'BLOCKED: projection and scalar matching unverified')
            print('Report:',report)
            models.append({'model':model,'status':result['status'],
                          'component_norms_gev_inverse':norms,
                          'n_samples':result['eft1_samples']})
        except (ValueError,TypeError,KeyError,OSError,RuntimeError,NotImplementedError) as err:
            failures.append({'model':model,'error_type':type(err).__name__,'error':str(err)})
            print('BLOCKED',model,f'{type(err).__name__}: {err}',file=sys.stderr)
    summary={'status':'ExtraEFT1NumericalWilsonSourcesIntegrated' if not failures else 'ReviewRequired',
             'models':models,'failures':failures,
             'physical_extra_c5_implemented':False,
             'production_files_modified':False}
    (args.output/'stage18_summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return 0 if not failures else 2

if __name__=='__main__':
    raise SystemExit(main())
