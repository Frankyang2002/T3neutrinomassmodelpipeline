"""Stage19: certify and transport selected extra EFT1 Weinberg source for B/C(-1)."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from Numerical.experimental.GeneralT3ExtraLowEnergyTransport import transport_stage19_model


def read(path:Path) -> dict:
    if not path.is_file():raise FileNotFoundError(f'Required data missing: {path}')
    data=json.loads(path.read_text(encoding='utf-8-sig'))
    if not isinstance(data,dict):raise ValueError(f'Invalid JSON object: {path}')
    return data


def main()->int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',type=Path,default=Path('docs/general_t3_threshold_benchmark.json'))
    p.add_argument('--output-root',type=Path,default=Path('output'))
    p.add_argument('--stage17',type=Path,default=Path('output/general_t3_stage17'))
    p.add_argument('--stage18',type=Path,default=Path('output/general_t3_stage18'))
    p.add_argument('--out',type=Path,default=Path('output/general_t3_stage19'))
    a=p.parse_args()
    a.out.mkdir(parents=True,exist_ok=True)
    failures=[]; models=[]
    for model,base_name in (('B_m1','T3_B_alpha_m1'),('C_m1','T3_C_alpha_m1')):
        try:
            base=a.output_root/base_name/'data'
            result=transport_stage19_model(model,read(a.config),
                read(base/'uv_rgbeta_rge.json'),read(base/'eft1_rgbeta_rge.json'),
                read(a.stage17/f'{model}_extra_eft1_scalar_pair_beta.json'),
                read(a.stage18/f'{model}_extra_eft1_wilson_source_integral.json'),
                read(base/'eft1_after_F_wilson_seed.json'),
                read(base/'eft1_wilson_rge.json'),
                read(base/'final_weinberg_coefficient.json'))
            path=a.out/f'{model}_selected_extra_weinberg_transport.json'
            path.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
            x=result['selected_extra_mass_at_mu_low_ev']['frobenius_norm']
            print(f'{model}: {result["status"]}, projected source-only ||delta mnu(mu_low)||_F={x:.9e} eV')
            print('  SU(2) components checked:',result['projection_certificate']['independent_components_verified'])
            print('  incomplete two-loop calculation:',not result['complete_two_loop_result'])
            print('Report:',path)
            models.append({'model':model,'status':result['status'],
                           'source_only_low_mass_norm_ev':x,
                           'template_residual':result['integral_projection']['projection_relative_residual']})
        except Exception as err:
            print(f'BLOCKED {model}: {type(err).__name__}: {err}',file=sys.stderr)
            failures.append({'model':model,'type':type(err).__name__,'error':str(err)})
    summary={'status':'SourceProjectionAndLowEnergyTransportVerified' if not failures else 'ReviewRequired',
             'models':models,'failures':failures,
             'complete_two_loop_result':False,'production_files_modified':False}
    (a.out/'stage19_summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return 0 if not failures else 2

if __name__=='__main__':raise SystemExit(main())
