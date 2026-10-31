#!/usr/bin/env python3
"""Collect current full coupled UV reports as CSV; optional lightweight figures.

Only reads the dedicated workstream output. Never treats these quantities as
physical neutrino masses or Matchete-normalized couplings.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))

from Numerical.experimental.T3GeneralYukawaCoupledUV import OUTPUT,require


def summarize(out:Path=OUTPUT,plot:bool=True)->dict:
    source=out/'summary.json'
    require(source.is_file(),f'Missing completion summary: {source}')
    obj=json.loads(source.read_text(encoding='utf-8-sig'))
    require(obj.get('format')=='t3-all-remaining-large-stack-v1','Unrecognized completion report')
    rows=[]
    for item in obj.get('rows',[]):
        if item.get('status')!='CoupledGeneralYukawaRGBetaOnlyUVPassed':
            rows.append({'model':item.get('model'),'status':item.get('status'),'error':item.get('error','')})
            continue
        data=out/'models'/(item['model']+'.json')
        require(data.is_file(),f'Missing per-model trajectory: {data}')
        run=json.loads(data.read_text(encoding='utf-8-sig'))
        require(run.get('status')=='CoupledGeneralYukawaRGBetaOnlyUVPassed','Invalid model result')
        ctl=run.get('seeded_vs_zero_control') or {}
        rows.append({
            'model':run['model'],'status':'Passed',
            'extra_yukawa':run['extra_yukawa'] or '',
            'real_ode_size':run['real_ode_size'],
            'n_coupling_blocks':run['n_coupling_blocks'],'nfev':run['nfev'],
            'mu_initial_gev':run['mu_initial_gev'],'mu_final_gev':run['mu_final_gev'],
            'final_heavy_mass_max_singular_gev':run['final_heavy_mass_max_singular_gev'],
            'y1_relative_change':run['relative_changes']['y1'],
            'y2_relative_change':run['relative_changes']['y2'],
            'MF_relative_change':run['relative_changes']['MF'],
            'extra_final_seeded_norm':ctl.get('extra_final_seeded_frobenius',''),
            'extra_final_zero_control_norm':ctl.get('extra_final_zero_seed_frobenius',''),
            'beta_extra_at_zero_norm':ctl.get('beta_extra_at_zero_frobenius',''),
            'y1_endpoint_seed_control_difference':ctl.get('delta_ordinary_y1_endpoint_frobenius',''),
            'y2_endpoint_seed_control_difference':ctl.get('delta_ordinary_y2_endpoint_frobenius',''),
            'scalar_masses_sq_positive':run['final_scalar_mass_sq_positive'],
        })
    columns=[
        'model','status','extra_yukawa','real_ode_size','n_coupling_blocks',
        'nfev','mu_initial_gev','mu_final_gev','final_heavy_mass_max_singular_gev',
        'y1_relative_change','y2_relative_change','MF_relative_change',
        'extra_final_seeded_norm','extra_final_zero_control_norm','beta_extra_at_zero_norm',
        'y1_endpoint_seed_control_difference','y2_endpoint_seed_control_difference',
        'scalar_masses_sq_positive','error',
    ]
    dest=out/'comparison.csv'
    with dest.open('w',encoding='utf-8',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=columns)
        writer.writeheader()
        for row in rows:writer.writerow(row)
    if plot:
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
        except ImportError:
            print('INFO: matplotlib unavailable; CSV generated without figures')
        else:
            folder=out/'figures';folder.mkdir(exist_ok=True)
            for row in rows:
                if row.get('status')!='Passed':continue
                run=json.loads((out/'models'/(row['model']+'.json')).read_text(encoding='utf-8'))
                seq=run['sampled_trajectory']
                mu=[v['mu_gev'] for v in seq]
                fig,ax=plt.subplots(figsize=(7,4.2))
                names=['y1','y2']+([run['extra_yukawa']] if run['extra_yukawa'] else [])
                for name in names:
                    ax.plot(mu,[v['coupling_frobenius_or_absolute'][name] for v in seq],
                            label=name,marker='.',markersize=3)
                ax.set_xscale('log');ax.set_xlabel('Renormalization scale μ [GeV]')
                ax.set_ylabel('Frobenius norm of Yukawa [RGBeta coordinates]')
                ax.legend();ax.grid(True,alpha=.3)
                fig.tight_layout();fig.savefig(folder/(row['model']+'_uv_yukawas.png'),dpi=150)
                plt.close(fig)
    return {'status':'RGOnlyUVSummaryBuilt','models_in_summary':len(rows),
            'csv':str(dest),'figures_generated':bool(plot)}


def self_test():
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        tmp=Path(tmp);(tmp/'models').mkdir()
        (tmp/'summary.json').write_text(json.dumps({'format':'t3-all-remaining-large-stack-v1',
                                                    'rows':[{'model':'A_p0','status':'Blocked','error':'fixture'}]}))
        result=summarize(tmp,plot=False)
        require(result['models_in_summary']==1,'Fixture summary failed')
        require('fixture' in (tmp/'comparison.csv').read_text(),'Fixture CSV error missing')
    print('PASS: CSV aggregation, blocked-result persistence, no required matplotlib dependency')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--self-test',action='store_true')
    ap.add_argument('--no-plots',action='store_true')
    a=ap.parse_args()
    if a.self_test:self_test()
    else:
        try:
            res=summarize(plot=not a.no_plots)
            print(f'{res["status"]}: {res["models_in_summary"]} records; CSV {res["csv"]}')
        except (ValueError,OSError,KeyError) as exc:
            print(f'BLOCKED: {exc}',file=sys.stderr);raise SystemExit(3)
