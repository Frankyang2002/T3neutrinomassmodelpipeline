"""Narrow T1/T2/T3 gamma sensitivity for the saved zero-extra control.

Reuses the *same* perturbed-Yukawa C5 calculation as the frozen-64 campaign;
never mixes held-zero, boundary-zero or old gamma scan data.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

from scripts import T3Frozen64Sources as policy
from scripts import T3ZeroExtraReportSources as zero


def run(run_id='zero_extra_all16_low_v1', scan_id='gamma_095_105_v1',
        resume=False, mode='boundary_zero'):
    policy.require(bool(policy.RUN_ID.fullmatch(scan_id)), 'invalid scan ID')
    policy.require(mode in (*zero.MODES, 'all'), 'invalid --mode')
    selected_modes = zero.MODES if mode == 'all' else (mode,)
    from scripts import GenerateT3FrozenFullScalarGamma as old
    meta=zero.load_summary(run_id)
    root=zero.run_root(run_id)
    grid=zero.NARROW_GAMMA
    generator_hash=policy.digest(Path(__file__))
    compute_hash=policy.digest(Path(old.__file__))
    complete=0
    for mode in selected_modes:
        for model in policy.MODELS:
            src=zero.sources(run_id,meta,mode,model)
            cell=src['cell']
            csvpath=root/scan_id/mode/(policy.tag(model)+'.csv')
            metapath=csvpath.with_suffix('.json')
            identity={'model':model,'mode':mode,'benchmark':meta['benchmark'],
                      'gamma':list(grid),'trajectory_sha256':cell['trajectory_sha256'],
                      'neutrino_sha256':cell['neutrino_sha256'],
                      'source_summary_sha256':meta['summary_sha256'],
                      'source_settings_sha256':meta['settings_sha256'],
                      'generator_sha256':generator_hash,
                      'frozen_gamma_calculator_sha256':compute_hash}
            if csvpath.exists() or metapath.exists():
                policy.require(resume and csvpath.is_file() and metapath.is_file(),
                               'partial/old gamma scan: '+str(csvpath))
                saved=policy.read_json(metapath)
                policy.require(saved.get('source')==identity and
                               saved.get('csv_sha256')==policy.digest(csvpath) and
                               saved.get('failures')==0 and saved.get('rows')==33,
                               'gamma scan provenance mismatch; choose a new --scan-id')
                complete+=1
                continue
            # BC matching reads <source_root>/<model_tag>/*_eft1_quadratic_running.json.
            # Frozen-64's old.compute expects the model directory as its argument.
            model_folder = src['source_root'] / policy.tag(model)
            rows,failed=old.compute(cell,grid,model_folder,csvpath)
            saved={'source':identity,'csv_sha256':policy.digest(csvpath),
                   'rows':rows,'failures':failed}
            metapath.write_text(json.dumps(saved,indent=2)+'\n',encoding='utf-8')
            policy.require(rows==33 and failed==0,
                           f'{mode}/{model}: failed gamma points {failed}/33; inspect {csvpath}')
            complete+=1
            print(f'{mode}/{model}: 33 narrow-gamma evaluations',flush=True)
    result={'run_id':run_id,'scan_id':scan_id,'mode_selection':list(selected_modes),
            'complete':complete,'expected':len(selected_modes)*len(policy.MODELS),
            'gamma_grid':list(grid),'points_per_mode':len(policy.MODELS)*33}
    (root/scan_id/'summary.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    return result


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-id',default='zero_extra_all16_low_v1')
    p.add_argument('--scan-id',default='gamma_095_105_v1')
    p.add_argument('--resume',action='store_true')
    p.add_argument('--mode',choices=(*zero.MODES,'all'),default='boundary_zero',
                   help='Which zero-extra control to scan (default: boundary_zero)')
    a=p.parse_args(argv)
    print(json.dumps(run(a.run_id,a.scan_id,a.resume,a.mode),indent=2))
    return 0


if __name__=='__main__':
    raise SystemExit(main())
