"""C-owned CLI for full-Wilson simultaneous scalar matching and SMEFT evolution."""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))

from Numerical.parallel_t3.threshold.models import ALL_MODELS,normalize_model
from Numerical.parallel_t3.threshold.pipeline import run_case
from Numerical.parallel_t3.threshold.fixtures import make_fixture,make_all_fixtures


def parser():
    p=argparse.ArgumentParser(description='Workstream C: ALL heavy scalars matched together at MS')
    p.add_argument('--input',type=Path,help='Directory containing A wilson.json and B background.json')
    p.add_argument('--wilson',type=Path,help='Explicit A canonical full Wilson manifest')
    p.add_argument('--background',type=Path,help='Explicit B full mass/quartic manifest')
    p.add_argument('--all16-input-root',type=Path,help='Directories <slug>/wilson.json and background.json')
    p.add_argument('--make-fixture',help='Create synthetic A/B-like input for this model')
    p.add_argument('--make-fixtures-all16',action='store_true',help='Generate synthetic inputs for all 16')
    p.add_argument('--fixture-dir',type=Path,default=Path('output/parallel_t3/threshold/fixtures'))
    p.add_argument('--output',type=Path,default=Path('output/parallel_t3/threshold/candidate_v1'))
    p.add_argument('--mode',choices=('strict_one_loop','rge_improved'),default='rge_improved')
    p.add_argument('--mu-low',type=float,default=173.)
    p.add_argument('--vev',type=float,default=246.22)
    p.add_argument('--allow-projection-residual',action='store_true',
                   help='Allow SU2-breaking/higher-dimensional projection (always candidate)')
    p.add_argument('--match-only',action='store_true',
                   help='Skip SMEFT ODE for diagnostic use; C5_low is then unchanged proxy')
    p.add_argument('--overwrite',action='store_true')
    return p


def main(argv=None)->int:
    args=parser().parse_args(argv)
    files={}
    if args.make_fixture:
        model=normalize_model(args.make_fixture)
        files[model]=make_fixture(model,args.fixture_dir)
    elif args.make_fixtures_all16:
        files=make_all_fixtures(args.fixture_dir)
    elif args.all16_input_root:
        for model in ALL_MODELS:
            slug=model.replace('_-','_m').replace('_+','_p')
            files[model]=(args.all16_input_root/slug/'wilson.json',
                          args.all16_input_root/slug/'background.json')
    else:
        if args.input:
            w=args.input/'wilson.json';b=args.input/'background.json'
        else:
            if not args.wilson or not args.background:
                parser().error('Provide --input, --wilson plus --background, or --make-fixture')
            w,b=args.wilson,args.background
        model=json.loads(w.read_text(encoding='utf-8'))['model']
        files[normalize_model(model)]=(w,b)
    summaries={}
    for model,(w,b) in files.items():
        target=args.output/model.replace('_-','_m').replace('_+','_p')
        try:
            result=run_case(wilson_manifest=w,background_manifest=b,
                       output_dir=target,mu_low=args.mu_low,
                       vev=args.vev,mode=args.mode,
                       allow_projection_residual=args.allow_projection_residual,
                       match_only=args.match_only,overwrite=args.overwrite)
            summaries[model]={'status':result['status'],'result':str(target/'manifest.json'),
                   'mass_eigenstates':result['scalar_eigenstates_integrated_out'],
                   'mass_singular_values_eV':result['mass_singular_values_eV']}
            print(model,'OK',target/'manifest.json',flush=True)
        except Exception as e:
            summaries[model]={'status':'BLOCKED','error':f'{type(e).__name__}: {e}'}
            print(model,'BLOCKED',type(e).__name__,str(e),file=sys.stderr,flush=True)
    args.output.mkdir(parents=True,exist_ok=True)
    (args.output/'summary.json').write_text(json.dumps({'mode':args.mode,
          'requested_models':list(files),'models':summaries},indent=2)+'\n',encoding='utf-8')
    return 0 if all(x['status']!='BLOCKED' for x in summaries.values()) else 1


if __name__=='__main__':sys.exit(main())
