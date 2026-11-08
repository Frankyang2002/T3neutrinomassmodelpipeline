"""A(-4) matrix-free seven-term Wilson ODE; requires Chat A/B contract JSON.

Results are *never* described as matched SMEFT C5. Refuses overwritten outputs.
Use --allow-synthetic for the included deliberately unphysical fixture.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from RGE.running.intermediate.GeneralEFT1SevenTermFlow import (
    run_contracts, save_results, evolve, c5_projection,
)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--boundary', type=Path, required=True,
                   help='Chat A certified-master-Weyl boundary (or synthetic fixture)')
    p.add_argument('--couplings', type=Path, required=True,
                   help='Chat B full scale-dependent master coupling samples')
    p.add_argument('--output-root', type=Path, required=True)
    p.add_argument('--allow-synthetic', action='store_true')
    p.add_argument('--beta-mode', choices=('source_v1_unstarred','ward_covariant_candidate','fonseca_eq22_six_term'),
                   default='source_v1_unstarred',
                   help='Source, crossed C-star, or independent six-term alternative; NOT physical-certified')
    p.add_argument('--points', type=int, default=11)
    p.add_argument('--rtol', type=float, default=1e-8)
    p.add_argument('--atol', type=float, default=1e-18)
    p.add_argument('--source-decomposition', action='store_true',
                   help='Evolve each canonical operator source independently; cost scales with sources')
    p.add_argument('--convergence-check', action='store_true',
                   help='Repeat integration with tighter tolerance and compare endpoints')
    args = p.parse_args(argv)
    target = args.output_root / 'A_m4'
    summary = target / 'seven_term_summary.json'
    arrays = target / 'seven_term_trajectory.npz'
    if summary.exists() or arrays.exists():
        raise FileExistsError('Refusing to overwrite existing Batch34 outputs')
    bbytes = args.boundary.read_bytes()
    cbytes = args.couplings.read_bytes()
    b = json.loads(bbytes.decode('utf-8-sig'))
    coupling = json.loads(cbytes.decode('utf-8-sig'))
    kw = dict(allow_synthetic=args.allow_synthetic, points=args.points, rtol=args.rtol,
              atol=args.atol, beta_mode=args.beta_mode)
    result = run_contracts(b, coupling, source_decomposition=args.source_decomposition, **kw)
    if args.convergence_check:
        extra = run_contracts(b,coupling,**{**kw,'rtol':kw['rtol']/10,'atol':kw['atol']/10})
        import numpy as np
        result['solver']['convergence_relative_endpoint'] = float(np.linalg.norm(
            result['final']-extra['final'])/max(np.linalg.norm(extra['final']),1e-300))
        result['solver']['convergence_refined_nfev'] = int(extra['solver']['nfev'])
    save_results(result,str(arrays),str(summary))
    report=json.loads(summary.read_text())
    report['source_boundary_sha256']=sha256(bbytes).hexdigest()
    report['source_couplings_sha256']=sha256(cbytes).hexdigest()
    report['source_boundary_path']=str(args.boundary)
    report['source_couplings_path']=str(args.couplings)
    summary.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('RGE beta mode:', report['beta_mode'], '; physically certified: NO')
    print(f'Wrote {arrays} and {summary}; generated zero-boundary coordinates: '
          f'{report["generated_initially_zero_coordinate_count"]}; matched SMEFT C5: NO')
    return 0

if __name__=='__main__':
    raise SystemExit(main())
