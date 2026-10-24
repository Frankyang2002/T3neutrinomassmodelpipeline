"""Compare existing T3-B-m1 fixed-order result with opt-in coupled C12/C5.

Writes a separate diagnostic JSON. Does not rerun optimisation or replace
production outputs.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from Numerical.orchestration import RunSmallLambdaRefit as base
from Numerical.orchestration.T3BIntegratedCoupledPipeline import compute_opt_in_coupled

MODEL = 'T3_dS1_2_dS2_2_dF_1_alpha_m1'


def load(path: Path) -> dict:
    if not path.is_file():
        raise FileNotFoundError(f'Required file absent: {path}')
    return json.loads(path.read_text(encoding='utf-8-sig'))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--mode', choices=['optimised', 'shared_a_m2'], default='optimised')
    ap.add_argument('--config', type=Path, help='Override evaluated configuration')
    ap.add_argument('--diagnostics', type=Path, help='Matching production running_diagnostics.json')
    ap.add_argument('--artifact-map', type=Path)
    ap.add_argument('--samples', type=int, default=120)
    ap.add_argument('--max-order', type=int, default=4)
    args = ap.parse_args()

    folder = (base.ROOT / 'output/studies/lambda_5e-9_refit/comparison'
              / args.mode / 'models' / MODEL)
    cfg_path = args.config or folder / 'shared_config.json'
    diag_path = args.diagnostics or folder / 'running_diagnostics.json'
    cfg, production = load(cfg_path), load(diag_path)
    rep = cfg['representation']
    if (rep['d_s1'],rep['d_s2'],rep['d_f'],rep['alpha'],bool(rep.get('shared_scalar'))) != (2,2,1,-1,False):
        raise ValueError('Configuration must be ordinary T3-B-m1')
    mapping = load(args.artifact_map) if args.artifact_map else None
    uv, eft, c5 = base.locate_artifacts(MODEL, mapping)
    result = compute_opt_in_coupled(
        config=cfg, uv_rgbeta=uv, eft1_rgbeta=eft, final_c5=c5,
        production=production, max_order=args.max_order, n_samples=args.samples)

    dest = (base.ROOT / 'output/studies/t3b_coupled_fixed_order_comparison'
            / args.mode / 'coupled_vs_fixed_order.json')
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        raise FileExistsError(f'Will not overwrite: {dest}')
    result['comparison_inputs'] = {
        'config':str(cfg_path), 'production':str(diag_path),
        'uv_rgbeta':str(uv), 'eft1_rgbeta':str(eft), 'final_c5':str(c5)}
    dest.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('Saved:',dest)
    print('Scalar threshold relative C5 difference:',result['coupled_at_ms']['relative_difference_from_fixed_order'])
    print('Low energy relative C5 difference:',result['low_energy']['relative_c5_difference_at_low_scale'])
    print('Reverse feedback relative:',result['coupled_at_ms']['relative_reverse_feedback'])
    print('Low-energy observable differences:',json.dumps(result['low_energy']['differences'],indent=2))
    print('Validated physical prediction:',result['validated_physical_prediction_available'])


if __name__ == '__main__':
    main()
