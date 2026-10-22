"""Compare staged RGE contributions for the exact 5e-10 A(-2) prepared seed.

Run from repository root:
    python -m Numerical.orchestration.DiagnoseLambda5e10UV

Uses saved UV/EFT1 RGBeta and final C5 artifacts; does not fit parameters.
Matches thresholds identically in every variant. Comparisons use one common
fully evolved low-energy charged-lepton basis (as in RunFittedStagedRGComparison).
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from Numerical.orchestration import RunSmallLambdaRefit as base
from Numerical.orchestration import RunLambda5e10Refit as study
from Numerical.orchestration.studies.fitted_rge.RunFittedStagedRGComparison import (
    evaluate_staged,
    summarise,
)
from Numerical.core.RGBetaEvaluator import load_rgbeta_payload
from Numerical.running.FinalC5TrajectoryAdapter import FinalC5TrajectoryEvaluator
from Numerical.fitting.OscillationFit import OscillationFitTarget

MODEL = 'T3_dS1_1_dS2_3_dF_2_alpha_m2'
OBS = (
    'delta_m21_sq_ev2', 'delta_m3l_sq_ev2',
    'sin2_theta12', 'sin2_theta13', 'sin2_theta23',
)
STAGES = ('000', '100', '110', '111')


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--artifact-map', type=Path, default=None,
                        help='Optional artifact-map JSON in the original study format')
    args = parser.parse_args()

    root = base.ROOT / 'output' / 'studies' / 'lambda_5e-10_refit'
    prepared_path = root / 'prepared' / (MODEL + '.json')
    output = root / 'diagnostics' / 'a_m2_exact_seed_uv_ablation'
    if not prepared_path.is_file():
        parser.error(f'Missing prepared seed: {prepared_path}')
    cfg = base.load(prepared_path)
    if cfg.get('representation', {}).get('alpha') != -2:
        raise ValueError('Unexpected A(-2) representation in seed')
    lam = float(cfg['base_state']['ordinary']['lambdaT3']['real'])
    if not 4.999e-10 <= lam <= 5.001e-10:
        raise ValueError(f'Expected lambdaT3=5e-10, found {lam}')

    mapping = base.load(args.artifact_map) if args.artifact_map else None
    artifacts = base.locate_artifacts(MODEL, mapping)
    uv_beta = load_rgbeta_payload(artifacts[0])
    eft_beta = load_rgbeta_payload(artifacts[1])
    c5_evaluator = FinalC5TrajectoryEvaluator(artifacts[2])
    target = OscillationFitTarget.from_json(base.TARGET)

    print('Evaluating four RGE configurations for exact A(-2) seed ...', flush=True)
    stages = evaluate_staged(cfg, uv_beta, eft_beta, c5_evaluator, target)
    record = {
        'model': MODEL,
        'lambdaT3': lam,
        'prepared_config': str(prepared_path),
        'artifacts': [str(p) for p in artifacts],
        'description': (
            '000 none; 100 UV only; 110 UV+EFT1; 111 all. '
            'All stages retain threshold matching, and all are expressed '
            'in the same low-energy charged-lepton basis taken from stage 111.'
        ),
        'stages': stages,
    }
    summarise(record)
    output.mkdir(parents=True, exist_ok=True)
    json_path = output / 'ablation.json'
    csv_path = output / 'ablation.csv'
    json_path.write_text(json.dumps(record, indent=2, allow_nan=False) + '\n', encoding='utf-8')

    with csv_path.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.writer(stream)
        writer.writerow(['observable', '000_frozen', '100_uv_only', '110_uv_eft1',
                         '111_full', 'uv_percent', 'eft1_percent',
                         'weinberg_percent', 'full_percent'])
        for name in OBS:
            writer.writerow([
                name, *(stages[s]['prediction'][name] for s in STAGES),
                *(record['increment_percent'][x][name]
                  for x in ('UV', 'EFT1', 'Weinberg')),
                record['total_percent'][name],
            ])
    print('\nConfiguration chi-square:')
    for s in STAGES:
        print(f"  {s}: chi2={stages[s]['chi2']:.8g}")
    print('\nObservable comparison (000 -> 100 isolates the UV running increment):')
    for name in OBS:
        v0 = stages['000']['prediction'][name]
        v1 = stages['100']['prediction'][name]
        vf = stages['111']['prediction'][name]
        pct = record['increment_percent']['UV'][name]
        print(f'  {name:23s} frozen={v0:.9g} uv_only={v1:.9g} full={vf:.9g}'
              f' uv_delta={pct:+.4g}%' if pct is not None else
              f'  {name:23s} frozen={v0:.9g} uv_only={v1:.9g} full={vf:.9g}')
    print(f'\nWrote: {json_path}\n       {csv_path}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
