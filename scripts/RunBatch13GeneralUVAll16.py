"""Evolve all 16 candidate general UV RGBeta systems: 1e7 -> 1e5 GeV.

Intended path: scripts/RunBatch13GeneralUVAll16.py
Run from repository root after RunBatch13UVExportsAll16.py.
This is NOT a physical matching calculation. Illustrative seeds are explicitly
recorded, and Matchete/RGBeta normalisation is NOT certified.
"""
from __future__ import annotations

import argparse
import json
import sys
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Numerical.experimental.Batch13GeneralUVRunning import (
    illustrative_initial_state, run_uv, serialise_state, validate_uv_export,
)

MODELS = ('A_+0','A_-2','A_-4','B_+1','B_-1','B_-3',
          'C_+1','C_-1','C_-3','D_+0','D_+2','D_-2',
          'E_+0','E_+2','E_-2','E_-4')


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exports', type=Path,
                        default=ROOT / 'output/batch13/all16_uv')
    parser.add_argument('--out', type=Path,
                        default=ROOT / 'output/batch13/general_uv_running')
    parser.add_argument('--mu-uv', type=float, default=1e7)
    parser.add_argument('--mu-f', type=float, default=1e5)
    parser.add_argument('--model', choices=MODELS, help='Run one model only')
    parser.add_argument('--rtol', type=float, default=1e-7)
    parser.add_argument('--atol', type=float, default=1e-10)
    args = parser.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    summary = {'status': 'CandidateUVNumerics',
               'physical_neutrino_mass_prediction': False,
               'matchete_rgbeta_normalisation_certified': False,
               'initial_condition_status': 'IllustrativeRGBetaNativeNotPhysicalBenchmark',
               'mu_uv_gev': args.mu_uv, 'mu_f_gev': args.mu_f,
               'models': {}}
    for model in ((args.model,) if args.model else MODELS):
        try:
            path = args.exports / f'{model}_UV.json'
            payload = json.loads(path.read_text(encoding='utf-8-sig'))
            validate_uv_export(payload, model)
            initial = illustrative_initial_state(payload, model, mu_gev=args.mu_uv)
            result = run_uv(payload, model, initial, mu_final_gev=args.mu_f,
                            rtol=args.rtol, atol=args.atol)
            endpoint = {k: v for k,v in result.items() if k != 'state'}
            endpoint['initial'] = serialise_state(initial)
            endpoint['final'] = serialise_state(result['state'])
            (args.out / f'{model}_UV_running.json').write_text(
                json.dumps(endpoint, indent=2, allow_nan=False), encoding='utf-8')
            summary['models'][model] = {'passed': True, 'nfev': result['nfev'],
                                        'n_couplings': result['n_couplings']}
            print(f'[{model}] PASSED | {result["n_couplings"]} couplings, '
                  f'{result["n_parameters_real"]} real ODE variables, '
                  f'nfev={result["nfev"]}', flush=True)
        except Exception as exc:
            summary['models'][model] = {'passed': False, 'error': str(exc)}
            (args.out / f'{model}_UV_running.error.log').write_text(
                traceback.format_exc(), encoding='utf-8')
            print(f'[{model}] FAILED: {exc}', flush=True)
    passed = sum(bool(x['passed']) for x in summary['models'].values())
    summary['passed'] = passed
    summary['failed'] = len(summary['models']) - passed
    (args.out / 'summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    print(f'Summary: {args.out / "summary.json"} | passed={passed} '
          f'failed={summary["failed"]}', flush=True)
    return 0 if summary['failed'] == 0 else 1


if __name__ == '__main__':
    raise SystemExit(main())
