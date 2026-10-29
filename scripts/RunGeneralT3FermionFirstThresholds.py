"""Stage 16: real two-segment renormalisable running and legacy baseline C5."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from Numerical.experimental.GeneralT3HierarchicalTrajectory import run_hierarchical_stage16


def load(path: Path):
    if not path.is_file():
        raise FileNotFoundError(f'Missing required file: {path}')
    value = json.loads(path.read_text(encoding='utf-8-sig'))
    if not isinstance(value, dict):
        raise ValueError(f'Expected JSON object: {path}')
    return value


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config', type=Path, default=Path('docs/general_t3_threshold_benchmark.json'))
    ap.add_argument('--rgbeta-root', type=Path, default=Path('output'))
    ap.add_argument('--output', type=Path, default=Path('output/general_t3_stage16'))
    ap.add_argument('--renormalisable-only', action='store_true',
                    help='Skip pre-existing baseline Weinberg matching and low-energy evolution')
    ns = ap.parse_args()
    ns.output.mkdir(parents=True, exist_ok=True)
    failures = []
    summaries = []
    for model, folder in [('B_m1', 'T3_B_alpha_m1'), ('C_m1', 'T3_C_alpha_m1')]:
        base = ns.rgbeta_root / folder / 'data'
        try:
            result = run_hierarchical_stage16(
                model, load(ns.config),
                load(base/'uv_rgbeta_rge.json'),
                load(base/'eft1_rgbeta_rge.json'),
                final_weinberg=None if ns.renormalisable_only else base/'final_weinberg_coefficient.json',
            )
            path = ns.output/f'{model}_hierarchical_trajectory.json'
            path.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
            print(model, result['status'], 'UV nfev=', result['uv_nfev'],
                  'EFT1 nfev=', result['eft1_nfev'])
            print('  extra C5:', result['extra_c5_status'])
            if result['legacy_baseline_neutrino_mass_at_mu_low_ev'] is not None:
                print('  legacy baseline ||mnu(mu_low)||_F =',
                      result['legacy_baseline_neutrino_mass_at_mu_low_ev']['frobenius_norm'], 'eV')
            print('Report:', path)
            summaries.append({'model': model, 'status': result['status'],
                              'uv_nfev': result['uv_nfev'],
                              'eft1_nfev': result['eft1_nfev'],
                              'baseline_mnu_norm_ev': None if result['legacy_baseline_neutrino_mass_at_mu_low_ev'] is None
                                      else result['legacy_baseline_neutrino_mass_at_mu_low_ev']['frobenius_norm'],
                              'extra_c5_status': result['extra_c5_status']})
        except (ValueError, TypeError, KeyError, OSError, RuntimeError, NotImplementedError) as err:
            failures.append({'model': model, 'error': str(err)})
            print('BLOCKED', model, err, file=sys.stderr)
    summary = {'status': 'RenormalisablePathsPassed' if not failures else 'ReviewRequired',
               'models': summaries, 'failures': failures,
               'extended_wilson_matching_implemented': False,
               'production_files_modified': False}
    (ns.output/'stage16_summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    return 0 if not failures else 2

if __name__ == '__main__':
    raise SystemExit(main())
