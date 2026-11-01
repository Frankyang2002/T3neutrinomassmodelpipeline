"""Experimental B(-1)/C(-1) general running and scalar mass-basis handoff.

This uses the same canonical UV/EFT1 numerical runner and improved benchmark
inputs as the 14-model study. It diagonalizes the mixed scalar quadratic
form at the EFT1 endpoint, but deliberately does NOT claim C5 matching.
No archived C5 is imported. No GitHub files are modified.

Run: python -X utf8 -m scripts.RunGeneralT3BCNormal --all
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np

from Numerical.orchestration.GeneralRunningBenchmark import (
    load_canonical_general_exports, run,
)
from Numerical.experimental.GeneralT3MixedScalarMassBasis import (
    MixedDoubletMassBasis, quadratic_mass_tensor_transport,
    verify_su2_covariance,
)
from scripts.RunGeneralT3ImprovedRunning import improved_overrides

ROOT = Path(__file__).resolve().parents[1]
MODELS = ('B_-1', 'C_-1')
OUT = ROOT / 'output/general_numerical_BC_normal'


def complex_value(raw, key):
    if isinstance(raw, dict) and set(raw) == {'real', 'imag'}:
        z = complex(float(raw['real']), float(raw['imag']))
    elif isinstance(raw, (float, int)):
        z = complex(raw)
    else:
        raise ValueError(f'{key}: invalid scalar {raw!r}')
    if not np.isfinite([z.real, z.imag]).all():
        raise ValueError(f'{key}: nonfinite')
    return z


def save(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def process(model: str, args) -> dict:
    tag = model.replace('-', 'm')
    uv, eft, provenance = load_canonical_general_exports(model, args.canonical_run_id)
    if provenance.get('mass_basis') is None:
        raise ValueError('Canonical export has no mass-basis metadata')
    running = run(uv, eft, mu_uv=args.mu_uv, mu_f=args.mu_f, mu_s=args.mu_s,
                  overrides=improved_overrides(model))
    if running['model'] != model or running['status'] != 'NumericalRenormalisableRunningOnly':
        raise RuntimeError('Unexpected running status')
    if not running['uv_solver']['success'] or not running['eft1_solver']['success']:
        raise RuntimeError('One or both RGE solvers failed')
    state = running['eft1_at_scalar_threshold']
    coupling = 'muT3GeneralS1S2'
    if coupling not in state:
        raise ValueError(f'EFT1 export lacks the required mixed bilinear {coupling}')
    mu2 = complex_value(state[coupling], coupling)
    m1 = float(state['mS1Sq'])
    m2 = float(state['mS2Sq'])
    basis = MixedDoubletMassBasis.build(model, m1, m2, mu2)
    tensor = quadratic_mass_tensor_transport(basis)
    covariance = verify_su2_covariance(basis)
    if basis.tachyonic:
        raise ValueError('Scalar mass eigenstate is tachyonic; threshold matching invalid')
    if basis.degenerate:
        raise ValueError('Scalar mass eigenstates degenerate; basis not unique')
    running['symbolic_provenance'] = provenance
    running['benchmark'] = 'improved_exploratory_BC_mixed_scalar'
    running['uv_benchmark_overrides'] = improved_overrides(model)
    running['physical_matching_performed'] = False
    running_path = args.output / 'running' / f'{tag}.json'
    save(running_path, running)
    output = {
        'model': model,
        'status': 'SuccessRunningAndScalarMassBasis_NoC5Matching',
        'benchmark': running['benchmark'],
        'scale_gev': args.mu_s,
        'basis': basis.info(),
        'quadratic_tensor_transport': tensor,
        'su2_covariance_max_absolute_error': covariance,
        'general_uv_eft1_running_json': str(running_path),
        'matching_status': 'BLOCKED_MixedScalarSequentialC5NotImplemented',
        'C5_computed': False,
        'neutrino_masses_computed': False,
        'physical_predictions_certified': False,
        'archived_C5_used': False,
    }
    save(args.output / f'{tag}_mass_basis_handoff.json', output)
    print(f'{model}: UV/EFT1 PASS; scalar masses [GeV]={basis.info()["masses_gev"]}; C5 NOT MATCHED', flush=True)
    return {'model': model, 'status': output['status'],
            'running_json': str(running_path),
            'mass_basis_json': str(args.output / f'{tag}_mass_basis_handoff.json')}


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    select = p.add_mutually_exclusive_group(required=True)
    select.add_argument('--all', action='store_true')
    select.add_argument('--model', choices=MODELS)
    p.add_argument('--canonical-run-id', default='full16_general_v1')
    p.add_argument('--output', type=Path, default=OUT)
    p.add_argument('--mu-uv', type=float, default=1e7)
    p.add_argument('--mu-f', type=float, default=1e5)
    p.add_argument('--mu-s', type=float, default=1e3)
    args = p.parse_args(argv)
    args.output = args.output.resolve()
    rows = []
    for model in MODELS if args.all else (args.model,):
        try:
            rows.append(process(model, args))
        except Exception as exc:
            print(f'{model}: FAILED {type(exc).__name__}: {exc}', file=sys.stderr, flush=True)
            rows.append({'model': model, 'status': 'Failed',
                         'error_type': type(exc).__name__, 'error': str(exc)})
    save(args.output / 'summary.json', {
        'scope': 'improved_uv_eft1_plus_mixed_scalar_mass_basis_only',
        'physical_predictions_certified': False,
        'direct_C5_matching_completed': False,
        'models': rows,
    })
    count = sum(row['status'].startswith('Success') for row in rows)
    print(f'B/C general running + mass basis: {count}/{len(rows)}; direct C5 matching NOT performed')
    return 0 if count == len(rows) else 1


if __name__ == '__main__':
    raise SystemExit(main())
