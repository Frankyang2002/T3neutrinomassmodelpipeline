"""Audit Batch14 EFT1 endpoints at the simultaneous scalar threshold.

Project path: scripts/RunBatch15ScalarThresholdAll16.py
Usage from repository root:
    python scripts/RunBatch15ScalarThresholdAll16.py --model A_-4
    python scripts/RunBatch15ScalarThresholdAll16.py

This produces SM *renormalisable* boundaries; NO Wilson C5, SMEFT running,
physical neutrino masses, or completed scalar matching is claimed.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import traceback

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Numerical.experimental.Batch15ScalarThresholdHandoff import (
    handoff, read_json, sha256_file,
)
from Numerical.experimental.Batch13GeneralUVRunning import REPRESENTATIONS

MODELS = ('A_+0','A_-2','A_-4','B_+1','B_-1','B_-3',
          'C_+1','C_-1','C_-3','D_+0','D_+2','D_-2',
          'E_+0','E_+2','E_-2','E_-4')


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--model', choices=MODELS)
    ap.add_argument('--uv-exports', type=Path, default=ROOT/'output/batch13/all16_uv')
    ap.add_argument('--uv-endpoints', type=Path, default=ROOT/'output/batch13/general_uv_running')
    ap.add_argument('--eft1', type=Path, default=ROOT/'output/batch14/eft1')
    ap.add_argument('--out', type=Path, default=ROOT/'output/batch15/scalar_threshold')
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    models = (args.model,) if args.model else MODELS
    summary = {
        'format': 't3-batch15-scalar-threshold-audit-summary-v1',
        'status': 'RenormalisableHandoffOnly',
        'physical_neutrino_mass_prediction': False,
        'models': {},
    }
    for model in models:
        try:
            uv_beta_path = args.uv_exports/f'{model}_UV.json'
            uv_final_path = args.uv_endpoints/f'{model}_UV_running.json'
            eft_beta_path = args.eft1/f'{model}_EFT1.json'
            eft_final_path = args.eft1/f'{model}_EFT1_running.json'
            uv_beta = read_json(uv_beta_path)
            eft_beta = read_json(eft_beta_path)
            eft_final = read_json(eft_final_path)
            report = handoff(model, uv_beta, eft_beta, uv_final_path, eft_final)
            report['input_sha256'] = {
                'uv_beta': sha256_file(uv_beta_path),
                'uv_endpoint': sha256_file(uv_final_path),
                'eft1_beta': sha256_file(eft_beta_path),
                'eft1_endpoint': sha256_file(eft_final_path),
            }
            (args.out/f'{model}_scalar_threshold_handoff.json').write_text(
                json.dumps(report, indent=2, allow_nan=False), encoding='utf-8')
            quad = report['scalar_threshold']['additional_quadratic_symbols']
            proxy_ok = report['scalar_threshold']['diagonal_proxy_threshold_hierarchy_passed']
            print(f'[{model}] AUDITED | {report["n_eft1_couplings"]} EFT1 couplings | '
                  f'quadratic={quad or "none"} | proxy hierarchy={proxy_ok} | physical C5=BLOCKED',
                  flush=True)
            summary['models'][model] = {
                'passed': True, 'physical_c5_blocked': True,
                'extra_quadratic': bool(quad), 'diagonal_proxy_hierarchy_passed': proxy_ok,
            }
        except Exception as exc:
            summary['models'][model] = {'passed': False, 'error': str(exc)}
            (args.out/f'{model}_error.log').write_text(traceback.format_exc(), encoding='utf-8')
            print(f'[{model}] FAILED: {exc}', flush=True)
    summary['passed'] = sum(bool(v['passed']) for v in summary['models'].values())
    summary['failed'] = len(models)-summary['passed']
    summary['with_additional_quadratic'] = sum(v.get('extra_quadratic',False)
                                                for v in summary['models'].values())
    (args.out/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(f'Summary: {args.out/"summary.json"} | audited={summary["passed"]} failed={summary["failed"]}')
    return 0 if summary['failed']==0 else 1


if __name__ == '__main__':
    raise SystemExit(main())
