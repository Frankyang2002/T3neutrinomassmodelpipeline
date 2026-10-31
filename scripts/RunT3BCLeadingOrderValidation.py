"""Run the existing Matchete direct-vs-sequential comparison, fail closed.

This deliberately does not claim validation of an RG-evolved physical mass spectrum.
It uses the existing repository's Workstream C v22 evaluator and source artifacts.
"""
from __future__ import annotations
import argparse
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    args = ap.parse_args()
    root = args.root.resolve()
    runner = root/'scripts'/'DirectVsSequentialOneLoop.experimental.py'
    destination = root/'Reports'/'validation'/'bc_leading_order_validation_result.json'
    native_output = root/'output'/'audits'/'t3_bc_direct_vs_sequential_full_mass_v22'/'summary.json'
    destination.parent.mkdir(parents=True, exist_ok=True)
    summary = {
        'format': 'bc-leading-order-validation-gate-v1',
        'utc_time': datetime.now(timezone.utc).isoformat(),
        'source': str(runner.relative_to(root)),
        'native_output': str(native_output.relative_to(root)),
        'matchete_direct_vs_sequential_test': 'NOT_RUN',
        'physical_rg_evolved_mass_certified': False,
        'full_renormalized_threshold_matching_certified': False,
        'notes': ['Compares frozen-parameter one-loop Matchete direct and F/Y/X sequential results.',
                  'Not an external Feynman-rule derivation or validation of the running five-mass numerical spectra.'],
    }
    failures = []
    if not runner.is_file():
        failures.append('Required existing repository script is absent: '+str(runner))
    else:
        for label, argv in [
            ('selftest', [sys.executable, str(runner), '--self-test']),
            ('comparison', [sys.executable, str(runner), '--root', str(root), '--output', str(native_output)]),
        ]:
            if failures:
                break
            process = subprocess.run(argv, cwd=root, capture_output=True, text=True, errors='replace')
            summary[label] = {'exit_code': process.returncode,
                              'stdout_tail': process.stdout[-10000:],
                              'stderr_tail': process.stderr[-10000:]}
            if process.returncode:
                failures.append(f'{label} returned exit code {process.returncode}')
        if not failures:
            try:
                native = json.loads(native_output.read_text(encoding='utf-8-sig'))
                summary['native_status'] = native.get('status')
                summary['model_results'] = {
                    k:{'pass':v.get('direct_vs_sequential_full_mass_agreement'),
                       'max_relative_residual':v.get('direct_vs_sequential_max_relative_residual'),
                       'matrix_relative_residual':v.get('direct_full_3x3_C5_vs_v18_relative_residual')}
                    for k,v in native.get('models',{}).items()
                }
                passed = (native.get('all_direct_vs_sequential_full_mass_agree') is True
                          and len(native.get('models',{})) == 2
                          and all(v.get('direct_vs_sequential_full_mass_agreement') is True for v in native['models'].values()))
                if not passed:
                    failures.append('Native comparison did not confirm both models passing')
            except (OSError, ValueError, TypeError) as exc:
                failures.append(f'Cannot parse native comparison report: {exc}')
    summary['matchete_direct_vs_sequential_test'] = 'FAIL_OR_BLOCKED' if failures else 'PASS_LIMITED_SCOPE'
    summary['failures'] = failures
    destination.write_text(json.dumps(summary, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')
    print(json.dumps({'result':summary['matchete_direct_vs_sequential_test'],
                      'report':str(destination),'failures':failures,
                      'physical_rg_evolved_mass_certified':False}, indent=2))
    return 1 if failures else 0


if __name__ == '__main__':
    raise SystemExit(main())
