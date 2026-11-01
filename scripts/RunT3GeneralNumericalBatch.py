"""Run the existing full-general symbolic and UV/EFT1 numerical pipeline per model.

No matching, physics expressions, or numerical integrators are reimplemented.
A_-2 is excluded by default because it has already been completed separately.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import subprocess
import sys
import traceback

ROOT = Path(__file__).resolve().parents[1]
MODELS = {
    'A': (-4, -2, 0),
    'B': (-3, -1, 1),
    'C': (-3, -1, 1),
    'D': (-2, 0, 2),
    'E': (-4, -2, 0, 2),
}
ALL_MODELS = tuple(f'{cls}_{a:+d}' for cls, alphas in MODELS.items() for a in alphas)


def run_command(command: list[str], log: Path) -> tuple[bool, int, str]:
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open('w', encoding='utf-8') as handle:
        handle.write('COMMAND: ' + subprocess.list2cmdline(command) + '\n\n')
        handle.flush()
        try:
            process = subprocess.run(command, cwd=ROOT, stdout=handle,
                                     stderr=subprocess.STDOUT, check=False)
            return process.returncode == 0, process.returncode, ''
        except OSError as exc:
            handle.write('\nLAUNCH ERROR: ' + str(exc) + '\n')
            return False, -1, str(exc)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--include-a-m2', action='store_true',
                   help='Also run A_-2 using a fresh batch-specific run ID')
    p.add_argument('--model', choices=ALL_MODELS, action='append',
                   help='Run only specified model(s); can be repeated')
    p.add_argument('--run-prefix', default='all16_general_numerical_v1',
                   help='Prefix for per-model canonical run IDs (max 35 chars)')
    p.add_argument('--output-dir', type=Path,
                   default=Path('output/general_numerical_batch_v1'))
    p.add_argument('--mu-uv', type=float, default=1e7)
    p.add_argument('--mu-f', type=float, default=1e5)
    p.add_argument('--mu-s', type=float, default=1e3)
    p.add_argument('--initial-json', type=Path,
                   help='Optional shared benchmark coupling overrides')
    p.add_argument('--wolframscript', default='wolframscript')
    p.add_argument('--list', action='store_true', help='Print selected models and exit')
    args = p.parse_args(argv)
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,34}', args.run_prefix):
        p.error('--run-prefix must be 1-35 letters/digits/underscores/hyphens')
    if not args.mu_uv > args.mu_f > args.mu_s > 0:
        p.error('Require --mu-uv > --mu-f > --mu-s > 0')
    if args.initial_json and not (ROOT / args.initial_json).is_file():
        p.error('Initial JSON does not exist')
    selected = list(dict.fromkeys(args.model)) if args.model else list(ALL_MODELS)
    if not args.model and not args.include_a_m2:
        selected.remove('A_-2')
    if args.list:
        print('\n'.join(selected))
        return 0
    output_dir = (ROOT / args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_file = output_dir / 'summary.json'
    rows: list[dict] = []
    for index, model in enumerate(selected, 1):
        token = model.replace('+', 'p').replace('-', 'm')
        run_id = f'{args.run_prefix}_{token}'
        study = f'general_batch_{token}'
        output_json = output_dir / f'{token}_running.json'
        symbolic_log = output_dir / 'logs' / f'{token}_symbolic.log'
        numerical_log = output_dir / 'logs' / f'{token}_numerical.log'
        row = {'model': model, 'run_id': run_id, 'symbolic_log': str(symbolic_log),
               'numerical_log': str(numerical_log), 'output_json': str(output_json)}
        print(f'[{index}/{len(selected)}] {model}: symbolic export', flush=True)
        symbolic = [sys.executable, '-X', 'utf8', '-m',
                    'scripts.RunT3GeneralSymbolicSafe', '--model', model,
                    '--run-id', run_id, '--study', study,
                    '--wolframscript', args.wolframscript]
        success, code, error = run_command(symbolic, symbolic_log)
        row['symbolic_exit_code'] = code
        if not success:
            row['status'] = 'SymbolicFailed'
            if error:
                row['error'] = error
        else:
            print(f'[{index}/{len(selected)}] {model}: numerical UV/EFT1', flush=True)
            numerical = [sys.executable, '-X', 'utf8', '-m',
                         'Numerical.orchestration.GeneralRunningBenchmark',
                         '--canonical-run-id', run_id, '--model', model,
                         '--output-json', str(output_json), '--mu-uv', str(args.mu_uv),
                         '--mu-f', str(args.mu_f), '--mu-s', str(args.mu_s)]
            if args.initial_json:
                numerical += ['--initial-json', str((ROOT / args.initial_json).resolve())]
            success, code, error = run_command(numerical, numerical_log)
            row['numerical_exit_code'] = code
            if success and output_json.is_file():
                try:
                    payload = json.loads(output_json.read_text(encoding='utf-8'))
                    if (payload.get('model') != model
                        or payload.get('status') != 'NumericalRenormalisableRunningOnly'
                        or not payload.get('uv_solver', {}).get('success')
                        or not payload.get('eft1_solver', {}).get('success')):
                        raise ValueError('Numerical result status/model/solvers inconsistent')
                    row['status'] = 'Completed'
                    row['uv_nfev'] = payload['uv_solver']['nfev']
                    row['eft1_nfev'] = payload['eft1_solver']['nfev']
                except (ValueError, OSError, KeyError, TypeError) as exc:
                    row['status'] = 'OutputInvalid'
                    row['error'] = str(exc)
            else:
                row['status'] = 'NumericalFailed' if not success else 'OutputMissing'
                if error:
                    row['error'] = error
        rows.append(row)
        summary_file.write_text(json.dumps({
            'generated_utc': datetime.now(timezone.utc).isoformat(),
            'scope': 'renormalisable_UV_and_EFT1_only',
            'selected_models': selected,
            'completed': sum(r['status'] == 'Completed' for r in rows),
            'processed': len(rows), 'models': rows,
        }, indent=2) + '\n', encoding='utf-8')
        print(f'[{index}/{len(selected)}] {model}: {row["status"]}', flush=True)
    completed = sum(r['status'] == 'Completed' for r in rows)
    print(f'Completed {completed}/{len(rows)}. Summary: {summary_file}')
    return 0 if completed == len(rows) else 1


if __name__ == '__main__':
    raise SystemExit(main())
