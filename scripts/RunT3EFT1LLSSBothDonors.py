#!/usr/bin/env python3
"""Create matching, title-free A(-4) and D(-2) EFT1 LLSS heatmap reports.

Project path: scripts/RunT3EFT1LLSSBothDonors.py

This wrapper calls the unchanged three-way physics runner separately per donor.
Use --a4-source/--d2-source to reuse completed saved results (replot only).
For donors without an input source it recomputes all three Wilson cases for
all 16 models, using the selected donor's own fitted Yukawas. No font override
is made: the plotter uses Matplotlib defaults, with user-controlled sizes.
"""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / 'scripts' / 'RunT3EFT1LLSSIsolation.py'
OUTPUT_ROOT = ROOT / 'output' / 'parallel_t3'
DONORS = (('A_-4', 'A_m4', 'a4_source'), ('D_-2', 'D_m2', 'd2_source'))
FONT_OPTIONS = {
    'axis_label_fontsize': 12.0,
    'bottom_tick_fontsize': 10.0,
    'top_tick_fontsize': 13.0,
    'y_tick_fontsize': 11.0,
    'cell_fontsize': 8.0,
    'colorbar_fontsize': 10.0,
}


def _existing_result(path: Path, donor: str) -> None:
    csv_path = path / 'all16_three_way.csv'
    if not csv_path.is_file():
        raise ValueError(f'{donor}: missing saved three-way data: {csv_path}')
    with csv_path.open(newline='', encoding='utf-8-sig') as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 16 or len({r.get('model') for r in rows}) != 16:
        raise ValueError(f'{donor}: expected exactly 16 distinct saved models, got {len(rows)}')
    for row in rows:
        if row.get('reference') and row['reference'] != donor:
            raise ValueError(f'{donor}: source CSV contains a different Yukawa donor')
    meta_path = path / 'provenance.json'
    if meta_path.is_file():
        meta = json.loads(meta_path.read_text(encoding='utf-8-sig'))
        if meta.get('reference') not in (None, donor):
            raise ValueError(f'{donor}: provenance describes {meta["reference"]}, not {donor}')


def _style_args(args: argparse.Namespace) -> list[str]:
    options = []
    for key, default in FONT_OPTIONS.items():
        value = getattr(args, key)
        value = args.font_size if value is None and args.font_size is not None else (default if value is None else value)
        if not 0 < float(value) <= 72:
            raise ValueError(f'Invalid font size {key}={value}')
        options += ['--' + key.replace('_', '-'), str(value)]
    if args.figure_width <= 0 or args.row_height <= 0 or not 72 <= args.dpi <= 1200:
        raise ValueError('Invalid figure dimensions or DPI')
    options += ['--figure-width', str(args.figure_width), '--row-height', str(args.row_height), '--dpi', str(args.dpi)]
    return options


def build_jobs(args: argparse.Namespace) -> list[dict]:
    if not RUNNER.is_file():
        raise FileNotFoundError(f'Missing existing runner {RUNNER}')
    out = (ROOT / args.output).resolve()
    if out == OUTPUT_ROOT.resolve() or OUTPUT_ROOT.resolve() not in out.parents:
        raise ValueError('--output must be below output/parallel_t3')
    if out.exists():
        raise ValueError(f'Output already exists: {out}; choose a new --output')
    style = _style_args(args)
    jobs = []
    for donor, tag, attr in DONORS:
        src = getattr(args, attr)
        source = (ROOT / src).resolve() if src is not None else None
        if source is not None:
            _existing_result(source, donor)
        dest = out / tag
        argv = [sys.executable, str(RUNNER)]
        if source is not None:
            argv += ['--plot-only', str(source)]
        else:
            argv += ['--reference', donor]
            if args.reuse_cached_batch1:
                argv.append('--reuse-cached-batch1')
            if args.debug:
                argv.append('--debug')
        argv += ['--output', str(dest), *style]
        jobs.append({'donor': donor, 'source': str(source) if source else None,
                     'directory': str(dest), 'command': argv,
                     'mode': 'replot_saved' if source else 'run_three_way'})
    return jobs


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--a4-source', type=Path, help='Existing completed A(-4) all16_three_way.csv directory')
    parser.add_argument('--d2-source', type=Path, help='Existing completed D(-2) all16_three_way.csv directory')
    parser.add_argument('--output', type=Path, required=True,
                        help='NEW parent directory containing A_m4 and D_m2 report subdirectories')
    parser.add_argument('--reuse-cached-batch1', action='store_true',
                        help='For a donor without saved results, reuse hashed cached Batch 1 outputs')
    parser.add_argument('--font-size', type=float,
                        help='Set all text sizes together; specific font-size flags override this')
    for key in FONT_OPTIONS:
        parser.add_argument('--' + key.replace('_', '-'), type=float, default=None)
    parser.add_argument('--figure-width', type=float, default=16.0)
    parser.add_argument('--row-height', type=float, default=0.48)
    parser.add_argument('--dpi', type=int, default=300)
    parser.add_argument('--debug', action='store_true')
    args = parser.parse_args(argv)
    try:
        jobs = build_jobs(args)
    except (ValueError, FileNotFoundError) as exc:
        parser.error(str(exc))
    out = (ROOT / args.output).resolve()
    out.mkdir(parents=True)
    summary = {'schema': 't3-eft1-llss-both-donors-report-v1',
               'figure_title': None, 'font_family': 'Matplotlib default',
               'common_frozen_denominator': True,
               'comparison_running': '100*(O_full-O_frozen)/O_frozen',
               'comparison_mixing_off': '100*(O_LLSS_off-O_frozen)/O_frozen',
               'jobs': []}
    success = True
    for job in jobs:
        print(f'\n{job["donor"]} donor: {job["mode"]}', flush=True)
        process = subprocess.run(job['command'], cwd=ROOT, check=False)
        record = {k: v for k, v in job.items() if k != 'command'}
        record['exit_code'] = int(process.returncode)
        summary['jobs'].append(record)
        (out / 'both_donors_summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
        if process.returncode:
            success = False
            print(f'{job["donor"]}: calculation/report returned {process.returncode}', flush=True)
    for job in summary['jobs']:
        report = Path(job['directory']) / 'reports'
        print(f'{job["donor"]}: figures={Path(job["directory"]) / "plots"}, report={report}', flush=True)
    return 0 if success else 1


if __name__ == '__main__':
    raise SystemExit(main())
