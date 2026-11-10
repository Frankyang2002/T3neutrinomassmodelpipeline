"""Plot new generalised running comparisons with the ORIGINAL thesis heatmap maker.

Converts RunT3GeneralisedRunningComparisons.py's wide running_comparisons.csv
into the long CSV accepted by scripts/PlotT3OscillationHeatmaps.py.
No physics is recalculated, and the original plotting code is not modified.
Failed runs become BLANK cells (NaN), NEVER numerical zero; details saved in JSON.

Intended project path: scripts/PlotT3GeneralisedOriginalHeatmaps.py
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
PLOTTER = ROOT / 'scripts' / 'PlotT3OscillationHeatmaps.py'
STAGES = ('UV', 'EFT1', 'SMEFT')
OBSERVABLES = {
    'sin2_theta12': 'sin2_theta12',
    'sin2_theta13': 'sin2_theta13',
    'sin2_theta23': 'sin2_theta23',
    'dm21_sq_eV2': 'delta_m21_sq_ev2',
    'dm31_sq_eV2': 'delta_m3l_sq_ev2',
}
MODELS = (
    'A_-4', 'A_-2', 'A_+0',
    'B_-3', 'B_-1', 'B_+1',
    'C_-3', 'C_-1', 'C_+1',
    'D_-2', 'D_+0', 'D_+2',
    'E_-4', 'E_-2', 'E_+0', 'E_+2',
)


def convert(source: Path, destination: Path, reference: str) -> tuple[int, list[dict]]:
    """Produce original plotter input while leaving missing data unrepresented.

    The repository's original plotter populates absent cells with NaN and draws
    those cells blank. Missing/blocked runs are never assigned numerical zero.
    """
    if not source.is_file():
        raise FileNotFoundError(f'Missing computed stage-freeze CSV: {source}')
    with source.open(newline='', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        expected = {'reference', 'model', 'stage_frozen', 'status'} | {
            k + '_pct' for k in OBSERVABLES
        }
        missing = expected - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f'Missing required computed columns: {sorted(missing)}')
        raw = [r for r in reader if r['reference'] == reference]

    lookup = {}
    for row in raw:
        key = (row['model'], row['stage_frozen'])
        if key in lookup:
            raise ValueError(f'Duplicate running calculation: {key}')
        lookup[key] = row

    long_rows = []
    omitted = []
    for model in MODELS:
        for stage in STAGES:
            key = (model, stage)
            row = lookup.get(key)
            if row is None:
                omitted.append({'reference':reference,'model':model,'stage':stage,
                                'reason':'No computed counterfactual',
                                'observables':list(OBSERVABLES)})
                continue
            if row['status'] != 'VALID_EXPLORATORY':
                omitted.append({'reference':reference,'model':model,'stage':stage,
                                'reason':row.get('error') or row['status'],
                                'observables':list(OBSERVABLES)})
                continue
            for source_key,target_key in OBSERVABLES.items():
                value = (row.get(source_key+'_pct') or '').strip()
                try:
                    numeric=float(value)
                    valid=math.isfinite(numeric)
                except ValueError:
                    valid=False
                if not valid:
                    omitted.append({'reference':reference,'model':model,'stage':stage,
                                    'reason':'Undefined/nonfinite percentage',
                                    'observables':[source_key]})
                    continue
                long_rows.append({'model':model,'switch':stage,
                    'observable':target_key,
                    'relative_change_percent':repr(numeric),
                    'comparison_kind':'independently_evolved_running_stage'})
    if not long_rows:
        raise ValueError(f'No finite running comparisons for {reference}')
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=(
            'model','switch','observable','relative_change_percent','comparison_kind'))
        writer.writeheader()
        writer.writerows(long_rows)
    (destination.parent/'missing_running_cells.json').write_text(
        json.dumps(omitted,indent=2)+'\n',encoding='utf-8')
    return len(long_rows),omitted


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--running', type=Path, required=True,
                   help='Directory from RunT3GeneralisedRunningComparisons.py')
    p.add_argument('--reference', choices=('D_-2', 'A_-4', 'both'), default='D_-2')
    p.add_argument('--output', type=Path, required=True,
                   help='Fresh directory for original-style figures and converted CSV')
    args = p.parse_args()
    if not PLOTTER.is_file():
        p.error(f'Original repository figure-maker not found: {PLOTTER}')
    source = args.running / 'running_comparisons.csv'
    refs = ('D_-2', 'A_-4') if args.reference == 'both' else (args.reference,)
    out = args.output.resolve()
    if out.exists() and any(out.iterdir()):
        p.error(f'Output directory is not empty: {out}')
    for ref in refs:
        target = out / ref
        input_csv = target / 'frozen_stage_observable_changes.csv'
        count, omitted = convert(source, input_csv, ref)
        print(f'{ref}: {count} valid observable cells; {len(omitted)} missing entries, drawn blank',flush=True)
        subprocess.run((sys.executable, str(PLOTTER),
                        '--input', str(input_csv),
                        '--plot-dir', str(target / 'plots'),
                        '--csv-dir', str(target / 'csv')),
                       check=True, cwd=str(ROOT))
        print(f'Figures: {target / "plots"}', flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
