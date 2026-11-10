#!/usr/bin/env python3
"""Read-only audit of the selected original 16-model T3 thesis dataset.

Only --write-lock writes a provenance JSON beside the dataset; no physics is run.
Execute from the project root: python -X utf8 scripts/LockT3ThesisDataset.py
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'Reports/output/thesis_qualitative/representation_running'
DETAIL = BASE / 'csv/frozen_stage_observable_changes.csv'
MASS = BASE / 'csv/oscillation_heatmap_masssq_table.csv'
MIXING = BASE / 'csv/oscillation_heatmap_mixing_table.csv'
LOCK = BASE / 'THESIS_DATASET_LOCK.json'
MODELS = ('A_-4', 'A_-2', 'A_+0', 'B_-3', 'B_-1', 'B_+1',
          'C_-3', 'C_-1', 'C_+1', 'D_-2', 'D_+0', 'D_+2',
          'E_-4', 'E_-2', 'E_+0', 'E_+2')
STAGES = {'UV': '0111', 'EFT1': '1011', 'SMEFT': '1110'}
OBSERVABLES = ('m_light_eV', 'm_middle_eV', 'm_heavy_eV',
               'delta_m21_sq_ev2', 'delta_m3l_sq_ev2',
               'sin2_theta12', 'sin2_theta13', 'sin2_theta23')
COLUMN_OBS = {
    MASS: ('delta_m21_sq_ev2', 'delta_m3l_sq_ev2'),
    MIXING: ('sin2_theta12', 'sin2_theta13', 'sin2_theta23'),
}


def fail(message: str) -> None:
    raise RuntimeError('Thesis data lock audit FAILED: ' + message)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        fail('Missing required file: ' + str(path.relative_to(ROOT)))
    with path.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))


def number(s: str, description: str) -> float:
    try:
        x = float(s)
    except (ValueError, TypeError):
        fail('Invalid number at ' + description + ': ' + repr(s))
    if not math.isfinite(x):
        fail('Nonfinite number at ' + description)
    return x


def audit() -> dict:
    data = rows(DETAIL)
    expected = {'model', 'switch', 'on_case', 'off_case', 'observable',
                'full_value', 'frozen_value', 'relative_change_percent',
                'comparison_kind'}
    if not data or not expected.issubset(data[0]):
        fail('Detailed CSV has missing columns or rows')
    if len(data) != len(MODELS) * 4 * len(OBSERVABLES):
        fail(f'Expected {len(MODELS) * 4 * len(OBSERVABLES)} detailed rows; got {len(data)}')
    index = {}
    counts = Counter()
    missing_counts = Counter()
    for r in data:
        model, stage, observable = r['model'], r['switch'], r['observable']
        key = (model, stage, observable)
        if key in index:
            fail('Duplicate row ' + str(key))
        if model not in MODELS or observable not in OBSERVABLES:
            fail('Unexpected model or observable: ' + str(key))
        expected_off = STAGES.get(stage, '1101' if stage == 'Wilson_conditional' else None)
        if expected_off is None or r['on_case'] != '1111' or r['off_case'] != expected_off:
            fail('Wrong stage switch for ' + str(key))
        expected_kind = ('conditional_wilson_log_reconstruction' if stage == 'Wilson_conditional'
                         else 'independently_evolved_running_stage')
        if r['comparison_kind'] != expected_kind:
            fail('Incorrect comparison kind ' + str(key))
        # The original CSV contains genuinely unavailable entries in the
        # conditional Wilson reconstruction.  They are missing data, NOT zeros.
        # The three plotted running stages must still have finite values.
        required = ('full_value', 'frozen_value', 'relative_change_percent')
        missing = [name for name in required if r[name] is None or not r[name].strip()]
        counts[stage] += 1
        if missing:
            if stage != 'Wilson_conditional':
                fail('Missing plotted value(s) ' + str(missing) + ' for ' + str(key))
            # Nonempty fields must nevertheless be valid finite numbers.
            for name in required:
                if name not in missing:
                    number(r[name], str(key) + ' ' + name)
            missing_counts[stage] += 1
            index[key] = None
            continue
        full = number(r['full_value'], str(key) + ' full')
        off = number(r['frozen_value'], str(key) + ' off')
        pct = number(r['relative_change_percent'], str(key) + ' pct')
        if off == 0:
            fail('Zero denominator in ' + str(key))
        computed = 100 * (full - off) / off
        if not math.isclose(pct, computed, abs_tol=1e-8, rel_tol=1e-9):
            fail('Percentage arithmetic mismatch ' + str(key))
        index[key] = pct
    for stage in (*STAGES, 'Wilson_conditional'):
        if counts[stage] != 16 * len(OBSERVABLES):
            fail(f'Incomplete {stage} stage: {counts[stage]}')
    for path, obs_list in COLUMN_OBS.items():
        matrix = rows(path)
        if len(matrix) != len(MODELS):
            fail('Heatmap matrix not 16 models: ' + str(path))
        cols = list(matrix[0])
        if len(cols) != 1 + 3 * len(obs_list):
            fail('Heatmap matrix column count changed: ' + str(path))
        for model, r in zip(MODELS, matrix):
            if r[cols[0]] != model:
                fail('Heatmap model ordering mismatch: ' + model)
            for j, stage in enumerate(STAGES):
                for k, observable in enumerate(obs_list):
                    stored = number(r[cols[1 + j * len(obs_list) + k]],
                                    f'{path.name} {model} {stage} {observable}')
                    expected_pct = index[(model, stage, observable)]
                    # Plotting CSV is rounded to roughly 10 significant figures.
                    if not math.isclose(stored, expected_pct, abs_tol=1e-6,
                                        rel_tol=5e-8):
                        fail(f'Heatmap does not match detailed CSV: {model} {stage} {observable}: '
                             f'{stored} != {expected_pct}')
    # One independently known checkpoint from the original published source CSV.
    if not math.isclose(index[('A_-4', 'UV', 'delta_m21_sq_ev2')], 15.954774428199116,
                        abs_tol=1e-8):
        fail('Original A(-4) UV checkpoint differs: wrong dataset?')
    return {
        'schema': 't3-original16-thesis-dataset-lock-v1',
        'scientific_status': 'original_exploratory_conditional_not_physical_matching_certification',
        'source': 'original_16_model_representation_running',
        'producer': 'unidentified_upstream_driver; plotting-only script verified',
        'canonical_dataset_root': str(BASE.relative_to(ROOT)).replace('\\', '/'),
        'full_running_case': '1111',
        'stage_off': {**STAGES, 'Wilson_conditional': '1101'},
        'wilson_comparison_kind': 'conditional_wilson_log_reconstruction',
        'relative_percent_definition': '100*(O_full-O_stage_off)/O_stage_off',
        'models': list(MODELS), 'detail_rows': len(data),
        'verified_stage_rows': dict(sorted(counts.items())),
        'unavailable_conditional_wilson_rows': missing_counts['Wilson_conditional'],
        'file_sha256': {
            str(p.relative_to(ROOT)).replace('\\', '/'): digest(p)
            for p in (DETAIL, MASS, MIXING)
        },
        'do_not_mix_with': [
            'scripts/RunT3ThesisTwoModelPilot.py',
            'scripts/RunT3GeneralisedRunningComparisons.py',
            'output/parallel_t3/thesis_original_style_figures_v2',
            'output/parallel_t3/thesis_original_style_figures_v3',
        ],
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--write-lock', action='store_true',
                    help='Write a JSON SHA256 manifest alongside the ORIGINAL results')
    ap.add_argument('--check-lock', action='store_true',
                    help='Also verify hashes against the saved lock JSON')
    args = ap.parse_args()
    manifest = audit()
    if args.check_lock:
        if not LOCK.is_file():
            fail('No lock exists yet; first run with --write-lock')
        previous = json.loads(LOCK.read_text(encoding='utf-8'))
        if previous.get('file_sha256') != manifest['file_sha256']:
            fail('Selected source file hashes have changed since lock creation')
    if args.write_lock:
        LOCK.write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
        print('Saved:', LOCK.relative_to(ROOT))
    print('PASS: original 16-model dataset, 512 rows, 3 plotted running switches,')
    print('      unavailable conditional Wilson entries:',
          manifest['unavailable_conditional_wilson_rows'], '(kept as missing, not zero)')
    print('      separate conditional Wilson switch, percentage definitions,')
    print('      both original heatmap matrices, and A(-4) checkpoint.')
    for path, hashvalue in manifest['file_sha256'].items():
        print('SHA256', hashvalue, path)
    print('NOTE: dataset provenance audit only; no EFT calculation or physics certification.')


if __name__ == '__main__':
    main()
