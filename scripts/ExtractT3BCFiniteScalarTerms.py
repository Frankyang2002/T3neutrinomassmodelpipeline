"""Audit pole/log/constant sectors from saved native Y/X Matchete outputs.

NO Wolfram run, NO modification of original matching data.  These exports
are formal pole-stripped candidates, NOT physical matching coefficients.
"""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path

from Numerical.experimental.T3BCMatcheteFiniteTerms import (
    AuditError, split_inputform,
)

ROOT = Path(__file__).resolve().parents[1]
MODELS = ('B_-1', 'C_-1')
Y_DIR = Path('output/general_t3_mass_eigen_matchete/eft1_scalar_Y_only_canonical_probe')
X_DIR = Path('output/general_t3_mass_eigen_matchete/eft1_scalar_xy_fresh_kernel_v3')
DEFAULT_OUTPUT = Path('output/audits/t3_bc_scalar_finite_term_extraction_v1')
FIELDS = ('model', 'threshold', 'sector', 'terms_total', 'pole_terms',
          'log_terms', 'constant_terms', 'physical_matching_certified')


def _read(path: Path) -> str:
    if not path.is_file():
        raise AuditError(f'Missing file: {path}')
    return path.read_text(encoding='utf-8')


def _json(path: Path) -> dict:
    data = json.loads(_read(path))
    if not isinstance(data, dict):
        raise AuditError(f'Not a JSON object: {path}')
    return data


def process(model: str, root: Path, y_dir: Path, x_dir: Path, out_dir: Path) -> list[dict]:
    if model not in MODELS:
        raise AuditError(f'Unexpected model: {model}')
    short = model.replace('_-', '_m')
    y_stem = f'{short}_final_y_only_canonical_match'
    x_stem = f'{short}_final_x_from_saved_y_match'
    yd, xd = root / y_dir, root / x_dir
    y_cert = _json(yd / (y_stem + '.json'))
    x_cert = _json(xd / (x_stem + '_x_followup.json'))
    if (y_cert.get('model') != model or
            y_cert.get('status') != 'YOnlyIndependentCanonicalRawEFT2NotCertified' or
            y_cert.get('tree_matching_completed') is not True or
            y_cert.get('one_loop_matching_completed') is not True or
            y_cert.get('raw_EFT2_operator_sector_counts', {}).get('LLHH', 0) <= 0):
        raise AuditError(f'{model}: Y-only source certificate missing or invalid')
    if (x_cert.get('model') != model or
            x_cert.get('x_tree_matching_completed') is not True or
            x_cert.get('rebased_check_passed') is not True or
            x_cert.get('source_is_y_tree_EOMEFT') is not True):
        raise AuditError(f'{model}: X source/handoff checks unavailable')
    log = _read(xd / (x_stem + '.log'))
    if (log.count('Matching pipeline completed successfully.') < 2 or
            'x-followup-one-loop-matching-returned' not in log or
            'x-followup-matching-exports-returned' not in log):
        raise AuditError(f'{model}: no independently evidenced completed X matching')
    # Original v3 report is left at an in-progress checkpoint due to a
    # verified post-export AssociateTo bug.  It must NOT be treated as success
    # without the separate Matchete log and the complete operator exports.
    selectors = (
        ('Y', 'LLHH', yd / (y_stem + '_raw_llhh.inputform.txt')),
        ('Y', 'LLHX', yd / (y_stem + '_raw_llhx.inputform.txt')),
        ('Y', 'LLXX', yd / (y_stem + '_raw_llxx.inputform.txt')),
        ('X', 'LLHH', xd / (x_stem + '_x_followup_x_raw_llhh.inputform.txt')),
    )
    staged: list[tuple[object, Path, dict]] = []
    for stage, sector, path in selectors:
        source = _read(path)
        split = split_inputform(source, stage=stage, sector=sector)
        expected = y_cert.get('raw_EFT2_operator_sector_counts', {}).get(sector) if stage == 'Y' else 36
        if expected is None or len(split.terms) != expected:
            raise AuditError(f'{model} {stage}/{sector}: expected {expected}, found {len(split.terms)}')
        row = {'model': model, **split.summary()}
        row['source_file'] = str(path.relative_to(root)).replace('\\', '/')
        staged.append((split, path, row))
    # Fail closed *before* writing any output for this model.
    out_model = out_dir / short
    out_model.mkdir(parents=True, exist_ok=True)
    for split, path, row in staged:
        stem = f'{split.stage}_{split.sector}'
        files = {
            'pole_residue': split.pole_residue,
            'pole_free_unrenormalized': split.pole_free,
            'at_mu_equal_heavy_mass_unrenormalized': split.scale_equal_mass,
        }
        for name, expr in files.items():
            target = out_model / f'{stem}_{name}.inputform.txt'
            target.write_text(expr + '\n', encoding='utf-8')
            row[name + '_file'] = str(target.relative_to(out_dir)).replace('\\', '/')
        (out_model / f'{stem}_audit.json').write_text(
            json.dumps(row, indent=2, ensure_ascii=False) + '\n', encoding='utf-8'
        )
    return [r for _, _, r in staged]


def run(root: Path, models: tuple[str, ...], y_dir: Path, x_dir: Path,
        output: Path) -> list[dict]:
    rows = []
    # No source writes; output kept under audit directory.
    for model in models:
        model_rows = process(model, root, y_dir, x_dir, output)
        rows.extend(model_rows)
        for row in model_rows:
            print(f'{model} {row["threshold"]}/{row["sector"]}: '
                  f'pole={row["pole_terms"]}, log={row["log_terms"]}, '
                  f'constant={row["constant_terms"]}, '
                  'physical=NOT CERTIFIED', flush=True)
    with (output / 'summary.csv').open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    summary = {
        'format': 't3-bc-native-local-matchete-pole-log-split-v1',
        'status': 'SymbolicLocalAmplitudePartitionedNoRenormalizedPhysicalMatching',
        'models': list(models),
        'sector_results': rows,
        'one_loop_y_inherited_into_x_calculation': False,
        'original_fermion_first_eomeft_equivalence_proven': False,
        'genuine_matching_counterterms_calculated': False,
        'numerical_running_boundary_substitution_performed': False,
        'physical_sequential_matching_certified': False,
    }
    (output / 'summary.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    return rows


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument('--all', action='store_true')
    g.add_argument('--model', choices=MODELS)
    p.add_argument('--root', type=Path, default=ROOT)
    p.add_argument('--y-dir', type=Path, default=Y_DIR)
    p.add_argument('--x-dir', type=Path, default=X_DIR)
    p.add_argument('--output-dir', type=Path, default=DEFAULT_OUTPUT)
    args = p.parse_args()
    root = args.root.resolve()
    output = (args.output_dir if args.output_dir.is_absolute() else root / args.output_dir).resolve()
    if not str(output).startswith(str(root)):
        raise AuditError('Output must be under the project root')
    run(root, MODELS if args.all else (args.model,), args.y_dir, args.x_dir, output)
    print('Finite-term audit summary:', output / 'summary.csv')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
