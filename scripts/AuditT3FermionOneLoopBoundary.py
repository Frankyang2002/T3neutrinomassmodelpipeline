"""Read-only Worker A audit: locate fermion Matchete evidence and test tensor interface.

Run from repository root. Produces a NEW JSON file (exclusive create) only
if --output is supplied. Does not change solver, inputs, or historical outputs.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from Numerical.running.T3FermionOneLoopBoundary import (
    assemble_fermion_boundary, UncertifiedMatchingError,
)

PATTERNS = ('**/native_status.json', '**/one_loop_operator_ledger.json',
            '**/laurent_audit.json', '**/eft1_after_F_one_loop_inputform.txt')


def probe(root: Path, limit: int) -> dict:
    results = []
    for pattern in PATTERNS:
        matches = sorted(root.glob(pattern))
        for file in matches[:limit]:
            record = {'relative_path': str(file.relative_to(root)), 'bytes': file.stat().st_size}
            if file.suffix.lower() == '.json':
                try:
                    obj = json.loads(file.read_text(encoding='utf-8-sig'))
                    if isinstance(obj, dict):
                        for name in ('status', 'model', 'integrated_fields',
                                     'light_bsm_fields', 'numerical_finite_Wilson_matching_complete',
                                     'pole_cancelled', 'one_loop_MSbar_renormalized'):
                            if name in obj:
                                record[name] = obj[name]
                except (UnicodeError, ValueError) as exc:
                    record['parse_error'] = str(exc)
            results.append(record)
    return {'repository_root': str(root), 'evidence': results,
            'note': 'Discovery is not certification; inspect original Matchete exports and EFT1 soft subtraction.'}


def smoke() -> dict:
    # Synthetic tensor: checks interface shape, symmetry and guard, NOT physics.
    t = np.zeros((3,3,2,2,2,2), dtype=complex)
    t[0,0,0,0,0,0] = 1e-5
    parts = assemble_fermion_boundary('A_-4', t)
    assert parts.finite_one_loop is None
    try:
        _ = parts.total
    except UncertifiedMatchingError:
        blocked = True
    else:
        blocked = False
    assert blocked
    try:
        assemble_fermion_boundary('D_-2', t, finite_one_loop=t*0.01,
                                  provenance={'source_type':'raw_matchete_laurent'})
    except UncertifiedMatchingError:
        rejected = True
    else:
        rejected = False
    assert rejected
    return {'synthetic_tree_norm_gev_inverse':float(np.linalg.norm(parts.tree)),
            'uncertified_total_blocked':blocked,
            'raw_D_minus2_rejected':rejected,
            'status':'interface_smoke_only_not_matched_benchmark'}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, default=ROOT)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--max-per-pattern', type=int, default=30)
    args = parser.parse_args()
    if args.max_per_pattern <= 0:
        parser.error('--max-per-pattern must be positive')
    root = args.repo.resolve()
    if not root.is_dir():
        parser.error('repository directory does not exist')
    result = {'smoke':smoke(), 'source_discovery':probe(root,args.max_per_pattern)}
    serialized = json.dumps(result, indent=2, allow_nan=False)+'\n'
    if args.output:
        path = args.output.resolve()
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open('x', encoding='utf-8') as handle:
            handle.write(serialized)
    print(serialized)

if __name__ == '__main__':
    main()
