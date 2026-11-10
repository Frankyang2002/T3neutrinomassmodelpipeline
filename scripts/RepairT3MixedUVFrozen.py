#!/usr/bin/env python3
"""Recover B_-1/C_-1 UV-frozen counterfactuals without inventing Wilson terms.

Source-specific diagnosis: the saved mixed Matchete expansion contains terms
whose mixing spurion prefactor can vanish at the frozen UV benchmark. The
repository's project_ordered_term previously rejected an empty SINGLE term.
Here only an independently verified zero prefactor is allowed to return an
empty component list, which the existing _dense_operator represents as a
zero tensor. The complete boundary must still be nonzero.

This is a local runtime adapter for the exploratory numerical pipeline; the
repository physics files, historical results and GitHub remain unchanged.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import sys
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Numerical.running import GeneralT3OrderedTensor as ordered
from Numerical.running import GeneralT3MixedExactFermionWilson as mixed
from scripts import RunT3GeneralisedRunningComparisons as running

MODELS = ('B_-1', 'C_-1')
ORIGINAL_PROJECT = ordered.project_ordered_term
ORIGINAL_MIXED = mixed.mixed_exact_boundary
ZERO_PREF_TOL = 1e-15  # Dimensionless prefactor; matches component cutoff.
ZERO_TERMS: list[dict] = []


def project_allow_structural_zero(term, registry, model, *, max_entries=3000):
    try:
        return ORIGINAL_PROJECT(term, registry, model, max_entries=max_entries)
    except ValueError as exc:
        if model not in MODELS or str(exc) != 'No nonzero real-component Wilson terms':
            raise
        # IMPORTANT: do not assume that a missing CG, an invalid source,
        # or cancellation of nonzero terms represents a zero contribution.
        from RGE.running.intermediate.ScalarOnlyTensorAdapters import _term_scalar_expression
        from RGE.running.intermediate.IntermediateMatcheteParsing import parse_mathematica_scalar
        source = term['operator_template_unit_flavor_inputform']
        raw = _term_scalar_expression({'Index': term.get('source_index', 0), 'TermInputForm': source})
        factor = ordered._number(parse_mathematica_scalar(raw))
        if abs(factor) > ZERO_PREF_TOL:
            raise ValueError(
                f'{model} source {term.get("source_index")}: empty Wilson term '
                f'but nonzero prefactor {factor!r}; original error retained') from exc
        basis = ordered.scalar_basis(model)
        ZERO_TERMS.append({'model': model,
                           'source_index': int(term['source_index']),
                           'prefactor_real': float(factor.real),
                           'prefactor_imag': float(factor.imag)})
        return {
            'source_index': term['source_index'],
            'basis_scalar_fields': dict(basis.dimensions),
            'real_scalar_dimension_including_higgs': basis.total_dimension,
            'flavour_matrix_axes_external_indices': term.get('flavour_matrix_axes_external_indices'),
            'ordered_not_canonical_symmetrized': True,
            'cg_registry_used': [],
            'component_count': 0,
            'components': [],
        }


def mixed_require_nonzero(model, uv_state, scalar_state, *, repo_root, mu_f=1e5):
    boundary, info = ORIGINAL_MIXED(model, uv_state, scalar_state,
                                    repo_root=repo_root, mu_f=mu_f)
    norm = float(np.linalg.norm(boundary.raw))
    if not math.isfinite(norm) or norm <= 0.0:
        raise ValueError(f'{model}: entire mixed F-tree Wilson tensor vanishes; '
                         'cannot form reliable running comparison')
    return boundary, info


def load_json(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def write_json(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--old-running', type=Path, required=True)
    ap.add_argument('--tables', type=Path, required=True)
    ap.add_argument('--reference-results', nargs='+', type=Path, required=True)
    ap.add_argument('--reference', choices=('D_-2', 'A_-4'), default='D_-2')
    ap.add_argument('--output', type=Path, required=True,
                    help='New output directory below output/parallel_t3')
    ap.add_argument('--debug', action='store_true')
    args = ap.parse_args(argv)
    old = (ROOT / args.old_running).resolve()
    out = (ROOT / args.output).resolve()
    recovery = out / '_recomputed_BC_UV'
    previous_recovery = recovery / 'results.json'
    if out.exists() and not previous_recovery.is_file():
        ap.error('Output already exists without a completed recovery results.json; choose a NEW directory')
    if not (old / 'results.json').is_file():
        ap.error('Missing original running-comparison results.json')
    if not (ROOT / args.tables / 'provenance.json').is_file():
        ap.error('Missing completed tables/provenance.json')
    original = load_json(old / 'results.json')
    prior = original['rows']
    provenance = load_json(ROOT / args.tables / 'provenance.json')
    if args.reference not in provenance.get('reference_benchmarks', {}):
        ap.error('Reference absent from baseline table provenance')
    solver = provenance['wilson_solver']
    steps = int(provenance['wilson_steps'])
    recovered = []
    ZERO_TERMS.clear()
    # Reuse already recovered numerical runs after a plotting-only crash.
    # New runs bypass the original runner's optional diagnostic plot, which
    # assumes exactly three frozen stages even when only UV is selected.
    if previous_recovery.is_file():
        print(f'Reusing existing recovered counterfactuals: {previous_recovery}', flush=True)
    else:
        with patch.object(ordered, 'project_ordered_term', project_allow_structural_zero), \
             patch.object(mixed, 'mixed_exact_boundary', mixed_require_nonzero), \
             patch.object(running, 'all_models', lambda: MODELS), \
             patch.object(running, 'STAGES', ('UV',)), \
             patch.object(running, '_plot', lambda *a, **kw: None):
            opts = ['--tables', str(args.tables), '--reference-results',
                    *(str(x) for x in args.reference_results), '--reference', args.reference,
                    '--wilson-solver', solver, '--wilson-steps', str(steps),
                    '--output', str(recovery)]
            if args.debug:
                opts.append('--debug')
            running.main(opts)
    recovery_doc = load_json(recovery / 'results.json')
    recovered = recovery_doc['rows']
    keys = {(r['reference'], r['model'], r['stage_frozen']) for r in recovered}
    if keys != {(args.reference, m, 'UV') for m in MODELS}:
        raise RuntimeError(f'Unexpected recovered counterfactual inventory: {keys}')
    old_by_key = {(r['reference'], r['model'], r['stage_frozen']): r for r in prior}
    for r in recovered:
        if r.get('status') == 'VALID_EXPLORATORY':
            old_by_key[(r['reference'], r['model'], r['stage_frozen'])] = r
    merged = list(old_by_key.values())
    write_json(out / 'results.json', {'rows': merged, 'definition': original['definition']})
    col = ['reference', 'model', 'stage_frozen', 'status', 'chi2', 'full_chi2'] + [
        f'{key}_{part}' for key in running.KEYS for part in ('full', 'frozen', 'pct')] + ['error']
    with (out / 'running_comparisons.csv').open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=col, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(merged)
    valid = [r for r in recovered if r.get('status') == 'VALID_EXPLORATORY']
    write_json(out / 'recovery_manifest.json', {
        'scope': 'B_-1 and C_-1 UV-frozen only',
        'old_running': str(old),
        'reference': args.reference,
        'zero_prefactor_terms_retained': ZERO_TERMS,
        'recovered_valid': len(valid),
        'recovered_requested': len(MODELS),
        'all_others_from_original_results': True,
        'physical_matching_certified': False,
        'method': 'Structurally zero individual mixing-spurion terms are kept as empty component sets; complete F boundary must remain nonzero',
    })
    print(f'Recovered {len(valid)}/{len(MODELS)} UV-frozen B/C calculations')
    print(f'Merged CSV: {out / "running_comparisons.csv"}')
    if len(valid) != len(MODELS):
        for r in recovered:
            if r.get('status') != 'VALID_EXPLORATORY':
                print('BLOCKED', r['model'], r.get('error'))
    return 0 if len(valid) == len(MODELS) else 1


if __name__ == '__main__':
    raise SystemExit(main())
