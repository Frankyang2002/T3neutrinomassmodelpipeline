#!/usr/bin/env python3
"""Extract *raw* all-16 Matchete LLHH sectors; never infer physical C5.

Project destination: scripts/ExtractLOSpinorRawLLHH.py
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
import re
import sys


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--repo-root', type=Path, default=Path(__file__).resolve().parents[1])
    ap.add_argument('--native-run', type=Path,
                    default=Path('output/general_observable/lo_certification_spinor/native_hard_e515d82d118284ba_native01'))
    ap.add_argument('--run-id', default='inventory01')
    ap.add_argument('--check-only', action='store_true')
    args = ap.parse_args(argv)
    root = args.repo_root.resolve()
    sys.path.insert(0, str(root))
    from Numerical.lo_certification.spinor.boundary import (MODELS, BoundaryError, checked_write, require, tag)
    from Numerical.lo_certification.spinor.raw_llhh import (
        FORMAT, create_inventory, sha_bytes, validate_batch_directory)
    try:
        owned = (root / 'output/general_observable/lo_certification_spinor').resolve()
        source = (root / args.native_run).resolve() if not args.native_run.is_absolute() else args.native_run.resolve()
        require(source.is_relative_to(owned), 'Input must be an existing owned native spinor run')
        require(re.fullmatch(r'[A-Za-z0-9_-]{1,40}', args.run_id) is not None,
                'Invalid run ID')
        summary, inputs = validate_batch_directory(source)
        inventories = {m: create_inventory(m, inputs[m]) for m in MODELS}
        require(len(inventories) == 16, 'Missing model inventory')
        # All sources and counts validated before any write.
        output = owned / f'raw_llhh_{summary["source_release_summary_sha256"][:16]}_{args.run_id}'
        totals = Counter()
        models = []
        for m in MODELS:
            inv = inventories[m]
            totals['candidates'] += inv['selected_llhh_term_count']
            totals['PL'] += inv['PL_holomorphic_terms']
            totals['PR'] += inv['PR_conjugate_terms_excluded_from_PL']
            totals['pole'] += inv['explicit_epsilon_pole_terms']
            totals['nonpole'] += inv['nonpole_terms_not_yet_MSbar_subtracted']
            models.append({k: v for k, v in inv.items() if k != 'selected_terms'})
        final = {
            'format': FORMAT, 'status': 'ALL16_RAW_LLHH_INVENTORIED_NOT_PHYSICAL',
            'source_native_run': str(source.relative_to(root)),
            'source_native_summary_sha256': sha_bytes((source / 'summary.json').read_bytes()),
            'mu_F_gev': summary['mu_F_gev'],
            'integrated_fermion_multiplets_together': True,
            'both_scalar_multiplets_active_at_mu_F': True,
            'SM_fermions_retained': True,
            'scalar_threshold_mu_S_not_performed': True,
            'model_count': len(models),
            'totals': dict(totals),
            'model_records': models,
            'native_gamma_cc_to_master_weyl_factor': None,
            'MSbar_subtraction_certified': False,
            'genuine_one_loop_hard_coefficient_certified': False,
            'physical_C5': None,
            'ready_for_physical_EFT1_or_scalar_matching': False,
            'required_next_physics': [
                'Identify and remove tree/no-hbar contamination from native one-loop export',
                'Perform consistent pole/counterterm and EFT1 soft subtraction in selected renormalization scheme',
                'Independently fix native GammaCC to chosen master Weyl field and epsilon convention',
                'Complete native connected four-point Dirac/Majorana exchange plus CG certification',
            ],
        }
        if args.check_only:
            print(json.dumps(final, indent=2, sort_keys=True))
            return 0
        require(not output.exists(), f'Immutable inventory directory exists: {output}')
        output.mkdir(parents=True, exist_ok=False)
        for m in MODELS:
            dst = output / tag(m)
            checked_write(dst / 'terms.json', inventories[m])
            for chirality in ('PL', 'PR'):
                terms = [v['inputform'] for v in inventories[m]['selected_terms']
                         if v['chirality'] == chirality]
                dest = dst / f'{chirality}_raw_llhh.inputform.txt'
                with dest.open('x', encoding='utf8') as file:
                    file.write(' '.join(terms) if terms else '0')
                    file.write('\n')
        checked_write(output / 'summary.json', final)
        print(json.dumps({**final, 'output': str(output.relative_to(root))}, indent=2, sort_keys=True))
        return 0
    except (OSError, ValueError, KeyError, UnicodeError) as exc:
        print(f'BLOCKED: {type(exc).__name__}: {exc}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
