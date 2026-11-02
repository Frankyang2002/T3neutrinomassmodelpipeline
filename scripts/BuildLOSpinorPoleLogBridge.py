#!/usr/bin/env python3
"""All-16 native LLHH pole/constant/log witness, not physical matching."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import sys


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--native-run', type=Path, default=Path(
        'output/general_observable/lo_certification_spinor/native_hard_e515d82d118284ba_native01'))
    parser.add_argument('--inventory-run', type=Path, default=Path(
        'output/general_observable/lo_certification_spinor/raw_llhh_e515d82d118284ba_inventory01'))
    parser.add_argument('--run-id', default='bridge01')
    parser.add_argument('--check-only', action='store_true')
    args = parser.parse_args(argv)
    root = args.repo_root.resolve()
    sys.path.insert(0, str(root))
    from Numerical.lo_certification.spinor.boundary import MODELS, BoundaryError, checked_write, require, tag
    from Numerical.lo_certification.spinor.pole_log_bridge import audit_all, expressions_for
    try:
        require(re.fullmatch(r'[A-Za-z0-9_-]{1,40}', args.run_id) is not None, 'Unsafe run-id')
        owned = (root / 'output/general_observable/lo_certification_spinor').resolve()
        require(owned.is_dir(), 'Owned spinor directory missing')
        def located(p: Path) -> Path:
            resolved = (p if p.is_absolute() else root / p).resolve()
            require(resolved.is_relative_to(owned), 'Source path escaped owned spinor output')
            return resolved
        native, inventory = located(args.native_run), located(args.inventory_run)
        require(native != inventory, 'Native and inventory directories must be distinct')
        summary, models = audit_all(root, native, inventory)
        output = owned / f'pole_log_{summary["source_release_summary_sha256"][:16]}_{args.run_id}'
        require(not output.exists(), f'Immutable output directory already exists: {output}')
        if args.check_only:
            print(json.dumps(summary, indent=2, sort_keys=True, ensure_ascii=False))
            print('PASS: 16/16 original sources, exact native Laurent triplets; no output written')
            return 0
        # Entire batch validated before writes begin. Output never overlaps
        # existing native, inventory, tree, or any other worker directory.
        output.mkdir(parents=False, exist_ok=False)
        for model in MODELS:
            d = output / tag(model)
            checked_write(d / 'trios.json', models[model])
            for chirality in ('PL', 'PR'):
                for kind, text in expressions_for(models[model], chirality).items():
                    f = d / f'{chirality}_{kind}.inputform.txt'
                    f.write_text(text + '\n', encoding='utf8')
        checked_write(output / 'summary.json', summary)
        print('CREATED:', output.relative_to(root))
        print('VERIFIED:', summary['totals'])
        print('BLOCKED: MSbar counterterms, EFT1 soft subtraction, native spinor phase; physical C5 null')
        return 0
    except (BoundaryError, ValueError, KeyError, OSError, TypeError, AssertionError) as e:
        print(f'BLOCKED: {type(e).__name__}: {e}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
