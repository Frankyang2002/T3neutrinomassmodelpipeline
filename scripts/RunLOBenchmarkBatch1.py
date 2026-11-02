"""Restricted T3 A_-4 Batch 1 diagnostic; no physics input modified."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
from Numerical.lo_certification.benchmark.batch1 import AuditError, run


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root', type=Path, default=Path('.'))
    p.add_argument('--run-id', default='source01')
    p.add_argument('--check-only', action='store_true')
    p.add_argument('--require-physical', action='store_true')
    args = p.parse_args(argv)
    try:
        summary, path = run(args.project_root, args.run_id, args.check_only)
    except (AuditError, KeyError, ValueError, ImportError, OSError, AssertionError) as exc:
        print('BATCH1 SOURCE/BASIS AUDIT BLOCKED: ' + str(exc), file=sys.stderr)
        return 2
    print(json.dumps({
        'model': summary['model'],
        'source_validation': summary['source_validation'],
        'physical_certification': summary['physical_certification'],
        'local_matches_GitHub_snapshot': summary['local_matches_inspected_github_blobs'],
        'output': str(path) if path else None,
    }, indent=2))
    if args.require_physical:
        print('Physical C_LLSS is not certified; explicit gate fails.', file=sys.stderr)
        return 2
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
