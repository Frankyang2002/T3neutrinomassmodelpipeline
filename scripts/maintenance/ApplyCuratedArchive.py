"""Apply only the explicitly approved subset of a T3 repository cleanup plan.

This utility is local-only. Default invocation checks and previews, without
moving any files. --apply archives files outside the checkout using the
existing ExtremeRepositoryCleanup SHA-256 verified archive/restore functions.

Install alongside scripts/maintenance/ExtremeRepositoryCleanup.py.
Intended plan: output/cleanup_audit/curated_archive_plan.json.

The plan is a path/hash allow-list, not proof of redundant scientific content.
The remaining technical notes and all executable physics code are untouched.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path, PurePosixPath
import sys

import ExtremeRepositoryCleanup as core

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_PLAN = Path('output/cleanup_audit/curated_archive_plan.json')
SCHEMA = 't3-curated-archive-plan-v1'
CONFIRMATION = 'ARCHIVE_CURATED_ROOT_ONLY'


def validate_relative_path(value: str) -> None:
    p = PurePosixPath(value)
    if (not value or '\\' in value or p.is_absolute() or
            any(part in {'', '.', '..'} for part in p.parts) or
            p.as_posix() != value):
        raise ValueError(f'Invalid project-relative path: {value!r}')


def prepare(root: Path, selection_path: Path) -> dict:
    """Re-evaluate the live reference scan; reject stale, dirty or held files."""
    root = root.resolve()
    if not core.is_git_checkout_root(root):
        raise RuntimeError('Expected the actual T3 Git working-tree root')
    data = json.loads(selection_path.read_text(encoding='utf-8'))
    if data.get('schema') != SCHEMA:
        raise ValueError('Unrecognized curated plan schema')
    if data.get('options') != {'include_nested_notes': False, 'include_presentations': False}:
        raise ValueError('Curated plan was not built from the default safe-scope inventory')
    approval = data.get('archive_now')
    if not isinstance(approval, list) or not approval:
        raise ValueError('Curated archive approval is empty')
    head = (core.git_output(root, 'rev-parse', 'HEAD') or b'').decode('ascii', 'replace').strip()
    if head != data.get('source_git_head'):
        raise RuntimeError('Git HEAD differs from the uploaded inventory; rebuild the curated plan')

    live = core.build_plan(root, nested_notes=False, presentations=False)
    if live.get('git_head') != head:
        raise RuntimeError('Inconsistent Git HEAD during live reference scan')
    possible = {item['path']: item for item in live['archive_candidates']}
    approved = []
    seen = set()
    for item in approval:
        if not isinstance(item, dict):
            raise ValueError('Invalid approved archive item')
        relative = item.get('path', '')
        validate_relative_path(relative)
        if relative in seen:
            raise ValueError(f'Duplicate approved path: {relative}')
        seen.add(relative)
        allowed_root = {'root-export-archive', 'root-historical-document', 'root-stray-file'}
        kind = item.get('category')
        if ((kind in allowed_root and '/' in relative) or
                (kind == 'obsolete-backup' and not relative.startswith(
                    ('Lagrangian/interactions/', 'Lagrangian/model/'))) or
                kind not in allowed_root | {'obsolete-backup'}):
            raise ValueError('Only root artifacts and explicit source backups are approved')
        current = possible.get(relative)
        if current is None:
            raise RuntimeError(f'File is no longer eligible (modified, held or missing): {relative}')
        if any(current.get(key) != item.get(key) for key in ('path', 'category', 'bytes', 'sha256')):
            raise RuntimeError(f'File metadata/hash changed since the uploaded inventory: {relative}')
        full = (root / relative).resolve()
        if (not full.is_relative_to(root) or not full.is_file() or
                (root / relative).is_symlink()):
            raise RuntimeError(f'Unsafe or missing approved path: {relative}')
        if core.sha256(full) != item['sha256']:
            raise RuntimeError(f'File contents changed: {relative}')
        approved.append(current)

    return {'archive_candidates': approved, 'git_head': head,
            'counts': {'approved': len(approved), 'bytes': sum(i['bytes'] for i in approved),
                       'deferred': len(data.get('deferred_technical_notes', [])),
                       'held': len(data.get('held_unchanged', [])),
                       'current_original_eligible': len(live['archive_candidates'])}}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT)
    parser.add_argument('--plan', type=Path, default=DEFAULT_PLAN)
    parser.add_argument('--apply', action='store_true')
    parser.add_argument('--confirm', default='')
    parser.add_argument('--restore', type=Path, metavar='EXTERNAL_ARCHIVE')
    args = parser.parse_args(argv)
    root = args.root.resolve()
    if args.restore is not None:
        if args.apply:
            parser.error('--apply and --restore are mutually exclusive')
        restored = core.restore(root, args.restore)
        print(f'Restored {restored} files from {args.restore}')
        return 0
    if args.apply and args.confirm != CONFIRMATION:
        parser.error(f'--apply requires --confirm {CONFIRMATION}')
    selection = args.plan if args.plan.is_absolute() else root / args.plan
    plan = prepare(root, selection)
    print('Git HEAD:', plan['git_head'])
    print('Approved archive files:', plan['counts']['approved'])
    print('Approved bytes:', plan['counts']['bytes'])
    print('Scientific notes deferred:', plan['counts']['deferred'])
    print('Already reference-held files:', plan['counts']['held'])
    print('Current total eligible according to original utility:',
          plan['counts']['current_original_eligible'])
    if not args.apply:
        print('Preview only. No files changed.')
        return 0
    archive = core.apply_plan(root, plan)
    print('SHA-256 verified external archive:', archive)
    print('No code, scientific notes, configs, input/output data or GitHub content changed.')
    print('Verify tests and git status. Restore using --restore <archive-path>.')
    return 0


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (OSError, ValueError, RuntimeError) as exc:
        print('ERROR:', exc, file=sys.stderr)
        raise SystemExit(1) from exc
