"""Conservative, reversible T3 Reports cleanup.

Dry-run by default. Only two categories are eligible for quarantine:
  * LaTeX transient files (.aux, .log, .compile.log) under Reports/output.
  * Exact content duplicates of image/PDF files under Reports, provided an
    identical retained copy exists in a preferred location.

Never acts on raw output/, source files, JSON/CSV, .tex, HTML or thesis drafts.
This is NOT a complete report deletion tool: presentation copies may still be
required by dashboard builders and are therefore left alone unless identical.

Usage:
    python -m tools.CleanupReports
    python -m tools.CleanupReports --apply
    python -m tools.CleanupReports --restore
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REPORTS = ROOT / 'Reports'
QUARANTINE = ROOT / '.cleanup_quarantine' / 'Reports'
MANIFEST = ROOT / '.cleanup_quarantine' / 'manifest.json'
DIGEST_EXT = {'.png', '.jpg', '.jpeg', '.svg', '.pdf'}
TRANSIENT_EXT = {'.aux', '.log'}


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for data in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(data)
    return h.hexdigest()


def priority(path: Path) -> tuple[int, str]:
    rel = path.relative_to(ROOT).as_posix()
    if rel.startswith('Reports/output/full/'):
        return (0, rel)
    if rel.startswith('Reports/thesis_draft/'):
        return (1, rel)
    if rel.startswith('Reports/output/'):
        return (2, rel)
    if rel.startswith('Reports/figure_review/'):
        return (3, rel)
    return (1, rel)  # Unknown report areas should be kept.


def candidates():
    if not REPORTS.is_dir():
        raise FileNotFoundError(f'Reports directory not found: {REPORTS}')
    assets = defaultdict(list)
    transient = []
    for p in REPORTS.rglob('*'):
        if not p.is_file():
            continue
        relative = p.relative_to(ROOT).as_posix()
        if relative.startswith('Reports/output/') and (p.suffix.lower() in TRANSIENT_EXT):
            transient.append(p)
        if p.suffix.lower() in DIGEST_EXT:
            assets[(p.stat().st_size, digest(p))].append(p)
    entries = []
    for p in sorted(transient):
        entries.append({'path': p.relative_to(ROOT).as_posix(), 'reason': 'latex-transient', 'canonical': ''})
    for same in assets.values():
        if len(same) < 2:
            continue
        canonical = min(same, key=priority)
        for p in sorted(same):
            if p == canonical or p.is_relative_to(ROOT / 'Reports' / 'thesis_draft'):
                continue
            # Only remove duplicates from well-defined generated/review areas.
            rel = p.relative_to(ROOT).as_posix()
            if not (rel.startswith('Reports/output/') or rel.startswith('Reports/figure_review/')):
                continue
            entries.append({'path': rel, 'reason': 'byte-identical-duplicate',
                            'canonical': canonical.relative_to(ROOT).as_posix()})
    return sorted(entries, key=lambda x: x['path'])


def write_csv(entries, target):
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=['path', 'reason', 'canonical'])
        w.writeheader()
        w.writerows(entries)


def restore():
    if not MANIFEST.is_file():
        raise FileNotFoundError(f'No quarantine manifest found: {MANIFEST}')
    entries = json.loads(MANIFEST.read_text(encoding='utf-8'))
    for row in entries:
        source = QUARANTINE / Path(row['path']).relative_to('Reports')
        dest = ROOT / row['path']
        if not source.is_file():
            raise FileNotFoundError(source)
        if dest.exists():
            raise FileExistsError(f'Refusing to overwrite: {dest}')
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(source), str(dest))
    MANIFEST.unlink()
    print(f'Restored {len(entries)} files.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    operation = parser.add_mutually_exclusive_group()
    operation.add_argument('--apply', action='store_true', help='Move identified files into quarantine')
    operation.add_argument('--restore', action='store_true', help='Restore files from quarantine manifest')
    parser.add_argument('--manifest-csv', type=Path, default=ROOT / 'report_cleanup_candidates.csv')
    args = parser.parse_args()
    if args.restore:
        restore()
        return
    items = candidates()
    write_csv(items, args.manifest_csv)
    counts = defaultdict(int)
    for item in items:
        counts[item['reason']] += 1
    print(f'Candidates: {len(items)}; by category: {dict(counts)}')
    print('CSV:', args.manifest_csv)
    if not args.apply:
        print('DRY RUN: no files moved. Review the CSV, then run with --apply.')
        return
    if MANIFEST.exists():
        raise FileExistsError(f'Previous quarantine manifest exists: {MANIFEST}. Restore first.')
    QUARANTINE.mkdir(parents=True, exist_ok=True)
    moved = []
    try:
        for entry in items:
            src = ROOT / entry['path']
            dst = QUARANTINE / Path(entry['path']).relative_to('Reports')
            dst.parent.mkdir(parents=True, exist_ok=True)
            if dst.exists():
                raise FileExistsError(dst)
            shutil.move(str(src), str(dst))
            moved.append(entry)
    finally:
        MANIFEST.parent.mkdir(parents=True, exist_ok=True)
        MANIFEST.write_text(json.dumps(moved, indent=2) + '\n', encoding='utf-8')
    print(f'Quarantined {len(moved)} files. To undo: python -m tools.CleanupReports --restore')


if __name__ == '__main__':
    main()
