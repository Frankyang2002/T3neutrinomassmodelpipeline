"""Read-only audit of exact duplicate TeX documents in Reports/output/full.

Run from the repository root:
    python -m tools.AuditDuplicateTex
    python -m tools.AuditDuplicateTex --root Reports/output/full --manifest duplicate_tex_audit.csv

The script never modifies the scanned files. Hashes compare raw bytes, not
mathematical equivalence. Identical filenames or expressions are not sufficient
for classifying documents as duplicates.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
from collections import defaultdict
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def resolve_path(value: Path) -> Path:
    """Resolve relative CLI paths from the current working directory."""
    return value.expanduser().resolve()


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def audit(root: Path, manifest: Path) -> tuple[int, int, int]:
    if not root.is_dir():
        raise FileNotFoundError(f'Input directory not found: {root}')
    files = sorted(root.rglob('*.tex'))
    groups: dict[str, list[Path]] = defaultdict(list)
    for path in files:
        groups[digest(path)].append(path)
    duplicate_groups = [paths for paths in groups.values() if len(paths) > 1]
    rows = []
    for checksum, paths in sorted(groups.items(), key=lambda entry: entry[1][0].as_posix()):
        canonical = paths[0]
        for path in paths:
            rows.append({
                'path': str(path),
                'sha256': checksum,
                'bytes': path.stat().st_size,
                'group_size': len(paths),
                'classification': 'EXACT_DUPLICATE' if path != canonical else ('CANONICAL_DUPLICATE_GROUP' if len(paths) > 1 else 'UNIQUE'),
                'canonical_path': str(canonical),
            })
    manifest.parent.mkdir(parents=True, exist_ok=True)
    with manifest.open('w', newline='', encoding='utf-8-sig') as stream:
        writer = csv.DictWriter(stream, fieldnames=['path', 'sha256', 'bytes', 'group_size', 'classification', 'canonical_path'])
        writer.writeheader()
        writer.writerows(rows)
    redundant = sum(len(paths) - 1 for paths in duplicate_groups)
    print(f'TeX files scanned: {len(files)}')
    print(f'Exact duplicate groups: {len(duplicate_groups)}')
    print(f'Redundant copies (candidate only): {redundant}')
    print(f'Manifest: {manifest}')
    print('No files changed. Do not remove files until references and regeneration are checked.')
    return len(files), len(duplicate_groups), redundant


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=PROJECT_ROOT / 'Reports' / 'output' / 'full')
    parser.add_argument('--manifest', type=Path, default=PROJECT_ROOT / 'duplicate_tex_audit.csv')
    args = parser.parse_args(argv)
    try:
        audit(resolve_path(args.root), resolve_path(args.manifest))
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
