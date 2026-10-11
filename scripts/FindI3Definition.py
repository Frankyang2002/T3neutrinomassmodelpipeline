"""Locate the matching implementation of the I3/i3 master integral (read only).

Searches the local repository and explicitly supplied Mathematica paclet roots.
Do not assign a standard textbook normalisation to an unverified convention.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import re

PATTERN = re.compile(r'(?i)(?<![a-z0-9_])i_?3\b|(?i:masterintegral)', re.M)
SUFFIXES = {'.wl', '.m', '.wls', '.nb', '.py', '.json', '.txt', '.md'}
EXCLUDED = {'.git', '.venv', 'node_modules', 'output', 'Reports', 'report_cleanup_archive'}


def search_root(root: Path):
    for path in sorted(root.rglob('*')):
        if not path.is_file() or path.suffix.lower() not in SUFFIXES:
            continue
        if any(p in EXCLUDED for p in path.relative_to(root).parts):
            continue
        if path.stat().st_size > 10_000_000:
            continue
        try:
            lines = path.read_text(encoding='utf-8', errors='replace').splitlines()
        except OSError:
            continue
        for n, line in enumerate(lines, 1):
            if PATTERN.search(line):
                yield path, n, line.strip()[:240]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('roots', nargs='*', type=Path, default=[Path('.')])
    args = parser.parse_args()
    total = 0
    for root in args.roots:
        if not root.exists():
            print(f'[missing] {root}')
            continue
        for path, n, line in search_root(root.resolve()):
            total += 1
            print(f'{path}:{n}: {line}')
    print(f'Found {total} matching source lines. Search for an actual integral assignment/definition, not just a display reference.')


if __name__ == '__main__':
    main()
