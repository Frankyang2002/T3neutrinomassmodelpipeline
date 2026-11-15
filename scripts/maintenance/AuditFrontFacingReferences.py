#!/usr/bin/env python3
"""Read-only reference inventory for front-facing repository consolidation.

Run from the local repository root:
    python -X utf8 scripts/maintenance/AuditFrontFacingReferences.py

Reads audits/front_facing_review.csv and scans tracked small text files, including
Python imports and documentation links, for literal path/name references.
Writes audits/front_facing_reference_report.csv. Never edits scientific inputs.
A missing literal reference is NOT proof a physics file is unused.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
from pathlib import Path
import re
import subprocess

TEXT_EXTENSIONS = {
    '.py', '.wl', '.wls', '.m', '.ps1', '.psm1', '.sh', '.bat', '.cmd',
    '.md', '.rst', '.txt', '.json', '.jsonl', '.yaml', '.yml',
    '.toml', '.ini', '.cfg', '.csv', '.tsv', '.tex', '.bib', '.html',
}
SKIP_DIRS = ('output/', 'Reports/output/', '.cleanup_quarantine/', 'audits/')
MAX_BYTES = 3_000_000


def git_files(root: Path) -> list[str]:
    result = subprocess.run(
        ['git', '-C', str(root), '-c', 'core.quotepath=false',
         'ls-files', '-z'], capture_output=True, check=True,
    )
    return [p.decode('utf-8', errors='surrogateescape')
            for p in result.stdout.split(b'\0') if p]


def aliases(paths: list[str]) -> dict[str, set[str]]:
    """Conservative literal aliases; not an import or call graph."""
    result: dict[str, set[str]] = defaultdict(set)
    for rel in paths:
        p = Path(rel)
        result[rel.casefold()].add(rel)
        result[p.name.casefold()].add(rel)
        if p.suffix.lower() in {'.py', '.wl', '.wls', '.ps1'} and len(p.stem) >= 14:
            result[p.stem.casefold()].add(rel)
    return result


def scan(root: Path, tracked: list[str], candidates: list[str]):
    names = aliases(candidates)
    expression = re.compile('|'.join(
        re.escape(s) for s in sorted(names, key=len, reverse=True)
    ), re.I)
    found: dict[str, set[str]] = defaultdict(set)
    skipped = 0
    scanned = 0
    for rel in tracked:
        if rel.startswith(SKIP_DIRS) or Path(rel).suffix.lower() not in TEXT_EXTENSIONS:
            continue
        path = root / rel
        try:
            if not path.is_file() or path.stat().st_size > MAX_BYTES:
                skipped += 1
                continue
            content = path.read_bytes()
            if b'\0' in content:
                skipped += 1
                continue
            text = content.decode('utf-8-sig', errors='replace').replace('\\', '/')
        except (OSError, UnicodeError):
            skipped += 1
            continue
        scanned += 1
        for match in expression.finditer(text):
            for target in names[match.group(0).casefold()]:
                if rel != target:
                    found[target].add(rel)
    return found, scanned, skipped


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path.cwd(),
                        help='Local Git checkout (default: current directory)')
    args = parser.parse_args()
    root = args.root.resolve()
    infile = root / 'audits/front_facing_review.csv'
    outfile = root / 'audits/front_facing_reference_report.csv'
    if not infile.is_file():
        parser.error(f'Missing input: {infile}')
    with infile.open(encoding='utf-8-sig', newline='') as file:
        candidates = list(csv.DictReader(file))
    if not candidates or len({x['path'] for x in candidates}) != len(candidates):
        parser.error('Candidate manifest missing or has duplicate paths')
    tracked = git_files(root)
    found, scanned, skipped = scan(root, tracked, [x['path'] for x in candidates])
    fields = ['path', 'proposed_action', 'inbound_literal_references',
              'sample_referrers', 'note']
    with outfile.open('w', encoding='utf-8', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=fields)
        writer.writeheader()
        for row in candidates:
            sources = sorted(found.get(row['path'], set()))
            writer.writerow({
                'path': row['path'],
                'proposed_action': row['proposed_action'],
                'inbound_literal_references': len(sources),
                'sample_referrers': '; '.join(sources[:12]),
                'note': 'Partial literal scan only; no automatic deletion approval',
            })
    print(f'Candidates: {len(candidates)}')
    print(f'Tracked files: {len(tracked)}; scanned text files: {scanned}; skipped text files: {skipped}')
    print(f'Files with observed inbound references: {sum(bool(found.get(x["path"])) for x in candidates)}')
    print(f'REPORT: {outfile}')
    print('Nothing modified except the reference report; no archival decisions made.')


if __name__ == '__main__':
    main()
