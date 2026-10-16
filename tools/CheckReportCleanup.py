"""Read-only checks after TeX/PDF report quarantine.

Run: python -m tools.CheckReportCleanup
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE_EXTS = {'.py', '.ps1', '.md', '.txt', '.json', '.yaml', '.yml', '.toml'}
SKIP_DIRS = {'.git', '.venv', '__pycache__', '.cleanup_quarantine', 'output'}


def tracked(root: Path) -> list[str]:
    p = subprocess.run(['git', '-C', str(root), 'ls-files', '-z'], capture_output=True)
    if p.returncode != 0:
        raise RuntimeError('git ls-files failed: ' + p.stderr.decode(errors='replace'))
    return [x.decode('utf-8', errors='replace') for x in p.stdout.split(b'\0') if x]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root', type=Path, default=ROOT)
    ap.add_argument('--out', type=Path, default=ROOT / 'report_cleanup_check.csv')
    args = ap.parse_args()
    root = args.root.resolve()
    names = tracked(root)
    records: list[dict[str, str]] = []
    for rel in names:
        relp = Path(rel)
        if relp.parts[:2] == ('Reports', 'output'):
            path = root / relp
            category = ('PRESENTATION' if path.suffix.lower() in {'.png', '.pdf', '.html', '.svg', '.jpg', '.jpeg'}
                        else 'DOCUMENT' if path.suffix.lower() in {'.tex', '.txt', '.md'}
                        else 'DATA_OR_OTHER')
            records.append({'kind': 'tracked_report', 'path': rel, 'detail': category,
                            'exists': str(path.is_file())})
    missing_tex = {r['path'].replace('\\', '/') for r in records if r['path'].lower().endswith('.tex') and r['exists'] == 'False'}
    # Only explicit quoted references to a missing filename are flagged; dynamic paths need manual assessment.
    for p in root.rglob('*'):
        if not p.is_file() or p.suffix.lower() not in SOURCE_EXTS:
            continue
        if any(part in SKIP_DIRS or part.startswith('display_stage') for part in p.relative_to(root).parts):
            continue
        if p.stat().st_size > 1_000_000:
            continue
        try:
            body = p.read_text(encoding='utf-8-sig')
        except (OSError, UnicodeError):
            continue
        for m in re.finditer(r'[\w./\\-]+\.tex\b', body, flags=re.I):
            token = m.group(0).replace('\\', '/')
            hits = [q for q in missing_tex if q.endswith('/' + token) or q == token]
            if hits:
                records.append({'kind': 'source_reference_review', 'path': str(p.relative_to(root)),
                                'detail': token, 'exists': 'False'})
    out = args.out if args.out.is_absolute() else root / args.out
    with out.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=['kind', 'path', 'detail', 'exists'])
        writer.writeheader()
        writer.writerows(records)
    from collections import Counter
    counts = Counter((r['kind'], r['detail']) for r in records)
    print('Tracked report files:', sum(r['kind'] == 'tracked_report' for r in records))
    print('Tracked report paths currently missing:', sum(r['kind'] == 'tracked_report' and r['exists'] == 'False' for r in records))
    print('Potential explicit references to missing TeX:', sum(r['kind'] == 'source_reference_review' for r in records))
    print('Report categories:', {k[1]: v for k, v in counts.items() if k[0] == 'tracked_report'})
    print('CSV:', out)
    print('Note: this does not resolve dynamically constructed file paths.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
