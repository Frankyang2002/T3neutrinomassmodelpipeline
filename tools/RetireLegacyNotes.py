"""Reversibly retire legacy T3 notes after consolidated documentation review.

Run as: python -m tools.RetireLegacyNotes [--apply | --restore]
Preview is the default. No repository changes unless --apply is given.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
QUARANTINE = ROOT / '.cleanup_quarantine' / 'legacy_notes'
MANIFEST = QUARANTINE / 'manifest.json'
PREVIEW = ROOT / 'legacy_notes_retirement_preview.csv'

# Only documents examined in the documentation consolidation; do not touch README.md
# or any source, test, config, thesis, or numerical result files.
LEGACY = [
    'INFO_IMPLEMENTATION.md', 'INFO_LAGRANGIAN.md', 'INFO_NUMERICAL.md',
    'INFO_PIPELINE.md', 'INFO_RGE.md', 'INFO_T3_BENCHMARK_SELECTION_AND_NUFIT_FEASIBILITY.md',
    'INFO_T3_MULTIPLICITY.md', 'INFO_T3_Numerical_Inputs_Sources_and_Justification.md',
    'INFO_T3_NUMERICAL_RESULTS.md', 'INFO_TESTS.md',
    'README.txt', 'README_GammaControls.txt', 'README_INSTALL.txt',
    'README_STAGE8.md', 'README_STAGE10.md', 'README_STAGE11.md', 'README_STAGE12.md',
    'README_StudyA_RG_controls.txt', 'README_StudyB.txt',
    'Reports/notes/T3_Numerical_Findings.md',
    'docs/archive/INFO_LLSS_DIAGNOSTIC.md',
    'docs/archive/INFO_T3_BENCHMARK_VALIDATION.md',
    'docs/archive/INFO_T3_NUMERICAL_RESUME.md',
    'docs/archive/INFO_T3_NUMERICAL_SCENARIO_FIX.md',
]
CONSOLIDATED = ['README.md', 'docs/PHYSICS.md', 'docs/NUMERICAL.md', 'docs/RESULTS.md']


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def preview() -> list[dict[str, str]]:
    rows = []
    for relative in LEGACY:
        p = ROOT / relative
        rows.append({'path': relative, 'status': 'PRESENT' if p.is_file() else 'MISSING',
                     'bytes': str(p.stat().st_size) if p.is_file() else '',
                     'sha256': sha(p) if p.is_file() else ''})
    with PREVIEW.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=['path', 'status', 'bytes', 'sha256'])
        writer.writeheader()
        writer.writerows(rows)
    print(f'Legacy notes: {len(rows)}; present: {sum(r["status"] == "PRESENT" for r in rows)}')
    print(f'Preview: {PREVIEW}')
    return rows


def apply() -> None:
    if MANIFEST.exists():
        raise SystemExit(f'Existing manifest: {MANIFEST}. Restore first; refusing to overwrite.')
    absent = [p for p in CONSOLIDATED if not (ROOT / p).is_file()]
    if absent:
        raise SystemExit(f'Missing consolidated docs: {absent}. Install and review first.')
    rows = preview()
    if any(r['status'] != 'PRESENT' for r in rows):
        raise SystemExit('Some legacy notes are already missing. Refusing partial retirement.')
    QUARANTINE.mkdir(parents=True, exist_ok=True)
    record = {'created_utc': datetime.now(timezone.utc).isoformat(), 'entries': rows}
    # Manifest before moves ensures recovery even if execution is interrupted.
    MANIFEST.write_text(json.dumps(record, indent=2), encoding='utf-8')
    moved = 0
    for row in rows:
        source = ROOT / row['path']
        target = QUARANTINE / row['path']
        if sha(source) != row['sha256']:
            raise SystemExit(f'Changed during operation: {source}. Partial moves can be restored.')
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(source), str(target))
        moved += 1
    print(f'Moved {moved} legacy notes to {QUARANTINE}')
    print('Do not git add -A until review: retired notes are now missing from the working tree.')


def restore() -> None:
    if not MANIFEST.is_file():
        raise SystemExit(f'No restore manifest at {MANIFEST}')
    entries = json.loads(MANIFEST.read_text(encoding='utf-8'))['entries']
    for row in entries:
        target = ROOT / row['path']
        source = QUARANTINE / row['path']
        if target.exists():
            raise SystemExit(f'Collision, refusing overwrite: {target}')
        if not source.is_file() or sha(source) != row['sha256']:
            raise SystemExit(f'Missing or changed quarantine file: {source}')
    for row in entries:
        source = QUARANTINE / row['path']
        target = ROOT / row['path']
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(source), str(target))
    MANIFEST.unlink()
    print(f'Restored {len(entries)} notes to their original paths.')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_mutually_exclusive_group()
    actions.add_argument('--apply', action='store_true')
    actions.add_argument('--restore', action='store_true')
    args = parser.parse_args()
    if args.restore:
        restore()
    elif args.apply:
        apply()
    else:
        preview()


if __name__ == '__main__':
    main()
