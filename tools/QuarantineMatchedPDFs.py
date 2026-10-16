"""Reversibly quarantine generated PDFs having a same-directory TeX source.

From the repository root:
    python -m tools.QuarantineMatchedPDFs
    python -m tools.QuarantineMatchedPDFs --apply
    python -m tools.QuarantineMatchedPDFs --restore

Only Reports/output/full is scanned. No TeX or numerical data are altered.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
import shutil

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SOURCE = PROJECT_ROOT / 'Reports' / 'output' / 'full'
QUARANTINE = PROJECT_ROOT / '.cleanup_quarantine' / 'matched_report_pdfs'
MANIFEST = QUARANTINE / 'manifest.json'
PREVIEW = PROJECT_ROOT / 'matched_pdf_cleanup_preview.csv'


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def candidates() -> list[Path]:
    if not SOURCE.is_dir():
        raise FileNotFoundError(f'Missing source directory: {SOURCE}')
    return sorted(p for p in SOURCE.rglob('*.pdf') if p.with_suffix('.tex').is_file())


def preview(paths: list[Path]) -> None:
    with PREVIEW.open('w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(('pdf', 'tex', 'bytes', 'action'))
        for p in paths:
            w.writerow((p.relative_to(PROJECT_ROOT).as_posix(),
                        p.with_suffix('.tex').relative_to(PROJECT_ROOT).as_posix(),
                        p.stat().st_size, 'QUARANTINE_CANDIDATE'))
    print(f'{len(paths)} matched PDF/TeX pairs; {sum(p.stat().st_size for p in paths):,} PDF bytes')
    print(f'Preview: {PREVIEW}')


def apply(paths: list[Path]) -> None:
    if MANIFEST.exists():
        raise RuntimeError(f'Existing quarantine manifest: {MANIFEST}. Restore before a new apply.')
    if QUARANTINE.exists() and any(QUARANTINE.rglob('*.pdf')):
        raise RuntimeError(f'Quarantine already contains PDFs: {QUARANTINE}')
    entries = []
    QUARANTINE.mkdir(parents=True, exist_ok=True)
    try:
        for source in paths:
            relative = source.relative_to(SOURCE)
            destination = QUARANTINE / 'files' / relative
            if destination.exists():
                raise FileExistsError(destination)
            info = {'relative_pdf':relative.as_posix(), 'size':source.stat().st_size,
                    'sha256': digest(source)}
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(source), str(destination))
            entries.append(info)
    except Exception:
        for entry in reversed(entries):
            rel = Path(entry['relative_pdf'])
            from_path = QUARANTINE / 'files' / rel
            to_path = SOURCE / rel
            if from_path.is_file() and not to_path.exists():
                to_path.parent.mkdir(parents=True, exist_ok=True)
                shutil.move(str(from_path), str(to_path))
        raise
    document = {'created_utc':datetime.now(timezone.utc).isoformat(),
                'source': str(SOURCE), 'entries':entries}
    MANIFEST.write_text(json.dumps(document, indent=2) + '\n', encoding='utf-8')
    print(f'Quarantined {len(entries)} PDFs. Restore using --restore.')
    print(f'Manifest: {MANIFEST}')


def restore() -> None:
    if not MANIFEST.is_file():
        raise FileNotFoundError(f'No quarantine manifest: {MANIFEST}')
    data = json.loads(MANIFEST.read_text(encoding='utf-8'))
    if Path(data['source']).resolve() != SOURCE.resolve():
        raise RuntimeError('Manifest source location differs from this repository.')
    entries = data['entries']
    # Verify before moving any data; refuse overwrites and corrupted quarantine files.
    for item in entries:
        rel = Path(item['relative_pdf'])
        if rel.is_absolute() or '..' in rel.parts:
            raise ValueError('Unsafe relative path in manifest')
        src, dst = QUARANTINE / 'files' / rel, SOURCE / rel
        if dst.exists():
            raise FileExistsError(f'Restore target exists: {dst}')
        if not src.is_file() or src.stat().st_size != item['size'] or digest(src) != item['sha256']:
            raise RuntimeError(f'Quarantine file missing or modified: {src}')
    for item in entries:
        rel = Path(item['relative_pdf'])
        src, dst = QUARANTINE / 'files' / rel, SOURCE / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(src), str(dst))
    MANIFEST.unlink()
    print(f'Restored {len(entries)} PDFs to their original locations.')


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group()
    action.add_argument('--apply', action='store_true', help='Move eligible PDFs to quarantine')
    action.add_argument('--restore', action='store_true', help='Restore the prior quarantine')
    args = parser.parse_args()
    if args.restore:
        restore()
        return 0
    paths = candidates()
    preview(paths)
    if args.apply:
        apply(paths)
    else:
        print('Dry run only. Use --apply after reviewing the preview CSV.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
