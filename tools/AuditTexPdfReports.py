"""Read-only TeX/PDF report audit. Does not alter any project files.

From repository root:
    python -m tools.AuditTexPdfReports
    python -m tools.AuditTexPdfReports --root Reports/output/full/comparison
"""
from __future__ import annotations
import argparse
import csv
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
EXTENSIONS = {'.tex', '.pdf', '.aux', '.log', '.out', '.toc', '.fls', '.fdb_latexmk', '.synctex.gz', '.bbl', '.blg'}
BUILD_SUFFIXES = {'.aux', '.log', '.out', '.toc', '.fls', '.fdb_latexmk', '.synctex.gz', '.bbl', '.blg'}

def suffix(path: Path) -> str:
    name = path.name.lower()
    return '.synctex.gz' if name.endswith('.synctex.gz') else path.suffix.lower()

def classification(path: Path) -> tuple[str, str]:
    ext = suffix(path)
    if ext in BUILD_SUFFIXES:
        return 'REVIEW_BUILD_BYPRODUCT', 'Standard LaTeX build artifact; verify no diagnostic value before removal'
    if ext == '.tex':
        return 'PRESERVE_SOURCE_TEX', 'Potential analytical expression, report source, or generated source: not automatically removable'
    if ext == '.pdf':
        return 'REVIEW_PDF', 'May be compiled from TeX, but could be a unique scientific/reference document'
    return 'OTHER', 'Not a TeX/PDF file'

def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, default=Path('Reports/output/full'))
    p.add_argument('--manifest', type=Path, default=Path('tex_pdf_report_audit.csv'))
    args = p.parse_args(argv)
    source = (ROOT / args.root).resolve() if not args.root.is_absolute() else args.root.resolve()
    if not source.is_dir():
        p.error(f'Input directory does not exist: {source}')
    rows = []
    for item in sorted(source.rglob('*')):
        if not item.is_file() or suffix(item) not in EXTENSIONS:
            continue
        status, explanation = classification(item)
        sibling = item.with_suffix('.tex') if suffix(item) == '.pdf' else None
        rows.append({
            'path': item.relative_to(ROOT).as_posix() if item.is_relative_to(ROOT) else str(item),
            'extension': suffix(item), 'classification': status, 'bytes': item.stat().st_size,
            'same_stem_tex_exists': bool(sibling and sibling.is_file()),
            'reason': explanation,
        })
    destination = (ROOT / args.manifest).resolve() if not args.manifest.is_absolute() else args.manifest.resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open('w', newline='', encoding='utf-8') as output:
        writer = csv.DictWriter(output, fieldnames=['path','extension','classification','bytes','same_stem_tex_exists','reason'])
        writer.writeheader(); writer.writerows(rows)
    counts = Counter(r['classification'] for r in rows)
    print('Scanned:', source)
    print('Manifest:', destination)
    for kind, count in sorted(counts.items()):
        size = sum(r['bytes'] for r in rows if r['classification'] == kind)
        print(f'{kind}: {count} files, {size / 2**20:.2f} MiB')
    print('No files moved or deleted.')
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
