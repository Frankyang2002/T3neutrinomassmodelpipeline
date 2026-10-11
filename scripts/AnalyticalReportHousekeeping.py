"""Automatically archive temporary or superseded combined-report artifacts.

Nothing in ``output/`` (raw physics output) is touched. Compilable report
sources and finished PDFs are preserved. Archive is outside ``Reports/output``.
"""
from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
import json
import shutil
import uuid

BUILD_SUFFIXES = (
    '.aux', '.log', '.out', '.toc', '.fls', '.fdb_latexmk',
    '.synctex.gz', '.compile.log', '.bbl', '.blg', '.nav', '.snm',
)
BACKUP_SUFFIXES = ('.tex.bak', '.tex.before_math_definitions.bak', '.bak')
# Optional generated presentation files; actual Matchete export lives in output/.
SUPERSEDED_NAMES = {
    'C5_original_matchete.tex', 'C5_original_matchete.pdf',
    'C5_original.tex', 'C5_original.pdf', 'C5_common_kernel.txt',
}


def classify(path: Path, *, include_superseded: bool = True) -> str | None:
    name = path.name.lower()
    if name.endswith(BUILD_SUFFIXES):
        return 'latex_build'
    if name.endswith(BACKUP_SUFFIXES) or '.before_math_definitions.' in name:
        return 'backup'
    if include_superseded and path.name in SUPERSEDED_NAMES:
        return 'superseded_presentation'
    return None


def clean(root: Path, *, apply: bool = True, include_superseded: bool = True) -> dict:
    root = Path(root).resolve()
    if root.name != 'analytical' or root.parent.name != 'output' or root.parent.parent.name != 'Reports':
        raise ValueError('Only Reports/output/analytical is allowed')
    if not root.exists():
        return {'archived': 0, 'archive': None}
    archive_root = root.parents[2] / 'report_cleanup_archive'
    files = [(p, kind) for p in sorted(root.rglob('*')) if p.is_file()
             if (kind := classify(p, include_superseded=include_superseded))]
    if not files or not apply:
        return {'archived': 0 if apply else len(files), 'archive': None, 'preview': not apply}
    batch = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '-' + uuid.uuid4().hex[:8]
    archive = archive_root / batch
    for path, _ in files:
        target = archive / path.relative_to(root)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(path), str(target))
    manifest = {'source': str(root), 'archive': str(archive), 'files': [
        {'path': str(p.relative_to(root)), 'category': kind} for p, kind in files]}
    (archive / 'cleanup_manifest.json').write_text(json.dumps(manifest, indent=2)+'\n', encoding='utf-8')
    return {'archived': len(files), 'archive': str(archive)}
