"""Preview-first migration of existing full/ PNG figures into display/full_studies/.

Moves only PNG presentation files, preserving study/scenario/model hierarchy.
Never touches JSON, CSV, TXT, benchmarks or symbolic artifacts.
Run: python -m tools.MoveStudyFigures [--apply]
"""
from __future__ import annotations
import argparse
import hashlib
from pathlib import Path
import shutil
from tools.RedirectStudyFigures import FULL, DISPLAY, display_path


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--apply', action='store_true')
    args = p.parse_args(argv)
    if not FULL.is_dir():
        p.error(f'Missing numerical reports directory: {FULL}')
    files = sorted(FULL.rglob('*.png'))
    plan = [(src, display_path(src)) for src in files]
    conflicts = [(s,d) for s,d in plan if d.exists() and sha(s)!=sha(d)]
    originals = [(s,d) for s,d in plan if not d.exists()]
    duplicates = [(s,d) for s,d in plan if d.exists() and sha(s)==sha(d)]
    print(f'Full study PNGs: {len(plan)}; new destinations: {len(originals)}; identical already in display: {len(duplicates)}; conflicts: {len(conflicts)}')
    for s,d in plan:
        print(f'{"CONFLICT" if (s,d) in conflicts else "MOVE" if (s,d) in originals else "DROP DUPLICATE"}: {s.relative_to(FULL)} -> {d.relative_to(DISPLAY)}')
    if conflicts:
        print('Refusing all moves because some display PNGs differ; manually review before retrying.')
        return 2
    if not args.apply:
        print('PREVIEW only. Use --apply to move PNGs; no numerical files will be touched.')
        return 0
    for src,dst in originals:
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(src), str(dst))
    for src,_ in duplicates:
        src.unlink()
    # Prune only genuinely empty directories, including old figures folders.
    for folder in sorted(FULL.rglob('*'), key=lambda x: len(x.parts), reverse=True):
        if folder.is_dir():
            try:
                folder.rmdir()
            except OSError:
                pass
    print(f'Moved {len(originals)} PNGs; removed {len(duplicates)} exact duplicate PNGs; all non-PNG files retained.')
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
