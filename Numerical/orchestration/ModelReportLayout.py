"""Presentation-only layout utilities for symbolic T3 model reports.

Raw symbolic pipeline artifacts stay in their historical directories because
running/matching readers can reference those paths. This module never moves raw
physics data and never launches Matchete/RGBeta or ODE solvers.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import re
import shutil
from typing import Iterable

SCENARIOS = ('smallY_smallL', 'smallY_largeL', 'largeY_smallL', 'largeY_largeL')
MODEL = re.compile(r'^T3_dS1_\d+_dS2_\d+_dF_\d+_alpha_[mp]\d+$')
SYMBOLIC_DIRS = ('Lagrangian', 'RGE', 'GroupFactors', 'Matching')
SYMBOLIC_SUFFIXES = {'.tex', '.pdf'}

@dataclass(frozen=True)
class Candidate:
    scenario: str
    model: str
    relative: Path
    source: Path
    destination: Path
    sha256: str


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def candidates(report_root: Path, destination_root: Path) -> list[Candidate]:
    source_root = Path(report_root) / 'full' / 'comparison'
    found: list[Candidate] = []
    for scenario in SCENARIOS:
        root = source_root / scenario
        if not root.is_dir():
            continue
        for model_dir in sorted(root.iterdir()):
            if not model_dir.is_dir() or not MODEL.fullmatch(model_dir.name):
                continue
            for dirname in SYMBOLIC_DIRS:
                folder = model_dir / dirname
                if not folder.is_dir():
                    continue
                for path in sorted(folder.rglob('*')):
                    if not path.is_file() or path.suffix.lower() not in SYMBOLIC_SUFFIXES:
                        continue
                    relative = path.relative_to(model_dir)
                    found.append(Candidate(scenario, model_dir.name, relative, path,
                                           destination_root / model_dir.name / relative, digest(path)))
    return found


def migrate(report_root: Path, *, execute: bool = False, prune: bool = False) -> dict:
    """Consolidate reports by model; safe dry run by default.

    Conflict policy: distinct bytes at the same destination abort BEFORE any
    write. Identical copies are deduplicated. --prune deletes only duplicates
    after the destination content is verified, never other report files.
    """
    if prune and not execute:
        raise ValueError('prune requires execute=True')
    root = Path(report_root).resolve()
    destination_root = root / 'models'
    sources = candidates(root, destination_root)
    by_dest: dict[Path, list[Candidate]] = {}
    for row in sources:
        by_dest.setdefault(row.destination, []).append(row)
    conflicts = []
    for dest, rows in by_dest.items():
        hashes = {row.sha256 for row in rows}
        if dest.exists():
            hashes.add(digest(dest))
        if len(hashes) > 1:
            conflicts.append({'destination': str(dest), 'sources': [str(r.source) for r in rows],
                              'distinct_hashes': sorted(hashes)})
    result = {'status': 'Conflict' if conflicts else 'Ready', 'mode': 'execute' if execute else 'dry-run',
              'report_root': str(root), 'files_found': len(sources), 'unique_files': len(by_dest),
              'models': sorted({r.model for r in sources}), 'conflicts': conflicts,
              'copied': 0, 'existing': 0, 'pruned': 0,
              'note': 'Only .tex/.pdf symbolic reports; raw output/ is untouched.'}
    if conflicts:
        return result
    if execute:
        for dest, rows in sorted(by_dest.items(), key=lambda kv: str(kv[0])):
            if not dest.exists():
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(rows[0].source, dest)
                if digest(dest) != rows[0].sha256:
                    raise IOError(f'Post-copy hash mismatch: {dest}')
                result['copied'] += 1
            else:
                result['existing'] += 1
        if prune:
            for row in sources:
                if digest(row.destination) != row.sha256:
                    raise IOError(f'Unexpected destination change: {row.destination}')
                row.source.unlink()
                result['pruned'] += 1
            for scenario in SCENARIOS:
                for model in result['models']:
                    basedir = root / 'full' / 'comparison' / scenario / model
                    for folder in SYMBOLIC_DIRS:
                        d = basedir / folder
                        if d.is_dir():
                            for subdir in sorted((x for x in d.rglob('*') if x.is_dir()), reverse=True):
                                if not any(subdir.iterdir()):
                                    subdir.rmdir()
                            if not any(d.iterdir()):
                                d.rmdir()
        result['status'] = 'Complete'
        destination_root.mkdir(parents=True, exist_ok=True)
        manifest = destination_root / 'report_migration_manifest.json'
        manifest.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    return result
